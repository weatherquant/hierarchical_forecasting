##################################################
### HIERARCHICAL ANALYSIS
### AVIATION FORECAST NOISE ANALYSIS
########################################################
########################################################
########################################################
library(data.table)
library(stringr)
library(lubridate)
library(ggplot2)
# install.packages("cmdstanr", repos = c('https://stan-dev.r-universe.dev', getOption("repos")))
# library(cmdstanr)
# check_cmdstan_toolchain()

install.packages("MCMCglmm")
library(MCMCglmm)
library(coda)
library(dplyr)
library(purrr)
library(tibble)
library(lme4)

#########
# FUNCTION 1 — Parse TAF (INCLUDING TEMPO/INTER)
expand_taf <- function(file) {
  taf <- fread(file, header = FALSE, sep = "\n")
  out <- data.table()
  
  for (i in 1:nrow(taf)) {
    line <- taf$V1[i]
    
    date_str <- str_extract(line, "^\\d{8}")
    time_str <- str_extract(line, "(?<=\\s)\\d{6}(?=\\s)")
    
    issued <- ymd_hms(paste(date_str, time_str))
    if (is.na(issued))
      next
    
    # ✅ FIXED HAZARD DETECTION
    has_fg <- str_detect(line, "\\bFG\\b")
    has_ts <- str_detect(line, "\\bTS\\b|TSRA|VCTS|TSGR")
    
    vis_groups <- str_extract_all(line, "(?<!\\d)\\d{4}(?!\\d)")[[1]]
    
    vis_groups <- suppressWarnings(as.numeric(vis_groups))
    
    has_vis <- any(vis_groups <= 5000 &
                     vis_groups >= 0, na.rm = TRUE)
    
    cld_vals <- str_extract_all(line, "\\b(?:BKN|OVC)\\d{3}\\b")[[1]]
    
    cld_vals <- as.numeric(str_extract(cld_vals, "\\d{3}")) * 100
    
    has_cld1000 <- any(cld_vals < 1000, na.rm = TRUE)
    
    times <- seq(floor_date(issued, "hour"),
                 floor_date(issued + hours(12), "hour"),
                 by = "1 hour")
    
    dt <- data.table(
      Issued = issued,
      validTime = times,
      fg_fcst = has_fg,
      ts_fcst = has_ts,
      vis_fcst = has_vis,
      cld1000_fcst = has_cld1000
    )
    
    out <- rbind(out, dt)
  }
  
  return(out)
}

# FUNCTION 2 — Parse METAR → events
parse_metar <- function(file) {
  met <- fread(file, header = FALSE, sep = "\n")
  
  dt <- met[, {
    line <- V1
    
    # Extract datetime
    date_str <- str_extract(line, "^\\d{8}")
    time_str <- str_extract(line, "(?<=\\s)\\d{6}(?=\\s)")
    
    datetime <- suppressWarnings(ymd_hms(paste(date_str, time_str)))
    
    # If bad parse, return NA row (NOT return())
    tokens <- str_split(line, " ")[[1]]
    tokens <- tokens[tokens != "" & tokens != "////"]
    
    vis <- str_extract(line, "(?<=VIS:)\\d{4}")
    vis <- suppressWarnings(as.numeric(vis))
    
    wx <- paste(tokens[grep("FG|TS|RA|SH|DZ", tokens)], collapse = " ")
    
    cld <- tokens[grep("FEW|SCT|BKN|OVC", tokens)]
    cld_base <- suppressWarnings(min(as.numeric(str_sub(cld, 4, 6)) * 100, na.rm =
                                       TRUE))
    
    # always return SAME structure
    list(
      time = datetime,
      vis = vis,
      wx = wx,
      cld_base = cld_base
    )
    
  }, by = 1:nrow(met)]
  
  # Clean bad rows
  dt <- dt[!is.na(time)]
  
  # Align to hourly grid
  dt[, time := floor_date(time, "hour")]
  
  # Hazard definitions (Centralised)
  dt[, `:=`(
    event_fg = str_detect(wx, "FG"),
    event_ts = str_detect(wx, "TS"),
    event_vis_lt1000 = !is.na(vis) & vis < 1000,
    event_cld_lt1000 = !is.infinite(cld_base) & cld_base < 1000,
    event_cld_lt2500 = !is.infinite(cld_base) & cld_base < 2500
  )]
  
  # Define event
  dt[, event_any := event_fg |
       event_ts | event_vis_lt1000 | event_cld_lt1000]
  
  return(dt[event_any == TRUE])
  
}

# FUNCTION 3 — Compute lead time
compute_lead_times <- function(taf_dt, met_dt) {
  results <- data.table()
  
  # collapse TAF to issue-level
  taf_summary <- taf_dt[, .(
    fg = any(fg_fcst),
    ts = any(ts_fcst),
    vis = any(vis_fcst),
    cld1000 = any(cld1000_fcst)
  ), by = Issued]
  
  for (i in 1:nrow(met_dt)) {
    t_event <- met_dt$time[i]
    
    # identify hazard type for this METAR
    hazard <- if (met_dt$event_fg[i])
      "FG"
    else
      if (met_dt$event_ts[i])
        "TS"
    else
      if (met_dt$event_vis_lt1000[i])
        "VIS"
    else
      if (met_dt$event_cld_lt1000[i])
        "CLD"
    else
      NA
    
    if (is.na(hazard))
      next
    
    # select matching forecasts
    candidates <- taf_summary[Issued <= t_event &
                                Issued >= t_event - hours(24) &
                                (
                                  (hazard == "FG"  & fg) |
                                    (hazard == "TS"  & ts) |
                                    (hazard == "VIS" & vis) |
                                    (hazard == "CLD" & cld1000)
                                )]
    
    if (nrow(candidates) == 0)
      next
    
    first <- candidates[which.min(Issued)]
    last  <- candidates[which.max(Issued)]
    
    results <- rbind(
      results,
      data.table(
        event_time = t_event,
        hazard = hazard,
        first_issue = first$Issued,
        last_issue  = last$Issued,
        lead_first = as.numeric(difftime(t_event, first$Issued, units = "hours")),
        lead_last  = as.numeric(difftime(t_event, last$Issued, units = "hours")),
        signal_gap = as.numeric(difftime(last$Issued, first$Issued, units =
                                           "hours"))
      )
    )
  }
  
  return(results)
}

# FUNCTION 4 - Run Full Analysis
run_analysis <- function(taf_files, metar_files) {
  all_results <- data.table()
  
  for (k in seq_along(taf_files)) {
    cat("Processing:", taf_files[k], "\n")
    
    taf_dt <- expand_taf(taf_files[k])
    met_dt <- parse_metar(metar_files[k])
    
    res <- compute_lead_times(taf_dt, met_dt)
    
    if (nrow(res) == 0)
      next
    
    res[, aerodrome := toupper(gsub("TAF_|_2025\\.txt", "", taf_files[k]))]
    #res[, aerodrome := toupper(gsub("TAF_|_[0-9]{4}\\.txt", "", taf_files[k]))]
    
    all_results <- rbind(all_results, res)
  }
  
  return(all_results)
}

# FUNCTION 5 - EXPAND TAF WITH AMD
expand_taf_with_amd <- function(file) {
  taf <- fread(file, header = FALSE, sep = "\n")
  out <- data.table()
  
  for (i in 1:nrow(taf)) {
    line <- taf$V1[i]
    
    date_str <- str_extract(line, "^\\d{8}")
    time_str <- str_extract(line, "(?<=\\s)\\d{6}(?=\\s)")
    
    issued <- ymd_hms(paste(date_str, time_str))
    if (is.na(issued))
      next
    
    is_amd <- str_detect(line, "\\bAMD\\b")
    
    # ✅ FIXED TS
    has_fg <- str_detect(line, "\\bFG\\b")
    has_ts <- str_detect(line, "\\bTS\\b|TSRA|VCTS|TSGR")
    
    vis_groups <- str_extract_all(line, "(?<!\\d)\\d{4}(?!\\d)")[[1]]
    
    vis_groups <- suppressWarnings(as.numeric(vis_groups))
    
    has_vis <- any(vis_groups <= 1000 &
                     vis_groups >= 0, na.rm = TRUE)
    
    cld_vals <- str_extract_all(line, "\\b(?:BKN|OVC)\\d{3}\\b")[[1]]
    
    cld_vals <- as.numeric(str_extract(cld_vals, "\\d{3}")) * 100
    
    has_cld1000 <- any(cld_vals < 1000, na.rm = TRUE)
    
    times <- seq(floor_date(issued, "hour"),
                 floor_date(issued + hours(12), "hour"),
                 by = "1 hour")
    
    dt <- data.table(
      Issued = issued,
      validTime = times,
      is_amd = is_amd,
      fg_fcst = has_fg,
      ts_fcst = has_ts,
      vis_fcst = has_vis,
      cld1000_fcst = has_cld1000
    )
    
    out <- rbind(out, dt, fill = TRUE)
  }
  
  return(out)
}

# FUNCTION 6 - AMD-AWARE LEAD TIME CALCULATION
compute_amd_analysis <- function(taf_dt, met_dt) {
  results <- data.table()
  
  for (i in 1:nrow(met_dt)) {
    t_event <- met_dt$time[i]
    
    hazard <- if (met_dt$event_fg[i])
      "FG"
    else
      if (met_dt$event_ts[i])
        "TS"
    else
      if (met_dt$event_vis_lt1000[i])
        "VIS"
    else
      if (met_dt$event_cld_lt1000[i])
        "CLD"
    else
      NA
    
    if (is.na(hazard))
      next
    
    # CRITICAL: match forecasts at EVENT TIME
    candidates <- taf_dt[validTime >= t_event - hours(1) &
                           validTime <= t_event + hours(1) &
                           Issued <= t_event &
                           Issued >= t_event - hours(24)]
    
    candidates <- candidates[order(Issued, validTime), .SD[.N], by = Issued]
    
    if (nrow(candidates) < 2)
      next
    
    candidates <- candidates[order(Issued)]
    
    before <- candidates[.N - 1]
    after  <- candidates[.N]
    
    # forecast presence BEFORE (at event time)
    before_fcst <- (
      (hazard == "FG"  & before$fg_fcst) |
        (hazard == "TS"  & before$ts_fcst) |
        (hazard == "VIS" & before$vis_fcst) |
        (hazard == "CLD" & before$cld1000_fcst)
    )
    
    # forecast presence AFTER
    after_fcst <- (
      (hazard == "FG"  & after$fg_fcst) |
        (hazard == "TS"  & after$ts_fcst) |
        (hazard == "VIS" & after$vis_fcst) |
        (hazard == "CLD" & after$cld1000_fcst)
    )
    
    before_any <- before$fg_fcst |
      before$ts_fcst | before$vis_fcst | before$cld1000_fcst
    after_any  <- after$fg_fcst  |
      after$ts_fcst  | after$vis_fcst  | after$cld1000_fcst
    
    # timing distance to event
    time_diff_before <- abs(difftime(t_event, before$Issued, units = "hours"))
    time_diff_after  <- abs(difftime(t_event, after$Issued, units = "hours"))
    
    timing_improved <- time_diff_after < time_diff_before
    
    added_signal <- (!before_any & after_any)
    removed_signal <- (before_any & !after_any)
    
    # relaxed improvement (includes timing refinement)
    useful_update <- (added_signal |
                        timing_improved |
                        (!before_fcst & after_fcst))
    
    improved  <- (!before_fcst & after_fcst)
    degraded  <- (before_fcst & !after_fcst)
    unchanged <- (before_fcst == after_fcst)
    
    results <- rbind(
      results,
      data.table(
        event_time = t_event,
        hazard = hazard,
        taf_issue = after$Issued,
        lead_time = as.numeric(difftime(t_event, after$Issued, units = "hours")),
        is_amd = after$is_amd,
        useful_update = useful_update,
        improved = improved,
        degraded = degraded,
        unchanged = unchanged
      )
    )
  }
  
  return(results)
}

# FUNCTION 7 - RUN AMD ANALYSIS
run_amd_analysis <- function(taf_files, metar_files) {
  all_results <- data.table()
  
  for (k in seq_along(taf_files)) {
    cat("Processing AMD:", taf_files[k], "\n")
    
    taf_dt <- expand_taf_operational(taf_files[k])
    met_dt <- parse_metar(metar_files[k])
    
    res <- compute_amd_analysis(taf_dt, met_dt)
    
    if (nrow(res) == 0)
      next
    
    res[, aerodrome := toupper(gsub("TAF_|_2025\\.txt", "", taf_files[k]))]
    #res[, aerodrome := toupper(gsub("TAF_|_[0-9]{4}\\.txt", "", taf_files[k]))]
    
    all_results <- rbind(all_results, res)
  }
  
  return(all_results)
}

# FUNCTION 8 - Extract TEMPO / INTER / PROB forecasts
expand_taf_prob_tempo <- function(file) {
  taf <- fread(file, header = FALSE, sep = "\n")
  out <- data.table()
  
  for (i in 1:nrow(taf)) {
    line <- taf$V1[i]
    
    date_str <- str_extract(line, "^\\d{8}")
    time_str <- str_extract(line, "(?<=\\s)\\d{6}(?=\\s)")
    issued <- ymd_hms(paste(date_str, time_str))
    if (is.na(issued))
      next
    
    is_tempo <- str_detect(line, "\\bTEMPO\\b") |
      str_detect(line, "\\bINTER\\b")
    is_prob30 <- str_detect(line, "PROB30")
    is_prob40 <- str_detect(line, "PROB40")
    
    # ✅ FIXED TS
    has_fg <- str_detect(line, "\\bFG\\b")
    has_ts <- str_detect(line, "\\bTS\\b|TSRA|VCTS|TSGR")
    
    vis_groups <- str_extract_all(line, "(?<!\\d)\\d{4}(?!\\d)")[[1]]
    
    vis_groups <- suppressWarnings(as.numeric(vis_groups))
    
    has_vis <- any(vis_groups <= 1000 &
                     vis_groups >= 0, na.rm = TRUE)
    
    cld_vals <- str_extract_all(line, "\\b(?:BKN|OVC)\\d{3}\\b")[[1]]
    
    cld_vals <- as.numeric(str_extract(cld_vals, "\\d{3}")) * 100
    has_cld <- any(cld_vals < 1000, na.rm = TRUE)
    
    times <- seq(floor_date(issued, "hour"),
                 floor_date(issued + hours(12), "hour"),
                 by = "1 hour")
    
    dt <- data.table(
      Issued = issued,
      validTime = times,
      is_tempo = is_tempo,
      is_prob30 = is_prob30,
      is_prob40 = is_prob40,
      fg_fcst = has_fg,
      ts_fcst = has_ts,
      vis_fcst = has_vis,
      cld_fcst = has_cld
    )
    
    out <- rbind(out, dt, fill = TRUE)
  }
  
  return(out)
}

# FUNCTION 9 - When TEMP/INTER forecast exists → does hazard occur?
verify_tempo <- function(taf_dt, met_dt) {
  results <- data.table()
  
  # use ALL times
  met_lookup <- met_dt[, .(event_any = any(event_any)), by = time]
  
  for (i in 1:nrow(taf_dt)) {
    if (!taf_dt$is_tempo[i])
      next
    
    t <- taf_dt$validTime[i]
    
    obs_row <- met_lookup[time == t]
    if (nrow(obs_row) == 0)
      next
    
    fcst_any <- taf_dt$fg_fcst[i] | taf_dt$ts_fcst[i] |
      taf_dt$vis_fcst[i] | taf_dt$cld_fcst[i]
    
    obs_any <- obs_row$event_any
    
    results <- rbind(
      results,
      data.table(
        forecast = fcst_any,
        observed = obs_any,
        hit = fcst_any & obs_any,
        false_alarm = fcst_any & !obs_any,
        miss = !fcst_any & obs_any
      )
    )
  }
  
  return(results)
}

# FUNCTION 10 - Verify PROB forecasts
# If PROB30 used, does event occur ~30% of time?
verify_prob <- function(taf_dt, met_dt) {
  results <- data.table(prob = numeric(), observed = logical())
  
  # Filter only rows where PROB exists AND hazard is forecast
  taf_dt <- taf_dt[!is.na(prob_type) &
                     (fg_fcst | ts_fcst | vis_fcst | cld_fcst)]
  
  met_lookup <- met_dt[, .(event_any = any(event_any)), by = time]
  
  for (i in 1:nrow(taf_dt)) {
    t <- taf_dt$validTime[i]
    
    obs_row <- met_lookup[time == t]
    if (nrow(obs_row) == 0)
      next
    
    obs_any <- obs_row$event_any
    
    results <- rbind(results,
                     data.table(prob = taf_dt$prob_type[i], observed = obs_any))
  }
  
  return(results)
}

# FUNCTION 11 - TEMPO / INTER Run across aerodromes
run_tempo_analysis <- function(taf_files, metar_files) {
  all_results <- data.table()
  
  for (k in seq_along(taf_files)) {
    tempo_events <- extract_tempo_events(taf_files[k])
    met_dt <- parse_metar_full(metar_files[k])
    
    res <- verify_tempo_events(tempo_events, met_dt)
    
    if (nrow(res) == 0)
      next
    
    res[, aerodrome :=
          toupper(gsub("TAF_|_[0-9]{4}\\.txt", "", taf_files[k]))]
    
    all_results <- rbind(all_results, res, fill = TRUE)
  }
  
  return(all_results)
}

# FUNCTION 12 - PROB Run across aerodromes
run_prob_analysis <- function(taf_files, metar_files) {
  all_results <- data.table()
  
  for (k in seq_along(taf_files)) {
    prob_events <- extract_prob_events(taf_files[k])
    met_dt <- parse_metar_full(metar_files[k])
    
    res <- verify_prob_events(prob_events, met_dt)
    
    if (nrow(res) == 0)
      next
    
    res[, aerodrome :=
          toupper(gsub("TAF_|_[0-9]{4}\\.txt", "", taf_files[k]))]
    
    all_results <- rbind(all_results, res, fill = TRUE)
  }
  
  return(all_results)
}

# FUNCTION 13 - Create full METAR dataset
parse_metar_full <- function(file) {
  met <- fread(file, header = FALSE, sep = "\n")
  
  dt <- met[, {
    line <- V1
    
    date_str <- str_extract(line, "^\\d{8}")
    
    time4 <- str_match(line, "^\\d{8}\\s+(\\d{4})")[, 2]
    
    datetime <- as.POSIXct(paste0(
      substr(date_str, 1, 4),
      "-",
      substr(date_str, 5, 6),
      "-",
      substr(date_str, 7, 8),
      " ",
      substr(time4, 1, 2),
      ":",
      substr(time4, 3, 4),
      ":00"
    ),
    tz = "UTC")
    
    tokens <- str_split(line, " ")[[1]]
    tokens <- tokens[tokens != "" &
                       tokens != "////"]
    
    wind_match <- str_match(line, "\\b(\\d{3})(\\d{2})/(\\d{2})KT\\b")
    
    wind_dir <- suppressWarnings(as.numeric(wind_match[2]))
    
    wind_spd <- suppressWarnings(as.numeric(wind_match[3]))
    
    temp <- suppressWarnings(as.numeric(str_extract(line, "-?\\d+\\.\\d+(?=/)")))
    
    vis <- str_extract(line, "(?<=VIS:)\\d{4}")
    
    vis <- suppressWarnings(as.numeric(vis))
    
    wx <- paste(tokens[grep("FG|TS|RA|SH|DZ", tokens)], collapse = " ")
    
    cld <- tokens[grep("FEW|SCT|BKN|OVC", tokens)]
    
    cld_base <- suppressWarnings(min(as.numeric(str_sub(cld, 4, 6)) * 100, na.rm = TRUE))
    
    cloud_info <- extract_cloud_layers(line)
    
    wx_info <- extract_weather_phenomena(line)
    
    list(
      time = datetime,
      vis = vis,
      wx = wx,
      cld_base = cld_base,
      
      wind_dir = wind_dir,
      wind_spd = wind_spd,
      temp = temp,
      
      cloud_amount =
        cloud_info$cloud_amount,
      
      cloud_base =
        cloud_info$cloud_base,
      
      event_ra = wx_info$ra,
      event_shra = wx_info$shra,
      event_dz = wx_info$dz,
      event_br = wx_info$br,
      event_gr = wx_info$gr,
      event_gs = wx_info$gs
    )
    
  }, by = 1:nrow(met)]
  
  dt <- dt[!is.na(time)]
  
  dt[, time := floor_date(time, "hour")]
  
  dt[, `:=`(
    event_fg =
      str_detect(wx, "FG|BCFG|MIFG|PRFG"),
    
    event_ts =
      str_detect(wx, "\\bTS\\b|TSRA|VCTS|TSGR"),
    
    event_vis_lt1000 =
      !is.na(vis) &
      vis < 1000,
    
    event_cld_lt1000 =
      !is.infinite(cld_base) &
      cld_base < 1000
    
  )]
  
  dt[, event_any :=
       event_fg |
       event_ts |
       event_vis_lt1000 |
       event_cld_lt1000]
  
  dt
}

# FUNCTION 14 - Extract TEMPO time windows
extract_tempo_windows <- function(line) {
  matches <- str_match_all(line, "TEMPO\\s+(\\d{4})/(\\d{4})")[[1]]
  
  if (nrow(matches) == 0)
    return(NULL)
  
  data.table(start = matches[, 2], end   = matches[, 3])
}

# FUNCTION 15 - Extract PROB time windows
extract_prob_windows <- function(line) {
  matches <- str_match_all(line, "PROB(30|40)\\s+(\\d{4})/(\\d{4})")[[1]]
  
  if (nrow(matches) == 0)
    return(NULL)
  
  data.table(prob  = as.numeric(matches[, 2]),
             start = matches[, 3],
             end   = matches[, 4])
}

# FUNCTION 16 - Build time-aware TEMPO expansion
# To replace earlier TEMPO logic
expand_taf_tempo_precise <- function(file) {
  taf <- fread(file, header = FALSE, sep = "\n")
  out <- data.table()
  
  for (i in 1:nrow(taf)) {
    line <- taf$V1[i]
    
    date_str <- str_extract(line, "^\\d{8}")
    issue_time <- str_extract(line, "(?<=\\s)\\d{6}(?=\\s)")
    
    issued <- ymd_hms(paste(date_str, issue_time))
    if (is.na(issued))
      next
    
    # base hourly timeline
    times <- seq(floor_date(issued, "hour"),
                 floor_date(issued + hours(12), "hour"),
                 by = "1 hour")
    
    dt <- data.table(
      Issued = issued,
      validTime = times,
      hour = hour(times),
      # ✅ ADD THIS LINE
      tempo_active = FALSE,
      fg_fcst = FALSE,
      ts_fcst = FALSE,
      vis_fcst = FALSE,
      cld_fcst = FALSE
    )
    
    # check TEMPO window
    tw <- extract_tempo_window(line)
    
    if (!is.null(tw)) {
      # approximate conversion (simplified UTC day handling)
      start_hour <- as.numeric(substr(tw$start, 3, 4))
      end_hour   <- as.numeric(substr(tw$end, 3, 4))
      
      dt[hour >= start_hour &
           hour <= end_hour, tempo_active := TRUE]
      
      # hazards ONLY applied in TEMPO window
      dt[tempo_active == TRUE, `:=`(
        fg_fcst = str_detect(line, "\\bFG\\b"),
        ts_fcst = str_detect(line, "\\bTS\\b"),
        vis_fcst = any(as.numeric(
          str_extract_all(line, "\\b\\d{4}\\b")[[1]]
        ) < 1000, na.rm = TRUE),
        cld_fcst = any(as.numeric(str_extract_all(
          line, "\\d{3}"
        )[[1]]) * 100 < 1000, na.rm = TRUE)
      )]
    }
    
    out <- rbind(out, dt, fill = TRUE)
  }
  
  return(out)
}

# FUNCTION 17 - Build time-aware PROB expansion
expand_taf_prob_precise <- function(file) {
  taf <- fread(file, header = FALSE, sep = "\n")
  out <- data.table()
  
  for (i in 1:nrow(taf)) {
    line <- taf$V1[i]
    
    date_str <- str_extract(line, "^\\d{8}")
    issue_time <- str_extract(line, "(?<=\\s)\\d{6}(?=\\s)")
    
    issued <- ymd_hms(paste(date_str, issue_time))
    if (is.na(issued))
      next
    
    times <- seq(floor_date(issued, "hour"),
                 floor_date(issued + hours(12), "hour"),
                 by = "1 hour")
    
    dt <- data.table(
      Issued = issued,
      validTime = times,
      prob_type = NA_real_,
      fg_fcst = FALSE,
      ts_fcst = FALSE,
      vis_fcst = FALSE,
      cld_fcst = FALSE
    )
    
    pw <- extract_prob_window(line)
    
    if (!is.null(pw)) {
      start_hour <- as.numeric(substr(pw$start, 3, 4))
      end_hour   <- as.numeric(substr(pw$end, 3, 4))
      
      dt[, hour := hour(validTime)]
      
      dt[hour >= start_hour &
           hour <= end_hour, prob_type := pw$prob]
      
      # apply hazard ONLY in window
      dt[!is.na(prob_type), `:=`(
        fg_fcst = str_detect(line, "\\bFG\\b"),
        ts_fcst = str_detect(line, "\\bTS\\b"),
        vis_fcst = any(as.numeric(
          str_extract_all(line, "\\b\\d{4}\\b")[[1]]
        ) < 1000, na.rm = TRUE),
        cld_fcst = any(as.numeric(str_extract_all(
          line, "\\d{3}"
        )[[1]]) * 100 < 1000, na.rm = TRUE)
      )]
    }
    
    out <- rbind(out, dt, fill = TRUE)
  }
  
  return(out)
}

# FUNCTION 18 - Verify TEMPO precise
verify_tempo_precise <- function(taf_dt, met_dt) {
  met_lookup <- met_dt[, .(event_any = any(event_any)), by = time]
  
  results <- merge(taf_dt,
                   met_lookup,
                   by.x = "validTime",
                   by.y = "time",
                   all.x = TRUE)
  
  results[, event_any := ifelse(is.na(event_any), FALSE, event_any)]
  
  results[, forecast := fg_fcst | ts_fcst | vis_fcst | cld_fcst]
  
  results[, `:=`(
    hit = forecast & event_any,
    false_alarm = forecast & !event_any,
    miss = !forecast & event_any
  )]
  
  return(results)
}

# FUNCTION 19 - collapse TEMPO to window level
verify_tempo_events <- function(tempo_events, met_dt) {
  results <- data.table(forecast = logical(), observed = logical())
  
  for (i in 1:nrow(tempo_events)) {
    start <- tempo_events$start_time[i]
    end   <- tempo_events$end_time[i]
    
    obs <- met_dt[time >= start & time <= end]
    
    obs_any <- if (nrow(obs) == 0)
      FALSE
    else
      any(obs$event_any, na.rm = TRUE)
    
    results <- rbind(results,
                     data.table(forecast = tempo_events$forecast[i], observed = obs_any))
  }
  
  results[, `:=`(
    hit = forecast & observed,
    false_alarm = forecast & !observed,
    miss = !forecast & observed
  )]
  
  return(results)
}

# FUNCTION 20 - collapse TEMPO to window level
verify_prob_window <- function(taf_dt, met_dt) {
  results <- data.table(prob = numeric(), observed = logical())
  
  # ✅ Identify UNIQUE PROB windows
  prob_windows <- taf_dt[!is.na(prob_type)][, .(
    start_time = min(validTime),
    end_time   = max(validTime),
    prob       = unique(prob_type)[1]
  ), by = .(Issued, prob_type)]
  
  # ✅ IMPORTANT: grouping by Issued AND prob_type
  # ensures PROB30 and PROB40 handled separately
  # but still not enough if multiple windows exist (see below)
  
  for (i in 1:nrow(prob_windows)) {
    start <- prob_windows$start_time[i]
    end   <- prob_windows$end_time[i]
    
    obs <- met_dt[time >= start & time <= end]
    
    obs_any <- if (nrow(obs) == 0)
      FALSE
    else
      any(obs$event_any, na.rm = TRUE)
    
    results <- rbind(results,
                     data.table(prob = prob_windows$prob[i], observed = obs_any))
  }
  
  return(results)
}

# FUNCTION 21 - Extract individual PROB events
extract_prob_events <- function(file) {
  taf_blocks <- build_taf_blocks(file)
  
  results <- data.table()
  
  for (line in taf_blocks) {
    date_str <- str_extract(line, "^\\d{8}")
    issue_time <- str_extract(line, "(?<=\\s)\\d{6}(?=\\s)")
    
    issued <- ymd_hms(paste(date_str, issue_time))
    if (is.na(issued))
      next
    
    # ✅ extract ALL PROBs
    matches <- str_match_all(line, "PROB(30|40).*?(\\d{4})/(\\d{4})")[[1]]
    
    if (nrow(matches) == 0)
      next
    
    for (j in 1:nrow(matches)) {
      prob  <- as.numeric(matches[j, 2])
      start <- matches[j, 3]
      end   <- matches[j, 4]
      
      start_day <- as.numeric(substr(start, 1, 2))
      start_hour <- as.numeric(substr(start, 3, 4))
      
      end_day <- as.numeric(substr(end, 1, 2))
      end_hour <- as.numeric(substr(end, 3, 4))
      
      start_time <- update(issued, day = start_day, hour = start_hour)
      end_time   <- update(issued, day = end_day, hour = end_hour)
      
      results <- rbind(
        results,
        data.table(
          Issued = issued,
          start_time = start_time,
          end_time   = end_time,
          prob       = prob
        )
      )
    }
  }
  
  return(results)
}

# FUNCTION 22 - EXTRACT TEMP EVENTS
extract_tempo_events <- function(file) {
  taf <- fread(file, header = FALSE, sep = "\n")
  results <- data.table()
  
  for (i in 1:nrow(taf)) {
    line <- taf$V1[i]
    
    date_str <- str_extract(line, "^\\d{8}")
    issue_time <- str_extract(line, "(?<=\\s)\\d{6}(?=\\s)")
    
    issued <- ymd_hms(paste(date_str, issue_time))
    if (is.na(issued))
      next
    
    tws <- extract_tempo_windows(line)
    if (is.null(tws))
      next
    
    for (j in 1:nrow(tws)) {
      start_hour <- as.numeric(substr(tws$start[j], 3, 4))
      end_hour   <- as.numeric(substr(tws$end[j], 3, 4))
      
      start_time <- floor_date(issued, "day") + hours(start_hour)
      end_time   <- floor_date(issued, "day") + hours(end_hour)
      
      # ✅ FIXED TS
      has_fg <- str_detect(line, "\\bFG\\b")
      has_ts <- str_detect(line, "\\bTS\\b|TSRA|VCTS|TSGR")
      
      vis_groups <- str_extract_all(line, "(?<!\\d)\\d{4}(?!\\d)")[[1]]
      
      vis_groups <- suppressWarnings(as.numeric(vis_groups))
      
      has_vis <- any(vis_groups <= 1000 &
                       vis_groups >= 0, na.rm = TRUE)
      
      cld_vals <- str_extract_all(line, "\\b(?:BKN|OVC)\\d{3}\\b")[[1]]
      
      cld_vals <- as.numeric(str_extract(cld_vals, "\\d{3}")) * 100
      has_cld <- any(cld_vals < 1000, na.rm = TRUE)
      
      results <- rbind(
        results,
        data.table(
          Issued = issued,
          start_time = start_time,
          end_time   = end_time,
          forecast   = (has_fg | has_ts | has_vis | has_cld)
        ),
        fill = TRUE
      )
    }
  }
  
  return(results)
}

# FUNCTION 23 - Verify individual PROB events
verify_prob_events <- function(prob_events, met_dt) {
  results <- data.table(prob = numeric(), observed = logical())
  
  for (i in 1:nrow(prob_events)) {
    start <- prob_events$start_time[i]
    end   <- prob_events$end_time[i]
    
    obs <- met_dt[time >= start & time <= end]
    
    obs_any <- if (nrow(obs) == 0)
      FALSE
    else
      any(obs$event_any, na.rm = TRUE)
    
    results <- rbind(results,
                     data.table(prob = prob_events$prob[i], observed = obs_any))
  }
  
  return(results)
}

# FUNCTION 24 - BUILD TAF BLOCKS
build_taf_blocks <- function(file) {
  taf <- fread(file, header = FALSE, sep = "\n")
  
  blocks <- list()
  current_block <- ""
  
  for (i in 1:nrow(taf)) {
    line <- taf$V1[i]
    
    # detect new TAF start (date at beginning)
    if (str_detect(line, "^\\d{8}")) {
      if (nchar(current_block) > 0) {
        blocks <- append(blocks, current_block)
      }
      
      current_block <- line
      
    } else {
      # continuation line
      current_block <- paste(current_block, line)
    }
  }
  
  # append final
  if (nchar(current_block) > 0) {
    blocks <- append(blocks, current_block)
  }
  
  return(unlist(blocks))
}

# FUNCTION 25 - BUILD FULL DATASET
build_full_dataset <- function(taf_files, metar_files) {
  all_data <- data.table()
  
  for (k in seq_along(taf_files)) {
    cat("Processing:", taf_files[k], "\n")
    
    # Parse TAF (function, now with forecaster)
    taf_dt <- expand_taf_operational(taf_files[k])
    
    # Parse METAR (raw)
    met_dt_raw <- parse_metar_full(metar_files[k])
    
    # Ensure time is aligned to hour
    met_dt_raw[, time := floor_date(time, "hour")]
    
    # Collapse to one observation per hour
    met_dt <- met_dt_raw[, .(
      vis = median(vis, na.rm = TRUE),
      wind_spd = median(wind_spd, na.rm = TRUE),
      wind_dir = median(wind_dir, na.rm = TRUE),
      temp = median(temp, na.rm = TRUE),
      cloud_amount = median(cloud_amount, na.rm = TRUE),
      cloud_base = median(cloud_base, na.rm = TRUE),
      
      event_fg = any(event_fg),
      event_ts = any(event_ts),
      event_ra = any(event_ra),
      event_shra = any(event_shra),
      event_dz = any(event_dz),
      event_br = any(event_br),
      event_gr = any(event_gr),
      event_gs = any(event_gs)
      
    ), by = time]
    
    # Merge forecast + obs
    dt <- merge(taf_dt,
                met_dt,
                by.x = "validTime",
                by.y = "time",
                all.x = TRUE)
    
    print(names(dt))
    
    cat("\n")
    
    print(dt[!is.na(wind_spd_fcst)][1:10])
    
    # Remove rows without observations
    dt <- dt[!is.na(vis)]
    
    # Build error variables
    dt[, `:=`(
      fg_error =
        as.numeric(fg_fcst) -
        as.numeric(event_fg),
      
      ts_error =
        as.numeric(ts_fcst) -
        as.numeric(event_ts),
      
      wind_spd_error =
        wind_spd_fcst -
        wind_spd,
      
      wind_dir_error =
        ((wind_dir_fcst -
            wind_dir +
            180) %% 360) - 180,
      
      wind_dir_good =
        fifelse(vrb_fcst == TRUE &
                  wind_spd <= 5, TRUE, abs(((wind_dir_fcst -
                                               wind_dir +
                                               180) %% 360
                  ) - 180) <= 20),
      
      cloud_base_error =
        cloud_base_fcst -
        cloud_base,
      
      cloud_amount_error =
        cloud_amount_fcst -
        cloud_amount,
      
      ra_error =
        as.numeric(ra_fcst) -
        as.numeric(event_ra),
      
      shra_error =
        as.numeric(shra_fcst) -
        as.numeric(event_shra),
      
      dz_error =
        as.numeric(dz_fcst) -
        as.numeric(event_dz),
      
      br_error =
        as.numeric(br_fcst) -
        as.numeric(event_br),
      
      gr_error =
        as.numeric(gr_fcst) -
        as.numeric(event_gr),
      
      gs_error =
        as.numeric(gs_fcst) -
        as.numeric(event_gs)
    )]
    
    # Add aerodrome name
    dt[, aerodrome := toupper(gsub("TAF_|_2025\\.txt", "", taf_files[k]))]
    
    all_data <- rbind(all_data, dt, fill = TRUE)
  }
  
  return(all_data)
}

# FUNCTION 26 - CORRECT SIGNAL-BASED LEAD TIME
compute_operational_lead_times <- function(taf_dt, met_dt) {
  results_list <- list()
  idx <- 1
  
  for (i in 1:nrow(met_dt)) {
    t_event <- met_dt$start_time[i]
    hazard  <- met_dt$hazard[i]
    event_start <- met_dt$start_time[i]
    event_end   <- met_dt$end_time[i]
    
    # last issued TAF before event
    
    taf_before <- unique(taf_dt[Issued <= event_end &
                                  Issued >= event_start - hours(24), .(Issued)])
    
    if (nrow(taf_before) == 0)
      next
    
    setorder(taf_before, Issued)
    
    last_any <- taf_before[.N]
    
    # -----------------------
    # MATCH FORECAST SIGNAL
    # -----------------------
    
    if (hazard == "TS") {
      taf_signal <- taf_dt[ts_fcst == TRUE &
                             Issued <= event_end &
                             validTime >= event_start &
                             validTime <= event_end]
      
    } else if (hazard == "FG") {
      taf_signal <- taf_dt[fg_fcst == TRUE &
                             Issued <= event_end &
                             validTime >= event_start &
                             validTime <= event_end]
      
    } else if (hazard == "VIS") {
      taf_signal <- taf_dt[vis_fcst == TRUE &
                             Issued <= event_end &
                             validTime >= event_start &
                             validTime <= event_end]
      
    } else if (hazard == "CLD") {
      taf_signal <- taf_dt[cld1000_fcst == TRUE &
                             Issued <= event_end &
                             validTime >= event_start &
                             validTime <= event_end]
      
    }
    
    if (nrow(taf_signal) == 0) {
      results_list[[idx]] <- data.table(
        event_time = t_event,
        hazard = hazard,
        
        first_issue = as.POSIXct(NA),
        last_issue = as.POSIXct(NA),
        last_any_issue = last_any$Issued,
        
        lead_signal_first = NA_real_,
        lead_signal_last = NA_real_,
        lead_last_any_taf =
          as.numeric(difftime(t_event, last_any$Issued, units = "hours")),
        
        signal_window = NA_real_,
        no_forecast = TRUE,
        late_forecast = FALSE
      )
      
      idx <- idx + 1
      next
    }
    
    setorder(taf_signal, Issued)
    
    first_signal <- taf_signal[1]
    last_signal  <- taf_signal[.N]
    
    results_list[[idx]] <- data.table(
      event_time = t_event,
      hazard = hazard,
      
      first_issue = first_signal$Issued,
      last_issue = last_signal$Issued,
      last_any_issue = last_any$Issued,
      
      lead_signal_first =
        as.numeric(difftime(t_event, first_signal$Issued, units = "hours")),
      
      lead_signal_last =
        as.numeric(difftime(t_event, last_signal$Issued, units = "hours")),
      
      lead_last_any_taf =
        as.numeric(difftime(t_event, last_any$Issued, units = "hours")),
      
      signal_window =
        as.numeric(
          difftime(last_signal$Issued, first_signal$Issued, units = "hours")
        ),
      
      no_forecast = FALSE,
      late_forecast =
        first_signal$Issued > event_start &
        first_signal$Issued <= event_end
    )
    
    idx <- idx + 1
  }
  
  rbindlist(results_list, fill = TRUE)
}

# FUNCTION 27 - NEW RUNNER
run_operational_analysis <- function(taf_files, metar_files) {
  all_results <- data.table()
  
  for (k in seq_along(taf_files)) {
    cat("Processing:", taf_files[k], "\n")
    
    # ✅ STEP 1: expand TAF
    taf_dt <- expand_taf_operational(taf_files[k])
    
    # ✅ STEP 2: FULL METAR (no filtering)
    met_dt_full <- parse_metar_full(metar_files[k])
    
    # ✅ STEP 3: COLLAPSE INTO EVENTS
    events_dt <- collapse_events(met_dt_full)
    
    if (nrow(events_dt) == 0)
      next
    
    # ✅ STEP 4: COMPUTE LEAD TIMES USING EVENTS
    res <- compute_operational_lead_times(taf_dt, events_dt)
    
    if (nrow(res) == 0)
      next
    
    res[, aerodrome :=
          toupper(gsub("TAF_|_[0-9]{4}\\.txt", "", taf_files[k]))]
    
    all_results <- rbind(all_results, res, fill = TRUE)
  }
  
  return(all_results)
}

# FUNCTION 28 - COLLAPSE EVENTS TO ONE
collapse_events <- function(met_dt) {
  # Ensure sorted
  setorder(met_dt, time)
  
  # Create hazard-specific flags
  events_list <- list()
  
  hazards <- c("event_fg",
               "event_ts",
               "event_vis_lt1000",
               "event_cld_lt1000")
  
  names(hazards) <- c("FG", "TS", "VIS", "CLD")
  
  for (h in names(hazards)) {
    flag <- hazards[h]
    
    dt <- met_dt[, .(time, active = get(flag))]
    
    # Only keep rows with TRUE
    dt <- dt[active == TRUE]
    
    if (nrow(dt) == 0)
      next
    
    # Identify breaks in continuity (gap > 1 hour)
    dt[, gap := as.numeric(difftime(time, shift(time), units = "hours"))]
    dt[is.na(gap), gap := 0]
    
    # New event if gap > 1 hour
    dt[, new_event := (gap > 1)]
    
    dt[, event_id := cumsum(new_event) + 1]
    
    # Collapse to one row per event
    collapsed <- dt[, .(
      start_time = min(time),
      end_time   = max(time),
      duration_hours =
        as.numeric(difftime(max(time), min(time), units = "hours")) + 1
    ), by = event_id]
    
    collapsed[, hazard := h]
    
    events_list[[h]] <- collapsed
  }
  
  events <- rbindlist(events_list, fill = TRUE)
  
  return(events)
}

# FUNCTION 29
extract_taf_blocks <- function(line) {
  starts <- gregexpr(
    "FM\\d{6}|BECMG\\s+\\d{4}/\\d{4}|TEMPO\\s+\\d{4}/\\d{4}|INTER\\s+\\d{4}/\\d{4}|PROB30\\s+\\d{4}/\\d{4}|PROB40\\s+\\d{4}/\\d{4}",
    line,
    perl = TRUE
  )[[1]]
  
  if (starts[1] == -1) {
    return(data.table(type = "BASE", text = line))
  }
  
  starts <- c(1, starts)
  ends <- c(starts[-1] - 1, nchar(line))
  
  data.table(text = substring(line, starts, ends))[, type := fifelse(grepl("^FM", text), "FM", fifelse(
    grepl("^BECMG", text),
    "BECMG",
    fifelse(
      grepl("^TEMPO", text),
      "TEMPO",
      fifelse(
        grepl("^INTER", text),
        "INTER",
        fifelse(
          grepl("^PROB30", text),
          "PROB30",
          fifelse(grepl("^PROB40", text), "PROB40", "BASE")
        )
      )
    )
  ))]
}

# FUNCTION 30
extract_block_hazards <- function(txt) {
  # remove remarks
  txt <- strsplit(txt, "\\bRMK\\b")[[1]][1]
  
  ts_flag <- str_detect(txt, "\\bTS\\b|TSRA|VCTS|TSGR")
  
  fg_flag <- str_detect(txt, "FG|BCFG|MIFG|PRFG")
  
  vis_groups <- str_match_all(txt,
                              "\\b(\\d{4})\\s+(?:FG|BR|HZ|FU|DZ|RA|SHRA|TSRA|VCTS|FEW|SCT|BKN|OVC)")[[1]]
  
  vis_vals <- if (nrow(vis_groups) > 0) {
    as.numeric(vis_groups[, 2])
  } else {
    numeric(0)
  }
  
  vis_flag <- any(vis_vals < 5000, na.rm = TRUE)
  
  # -----------------------------
  # Cloud information
  # -----------------------------
  
  cloud_info <- extract_cloud_layers(txt)
  
  cld_flag <- cloud_info$low_cloud
  
  # -----------------------------
  # Wind
  # -----------------------------
  
  wind_match <- str_match(txt, "\\b(\\d{3}|VRB)(\\d{2})KT\\b")
  
  wind_dir_txt <- wind_match[2]
  
  vrb_fcst <- identical(wind_dir_txt, "VRB")
  
  wind_dir <- suppressWarnings(as.numeric(wind_dir_txt))
  
  wind_spd <- suppressWarnings(as.numeric(wind_match[3]))
  
  if (length(wind_dir) == 0)
    wind_dir <- NA_real_
  
  if (length(wind_spd) == 0)
    wind_spd <- NA_real_
  
  # -----------------------------
  # Weather phenomena
  # -----------------------------
  
  wx_info <- extract_weather_phenomena(txt)
  
  
  list(
    fg = fg_flag,
    ts = ts_flag,
    vis = vis_flag,
    cld = cld_flag,
    wind_dir = wind_dir,
    wind_spd = wind_spd,
    vrb_fcst = vrb_fcst,
    cloud_amount = cloud_info$cloud_amount,
    cloud_base = cloud_info$cloud_base,
    ra = wx_info$ra,
    shra = wx_info$shra,
    dz = wx_info$dz,
    br = wx_info$br,
    gr = wx_info$gr,
    gs = wx_info$gs
  )
}

# FUNCTION 31
decode_ddhh <- function(issued, ddhh) {
  day  <- as.numeric(substr(ddhh, 1, 2))
  hour <- as.numeric(substr(ddhh, 3, 4))
  
  yr  <- year(issued)
  mon <- month(issued)
  
  candidate <- as.POSIXct(sprintf("%04d-%02d-%02d %02d:00:00", yr, mon, day, hour), tz = "UTC")
  
  candidate
}

# FUNCTION 32
extract_taf_temperature <- function(txt, valid_start, valid_end) {
  temp_match <- str_match(txt, "RMK\\s+T\\s+([0-9\\s\\-]+)")
  
  if (is.na(temp_match[2])) {
    return(NULL)
  }
  
  temps <- suppressWarnings(as.numeric(strsplit(trimws(temp_match[2]), "\\s+")[[1]]))
  
  if (length(temps) == 0) {
    return(NULL)
  }
  
  temp_times <- valid_start + hours(seq_len(length(temps)) * 3)
  
  data.table(validTime = temp_times, temp_fcst = temps)
}

# FUNCTION 33 - TAF OPERATIONAL MASTER FUNCTION
expand_taf_operational <- function(file) {
  taf_blocks <- build_taf_blocks(file)
  
  out <- data.table()
  
  for (line in taf_blocks) {
    date_str <- str_extract(line, "^\\d{8}")
    
    issue_time <- str_extract(line, "(?<=\\s)\\d{6}(?=\\s)")
    
    issued <- suppressWarnings(ymd_hms(paste(date_str, issue_time), tz = "UTC"))
    
    if (is.na(issued))
      next
    
    # -----------------------------
    # AMENDMENT STATUS
    # -----------------------------
    
    is_amd <- str_detect(line, "\\bTAF\\s+AMD\\b")
    
    # -----------------------------
    # FORECASTER
    # -----------------------------
    forecaster <- str_match(line, "^\\d{8}\\s+\\d{6}\\s+(\\S+)")[, 2]
    
    if (is.na(forecaster))
      forecaster <- "UNKNOWN"
    
    # -----------------------------
    # VALIDITY PERIOD
    # -----------------------------
    
    valid_period <- str_match(line, "\\b(\\d{4})\\s+(\\d{4})\\b")
    
    if (any(is.na(valid_period)))
      next
    
    valid_start <- decode_ddhh(issued, valid_period[2])
    
    valid_end <- decode_ddhh(issued, valid_period[3])
    
    if (valid_end < valid_start) {
      valid_end <- valid_end %m+% months(1)
    }
    
    dt <- data.table(
      Issued = issued,
      
      validTime = seq(valid_start, valid_end, by = "1 hour"),
      
      forecaster = forecaster,
      
      is_amd = is_amd,
      
      forecast_type = "BASE",
      
      cloud_amount_fcst = NA_real_,
      cloud_base_fcst = NA_real_,
      
      ra_fcst = FALSE,
      shra_fcst = FALSE,
      dz_fcst = FALSE,
      br_fcst = FALSE,
      gr_fcst = FALSE,
      gs_fcst = FALSE,
      
      fg_fcst = FALSE,
      ts_fcst = FALSE,
      vis_fcst = FALSE,
      cld1000_fcst = FALSE,
      
      wind_dir_fcst = NA_real_,
      wind_spd_fcst = NA_real_,
      vrb_fcst = FALSE
    )
    
    # -----------------------------
    # BASE FORECAST
    # -----------------------------
    
    base_text <- strsplit(
      line,
      "FM\\d{6}|BECMG\\s+\\d{4}/\\d{4}|TEMPO\\s+\\d{4}/\\d{4}|INTER\\s+\\d{4}/\\d{4}|PROB30\\s+\\d{4}/\\d{4}|PROB40\\s+\\d{4}/\\d{4}"
    )[[1]][1]
    
    base_hz <- extract_block_hazards(base_text)
    
    dt[, `:=`(
      fg_fcst = base_hz$fg,
      ts_fcst = base_hz$ts,
      vis_fcst = base_hz$vis,
      cld1000_fcst = base_hz$cld,
      
      wind_dir_fcst = base_hz$wind_dir,
      wind_spd_fcst = base_hz$wind_spd,
      vrb_fcst = base_hz$vrb_fcst,
      
      cloud_amount_fcst = base_hz$cloud_amount,
      cloud_base_fcst = base_hz$cloud_base,
      ra_fcst = base_hz$ra,
      shra_fcst = base_hz$shra,
      dz_fcst = base_hz$dz,
      br_fcst = base_hz$br,
      gr_fcst = base_hz$gr,
      gs_fcst = base_hz$gs
    )]
    
    # ===================================================
    # FM GROUPS
    # ===================================================
    
    fm_matches <- str_match_all(line, "FM(\\d{6})")[[1]]
    
    if (nrow(fm_matches) > 0) {
      fm_pos <- gregexpr("FM\\d{6}", line)[[1]]
      
      for (j in seq_len(nrow(fm_matches))) {
        fm_code <- fm_matches[j, 2]
        
        start_time <- decode_ddhh(issued, substr(fm_code, 1, 4))
        
        text_start <- fm_pos[j]
        
        if (j < nrow(fm_matches)) {
          text_end <- fm_pos[j + 1] - 1
          
          next_start <- decode_ddhh(issued, substr(fm_matches[j + 1, 2], 1, 4))
          
          end_time <- next_start - hours(1)
          
        } else {
          text_end <- nchar(line)
          end_time <- valid_end
        }
        
        fm_text <- substr(line, text_start, text_end)
        
        hz <- extract_block_hazards(fm_text)
        
        dt[validTime >= start_time &
             validTime <= end_time, `:=`(
               forecast_type = "FM",
               
               fg_fcst = hz$fg,
               ts_fcst = hz$ts,
               vis_fcst = hz$vis,
               cld1000_fcst = hz$cld,
               
               wind_dir_fcst = hz$wind_dir,
               wind_spd_fcst = hz$wind_spd,
               vrb_fcst = hz$vrb_fcst,
               
               cloud_amount_fcst = hz$cloud_amount,
               cloud_base_fcst = hz$cloud_base,
               ra_fcst = hz$ra,
               shra_fcst = hz$shra,
               dz_fcst = hz$dz,
               br_fcst = hz$br,
               gr_fcst = hz$gr,
               gs_fcst = hz$gs
             )]
      }
    }
    
    # ===================================================
    # BECMG
    # ===================================================
    
    becmg_matches <- str_match_all(
      line,
      "BECMG\\s+(\\d{4})/(\\d{4})(.*?)($|TEMPO|INTER|PROB30|PROB40|FM|BECMG)"
    )[[1]]
    
    if (nrow(becmg_matches) > 0) {
      for (j in seq_len(nrow(becmg_matches))) {
        end_time <- decode_ddhh(issued, becmg_matches[j, 3])
        
        hz <- extract_block_hazards(becmg_matches[j, 4])
        
        dt[validTime >= end_time, `:=`(
          forecast_type = "BECMG",
          
          fg_fcst = hz$fg,
          ts_fcst = hz$ts,
          vis_fcst = hz$vis,
          cld1000_fcst = hz$cld,
          
          wind_dir_fcst = hz$wind_dir,
          wind_spd_fcst = hz$wind_spd,
          vrb_fcst = hz$vrb_fcst,
          
          cloud_amount_fcst = hz$cloud_amount,
          cloud_base_fcst = hz$cloud_base,
          ra_fcst = hz$ra,
          shra_fcst = hz$shra,
          dz_fcst = hz$dz,
          br_fcst = hz$br,
          gr_fcst = hz$gr,
          gs_fcst = hz$gs
        )]
      }
    }
    
    # ===================================================
    # TEMPO
    # ===================================================
    
    tempo_matches <- str_match_all(
      line,
      "TEMPO\\s+(\\d{4})/(\\d{4})(.*?)($|TEMPO|INTER|PROB30|PROB40|FM|BECMG)"
    )[[1]]
    
    if (nrow(tempo_matches) > 0) {
      for (j in seq_len(nrow(tempo_matches))) {
        start_time <- decode_ddhh(issued, tempo_matches[j, 2])
        
        end_time <- decode_ddhh(issued, tempo_matches[j, 3])
        
        hz <- extract_block_hazards(tempo_matches[j, 4])
        
        dt[validTime >= start_time &
             validTime <= end_time, `:=`(
               forecast_type = "TEMPO",
               
               fg_fcst = hz$fg,
               ts_fcst = hz$ts,
               vis_fcst = hz$vis,
               cld1000_fcst = hz$cld,
               
               wind_dir_fcst = hz$wind_dir,
               wind_spd_fcst = hz$wind_spd,
               vrb_fcst = hz$vrb_fcst,
               
               cloud_amount_fcst = hz$cloud_amount,
               cloud_base_fcst = hz$cloud_base,
               ra_fcst = hz$ra,
               shra_fcst = hz$shra,
               dz_fcst = hz$dz,
               br_fcst = hz$br,
               gr_fcst = hz$gr,
               gs_fcst = hz$gs
             )]
      }
    }
    
    # ===================================================
    # INTER
    # ===================================================
    
    inter_matches <- str_match_all(
      line,
      "INTER\\s+(\\d{4})/(\\d{4})(.*?)($|TEMPO|INTER|PROB30|PROB40|FM|BECMG)"
    )[[1]]
    
    if (nrow(inter_matches) > 0) {
      for (j in seq_len(nrow(inter_matches))) {
        start_time <- decode_ddhh(issued, inter_matches[j, 2])
        
        end_time <- decode_ddhh(issued, inter_matches[j, 3])
        
        hz <- extract_block_hazards(inter_matches[j, 4])
        
        dt[validTime >= start_time &
             validTime <= end_time, `:=`(
               forecast_type = "INTER",
               
               fg_fcst = hz$fg,
               ts_fcst = hz$ts,
               vis_fcst = hz$vis,
               cld1000_fcst = hz$cld,
               
               wind_dir_fcst = hz$wind_dir,
               wind_spd_fcst = hz$wind_spd,
               vrb_fcst = hz$vrb_fcst,
               
               cloud_amount_fcst = hz$cloud_amount,
               cloud_base_fcst = hz$cloud_base,
               ra_fcst = hz$ra,
               shra_fcst = hz$shra,
               dz_fcst = hz$dz,
               br_fcst = hz$br,
               gr_fcst = hz$gr,
               gs_fcst = hz$gs
             )]
      }
    }
    
    # ===================================================
    # PROB30 / PROB40
    # ===================================================
    
    prob_matches <- str_match_all(
      line,
      "PROB(30|40)\\s+(\\d{4})/(\\d{4})(.*?)($|TEMPO|INTER|PROB30|PROB40|FM|BECMG)"
    )[[1]]
    
    if (nrow(prob_matches) > 0) {
      for (j in seq_len(nrow(prob_matches))) {
        start_time <- decode_ddhh(issued, prob_matches[j, 3])
        
        end_time <- decode_ddhh(issued, prob_matches[j, 4])
        
        hz <- extract_block_hazards(prob_matches[j, 5])
        
        dt[validTime >= start_time &
             validTime <= end_time, `:=`(
               forecast_type =
                 paste0("PROB", prob_matches[j, 2]),
               
               fg_fcst = hz$fg,
               ts_fcst = hz$ts,
               vis_fcst = hz$vis,
               cld1000_fcst = hz$cld,
               
               wind_dir_fcst = hz$wind_dir,
               wind_spd_fcst = hz$wind_spd,
               vrb_fcst = hz$vrb_fcst,
               
               cloud_amount_fcst = hz$cloud_amount,
               cloud_base_fcst = hz$cloud_base,
               ra_fcst = hz$ra,
               shra_fcst = hz$shra,
               dz_fcst = hz$dz,
               br_fcst = hz$br,
               gr_fcst = hz$gr,
               gs_fcst = hz$gs
             )]
      }
    }
    
    out <- rbind(out, dt, fill = TRUE)
  }
  
  out
}

# FUNCTION 34 - EXTRACT TEMP AND VALIDATE NOISE MASTER FUNCTION
extract_taf_temperatures <- function(file) {
  taf_blocks <- build_taf_blocks(file)
  
  out <- data.table()
  
  for (line in taf_blocks) {
    date_str <- str_extract(line, "^\\d{8}")
    
    issue_time <- str_extract(line, "(?<=\\s)\\d{6}(?=\\s)")
    
    issued <- suppressWarnings(ymd_hms(paste(date_str, issue_time), tz = "UTC"))
    
    if (is.na(issued))
      next
    
    forecaster <- str_match(line, "^\\d{8}\\s+\\d{6}\\s+(\\S+)")[, 2]
    
    valid_period <- str_match(line, "\\b(\\d{4})\\s+(\\d{4})\\b")
    
    if (any(is.na(valid_period)))
      next
    
    valid_start <- decode_ddhh(issued, valid_period[2])
    
    # -------------------------
    # RMK T extraction
    # -------------------------
    
    t_match <- str_match(
      line,
      "RMK\\s+T\\s+([0-9\\-]+)\\s+([0-9\\-]+)\\s+([0-9\\-]+)\\s+([0-9\\-]+)"
    )
    
    if (all(is.na(t_match)))
      next
    
    temps <- as.numeric(t_match[2:5])
    
    temp_times <- valid_start +
      hours(c(3, 6, 9, 12))
    
    dt <- data.table(
      Issued = issued,
      forecaster = forecaster,
      validTime = temp_times,
      temp_fcst = temps
    )
    
    out <- rbind(out, dt, fill = TRUE)
  }
  
  out
}

# FUNCTION 35 - BUILD TEMP DATA MASTER FUNCTION
build_temperature_dataset <- function(taf_files, metar_files) {
  all_data <- data.table()
  
  for (k in seq_along(taf_files)) {
    taf_temp <- extract_taf_temperatures(taf_files[k])
    
    met <- parse_metar_full(metar_files[k])
    
    met <- met[, .(temp = median(temp, na.rm = TRUE)), by = time]
    
    dt <- merge(taf_temp,
                met,
                by.x = "validTime",
                by.y = "time",
                all.x = TRUE)
    
    dt[, temp_error :=
         temp_fcst - temp]
    
    dt[, aerodrome :=
         toupper(gsub("TAF_|_2025\\.txt", "", taf_files[k]))]
    
    all_data <- rbind(all_data, dt, fill = TRUE)
  }
  
  all_data
}

# FUNCTION 36 - CLOUD LAYER EXTRACTION
extract_cloud_layers <- function(txt) {
  cloud_match <- str_match_all(txt, "\\b(FEW|SCT|BKN|OVC)(\\d{3})\\b")[[1]]
  
  if (nrow(cloud_match) == 0) {
    return(list(
      cloud_amount = NA_real_,
      cloud_base = NA_real_,
      low_cloud = FALSE
    ))
  }
  
  amount_lookup <- c(
    FEW = 1,
    SCT = 2,
    BKN = 3,
    OVC = 4
  )
  
  cloud_amount <- max(amount_lookup[cloud_match[, 2]], na.rm = TRUE)
  
  cloud_base <- min(as.numeric(cloud_match[, 3]) * 100, na.rm = TRUE)
  
  low_cloud <- any(cloud_match[, 2] %in%
                     c("BKN", "OVC") &
                     (as.numeric(cloud_match[, 3]) * 100 < 1000))
  
  list(cloud_amount = cloud_amount,
       cloud_base = cloud_base,
       low_cloud = low_cloud)
}

# FUNCTION 37 - WEATHER PHENOMENA EXTRACTION
extract_weather_phenomena <- function(txt) {
  list(
    ra =
      grepl("\\bRA\\b", txt),
    
    shra =
      grepl("SHRA", txt),
    
    dz =
      grepl("\\bDZ\\b", txt),
    
    br =
      grepl("\\bBR\\b", txt),
    
    fg =
      grepl("\\bFG\\b", txt),
    
    ts =
      grepl("TS", txt),
    
    gr =
      grepl("\\bGR\\b", txt),
    
    gs =
      grepl("\\bGS\\b", txt)
    
  )
}

# FUNCTION 38 - BUILD COMPOSITE NOISE INDEX
build_composite_noise_index <- function(noise_data, temp_data) {
  temp_metric <- temp_data[, .(temp_mae =
                                 mean(abs(temp_error), na.rm = TRUE)), by = ForecasterID]
  
  wind_metric <- noise_data[, .(
    wind_spd_mae =
      mean(abs(wind_spd_error), na.rm = TRUE),
    
    wind_dir_mae =
      mean(abs(wind_dir_error), na.rm = TRUE)
    
  ), by = ForecasterID]
  
  cloud_metric <- noise_data[, .(
    cloud_base_mae =
      mean(abs(cloud_base_fcst -
                 cloud_base), na.rm = TRUE),
    
    cloud_amount_mae =
      mean(abs(cloud_amount_fcst -
                 cloud_amount), na.rm = TRUE)
    
  ), by = ForecasterID]
  
  hazard_metric <- noise_data[, .(
    fg_rate =
      mean(fg_fcst),
    
    ts_rate =
      mean(ts_fcst),
    
    vis_rate =
      mean(vis_fcst),
    
    cld_rate =
      mean(cld1000_fcst)
    
  ), by = ForecasterID]
  
  structure_metric <- structure_entropy
  
  out <- Reduce(
    function(x, y)
      merge(x, y, by = "ForecasterID", all = TRUE),
    
    list(
      temp_metric,
      
      wind_metric,
      
      cloud_metric,
      
      hazard_metric,
      
      structure_metric
      
    )
    
  )
  
  metric_cols <-
    
    names(out)[names(out) != "ForecasterID"]
  
  for (col in metric_cols) {
    out[, paste0(col, "_z") := scale(get(col))]
    
  }
  
  z_cols <-
    
    grep("_z$", names(out), value = TRUE)
  
  out[, CompositeNoiseIndex :=
        rowMeans(.SD, na.rm = TRUE), .SDcols = z_cols]
  
  out[order(-CompositeNoiseIndex)]
}


#######################################
# Load data
### READ TAF DATA
### AVS1 temp and wind
setwd(
  "C:/Users/jwest.INTERNAL/OneDrive - Bureau of Meteorology/Documents/Data/Aviation/AVS1"
)
taf_files <- list.files(pattern = "TAF")
metar_files <- list.files(pattern = "METAR")
taf_files <- sort(taf_files)
metar_files <- sort(metar_files)

# test
# test <- expand_taf_with_amd(taf_files[1])
# unique(test$forecaster)

# Run BUILD DATASET function
noise_data <- build_full_dataset(taf_files, metar_files)
temp_data <- build_temperature_dataset(taf_files, metar_files)

# Immediate check
summary(temp_data$temp_fcst)
summary(temp_data$temp)
summary(temp_data$temp_error)

# Anonymise Forecasters
set.seed(2025)
forecaster_lookup <- data.table(forecaster = sort(unique(noise_data$forecaster)))

forecaster_lookup[, ForecasterID := sample(sprintf("F%03d", seq_len(.N)))]

noise_data <- merge(noise_data,
                    forecaster_lookup,
                    by = "forecaster",
                    all.x = TRUE)
temp_data <- merge(temp_data, forecaster_lookup, by = "forecaster", all.x =
                     TRUE)
saveRDS(noise_data, file = "noise_data.rds")
saveRDS(temp_data, file = "temp_data.rds")
getwd()
##################################################
# 1. Prepare the modelling identifiers
# model_data <- copy(noise_data)
#
# model_data <- model_data[!is.na(ForecasterID) &
#                            ForecasterID != "UNKNOWN" &
#                            !is.na(aerodrome) &
#                            !is.na(Issued) &
#                            !is.na(validTime)]
#
# model_data[, `:=`(
#   lead_hours = as.numeric(difftime(validTime, Issued, units = "hours")),
#
#   taf_id = interaction(aerodrome, format(Issued, "%Y%m%d%H%M%S"), drop = TRUE),
#
#   verify_case_id = interaction(aerodrome, format(validTime, "%Y%m%d%H"), drop = TRUE),
#
#   issue_cycle = factor(floor(hour(Issued) / 3) * 3),
#
#   valid_hour = factor(hour(validTime)),
#   month = factor(month(validTime)),
#
#   ForecasterID = factor(ForecasterID),
#   aerodrome = factor(aerodrome),
#   # taf_id = factor(taf_id),
#   # verify_case_id = factor(verify_case_id),
#   forecast_type = factor(forecast_type),
#   is_amd = factor(is_amd)
# )]
#
# model_data <- model_data[lead_hours >= 0 &
#                            lead_hours <= 30]
#
# model_data[, lead_z :=
#              as.numeric(scale(lead_hours))]
#
# setorder(model_data, aerodrome, validTime, Issued)
#
# model_data[, `:=`(
#   forecast_sequence = seq_len(.N),
#   n_forecasts_for_case = .N,
#   is_latest_forecast = seq_len(.N) == .N
# ), by = verify_case_id]
#
# model_data[, sequence_z :=
#              as.numeric(scale(forecast_sequence))]

##################################################
##### FORECAST VINTAGE STRUCTURE
##################################################
analysis_data <- copy(noise_data)

analysis_data <- analysis_data[!is.na(ForecasterID) &
                                 ForecasterID != "UNKNOWN" &
                                 !is.na(Issued) &
                                 !is.na(validTime) &
                                 !is.na(aerodrome)]

analysis_data[, `:=`(
  lead_hours = as.numeric(difftime(validTime, Issued, units = "hours")),
  
  taf_id = interaction(aerodrome, format(Issued, "%Y%m%d%H%M%S"), drop = TRUE),
  
  verify_case_id = interaction(aerodrome, format(validTime, "%Y%m%d%H"), drop = TRUE),
  
  issue_cycle = factor(floor(hour(Issued) / 3) * 3, levels = seq(0, 21, 3)),
  
  valid_hour = factor(hour(validTime)),
  month = factor(month(validTime))
)]

analysis_data <- analysis_data[lead_hours >= 0 &
                                 lead_hours <= 30]

setorder(analysis_data, aerodrome, validTime, Issued)

analysis_data[, `:=`(
  forecast_sequence = seq_len(.N),
  n_forecasts_for_case = .N,
  is_latest_forecast = seq_len(.N) == .N,
  previous_forecaster = shift(ForecasterID),
  previous_issue = shift(Issued)
), by = verify_case_id]

analysis_data[, hours_since_previous :=
                as.numeric(difftime(Issued, previous_issue, units = "hours"))]

# Check that each forecast vintage is unique:
duplicate_check <- analysis_data[, .N, by = .(aerodrome, Issued, validTime)][N > 1]
duplicate_check

# Check TAF Amendments
table(noise_data$is_amd, useNA = "ifany")
# By aerodrome
noise_data[, .(
  n_rows = .N,
  amd_rows = sum(is_amd, na.rm = TRUE),
  amd_rate = mean(is_amd, na.rm = TRUE)
), by = aerodrome][order(-amd_rate)]

taf_amendment_summary <- unique(noise_data[, .(aerodrome, Issued, ForecasterID, is_amd)])[, .(
  n_tafs = .N,
  n_amd = sum(is_amd, na.rm = TRUE),
  amd_rate = mean(is_amd, na.rm = TRUE)
), by = aerodrome][order(-amd_rate)]

taf_amendment_summary


# 3. Build a long-format performance dataset
performance_data <- rbindlist(
  list(
    analysis_data[!is.na(wind_spd_error), .(
      ForecasterID,
      aerodrome,
      Issued,
      validTime,
      taf_id,
      verify_case_id,
      lead_hours,
      issue_cycle,
      valid_hour,
      month,
      forecast_sequence,
      n_forecasts_for_case,
      is_latest_forecast,
      is_amd,
      element = "wind_speed",
      signed_error = wind_spd_error,
      absolute_error = abs(wind_spd_error),
      accurate = as.integer(abs(wind_spd_error) <= 5)
    )],
    
    analysis_data[!is.na(wind_dir_error), .(
      ForecasterID,
      aerodrome,
      Issued,
      validTime,
      taf_id,
      verify_case_id,
      lead_hours,
      issue_cycle,
      valid_hour,
      month,
      forecast_sequence,
      n_forecasts_for_case,
      is_latest_forecast,
      is_amd,
      element = "wind_direction",
      signed_error = wind_dir_error,
      absolute_error = abs(wind_dir_error),
      accurate = as.integer(wind_dir_good)
    )],
    
    analysis_data[!is.na(cloud_base_error), .(
      ForecasterID,
      aerodrome,
      Issued,
      validTime,
      taf_id,
      verify_case_id,
      lead_hours,
      issue_cycle,
      valid_hour,
      month,
      forecast_sequence,
      n_forecasts_for_case,
      is_latest_forecast,
      is_amd,
      element = "cloud_base",
      signed_error = cloud_base_error,
      absolute_error = abs(cloud_base_error),
      accurate = as.integer(abs(cloud_base_error) <= 500)
    )],
    
    analysis_data[!is.na(cloud_amount_error), .(
      ForecasterID,
      aerodrome,
      Issued,
      validTime,
      taf_id,
      verify_case_id,
      lead_hours,
      issue_cycle,
      valid_hour,
      month,
      forecast_sequence,
      n_forecasts_for_case,
      is_latest_forecast,
      is_amd,
      element = "cloud_amount",
      signed_error = cloud_amount_error,
      absolute_error = abs(cloud_amount_error),
      accurate = as.integer(abs(cloud_amount_error) <= 1)
    )]
  ),
  use.names = TRUE,
  fill = TRUE
)

performance_data[, element := factor(element,
                                     levels = c("wind_speed", "wind_direction", "cloud_base", "cloud_amount"))]

performance_data[, lead_z :=
                   as.numeric(scale(lead_hours))]

performance_data[, sequence_z :=
                   as.numeric(scale(forecast_sequence))]

# Temperature is added separately
temp_performance <- copy(temp_data)

temp_performance <- temp_performance[!is.na(temp_error) &
                                       !is.na(ForecasterID) &
                                       ForecasterID != "UNKNOWN"]

temp_performance[, `:=`(
  lead_hours = as.numeric(difftime(validTime, Issued, units = "hours")),
  
  taf_id = interaction(aerodrome, format(Issued, "%Y%m%d%H%M%S"), drop = TRUE),
  
  verify_case_id = interaction(aerodrome, format(validTime, "%Y%m%d%H"), drop = TRUE),
  
  issue_cycle = factor(floor(hour(Issued) / 3) * 3, levels = seq(0, 21, 3)),
  
  valid_hour = factor(hour(validTime)),
  month = factor(month(validTime)),
  
  element = factor(
    "temperature",
    levels = c(
      "wind_speed",
      "wind_direction",
      "cloud_base",
      "cloud_amount",
      "temperature"
    )
  ),
  
  signed_error = temp_error,
  absolute_error = abs(temp_error),
  accurate = as.integer(abs(temp_error) <= 1)
)]

setorder(temp_performance, aerodrome, validTime, Issued)

temp_performance[, `:=`(
  forecast_sequence = seq_len(.N),
  n_forecasts_for_case = .N,
  is_latest_forecast = seq_len(.N) == .N
), by = verify_case_id]

temp_performance[, `:=`(lead_z = as.numeric(scale(lead_hours)), sequence_z = as.numeric(scale(forecast_sequence)))]

# Combine:
common_columns <- intersect(names(performance_data), names(temp_performance))

performance_data <- rbindlist(list(performance_data[, ..common_columns], temp_performance[, ..common_columns]),
                              use.names = TRUE,
                              fill = TRUE)

performance_data[, element := factor(element)]

# 4. Descriptive check:
# does accuracy improve as lead time contracts?
performance_data[, lead_bin := cut(
  lead_hours,
  breaks = c(0, 3, 6, 9, 12, 18, 24, 30),
  include.lowest = TRUE,
  right = FALSE
)]

lead_summary <- performance_data[, .(
  n = .N,
  accuracy = mean(accurate, na.rm = TRUE),
  mae = mean(absolute_error, na.rm = TRUE),
  bias = mean(signed_error, na.rm = TRUE)
), by = .(element, lead_bin)]

lead_summary

##################################################
##### HIERARCHICAL PERFORMANCE MODEL
##################################################
hier_data <- copy(performance_data)

# Ensure identifiers are factors
hier_data[, `:=`(
  ForecasterID = factor(ForecasterID),
  aerodrome = factor(aerodrome),
  taf_id = factor(taf_id),
  verify_case_id = factor(verify_case_id),
  element = factor(element),
  issue_cycle = factor(issue_cycle),
  valid_hour = factor(valid_hour),
  month = factor(month)
)]

# Forecaster-by-element identifier:
# this represents Stable Pattern Noise across forecast elements
hier_data[, forecaster_element := interaction(ForecasterID, element, drop = TRUE)]

# Binary response for MCMCglmm
hier_data[, accuracy_factor := factor(accurate,
                                      levels = c(0, 1),
                                      labels = c("not_accurate", "accurate"))]

# Keep complete model records
hier_data <- hier_data[!is.na(accuracy_factor) &
                         !is.na(lead_z) &
                         !is.na(sequence_z) &
                         !is.na(ForecasterID) &
                         !is.na(aerodrome) &
                         !is.na(taf_id) &
                         !is.na(verify_case_id)]

droplevels(hier_data)

hier_data[, .(
  n = .N,
  accuracy = mean(accurate),
  n_forecasters = uniqueN(ForecasterID),
  n_airports = uniqueN(aerodrome),
  n_tafs = uniqueN(taf_id),
  n_verification_cases = uniqueN(verify_case_id)
), by = element]

# 3. Check whether AMD is available across all elements
"is_amd" %in% names(hier_data)

if ("is_amd" %in% names(hier_data)) {
  hier_data[, is_amd := factor(is_amd,
                               levels = c(FALSE, TRUE),
                               labels = c("original", "amended"))]
  
  table(hier_data$element, hier_data$is_amd, useNA = "ifany")
}

install.packages("blme")

library(blme)
library(lme4)

# 3. Create a manageable difficulty block
fast_data <- copy(hier_data)

fast_data[, valid_date := as.Date(validTime)]

fast_data[, weather_block := interaction(aerodrome, valid_date, drop = TRUE)]

fast_data[, `:=`(
  ForecasterID = factor(ForecasterID),
  forecaster_element = factor(forecaster_element),
  aerodrome = factor(aerodrome),
  weather_block = factor(weather_block),
  element = factor(element),
  issue_cycle = factor(issue_cycle),
  valid_hour = factor(valid_hour),
  month = factor(month)
)]

fast_data <- fast_data[complete.cases(
  accurate,
  element,
  lead_z,
  sequence_z,
  issue_cycle,
  valid_hour,
  month,
  ForecasterID,
  forecaster_element,
  aerodrome,
  weather_block
)]

fast_data[, .(
  n = .N,
  n_forecasters = uniqueN(ForecasterID),
  n_weather_blocks = uniqueN(weather_block),
  accuracy = mean(accurate)
), by = element]

# 4. Fit the fast hierarchical model
# The priors are explicit:
# fixed effects: normal prior centred on zero;
# random-effect standard deviations: Wishart-type regularising prior;
# fixed-effect prior scale: variance =4=4=4, equivalent to SD =2=2=2 on the log-odds scale.
set.seed(20261015)

system.time({
  fit_accuracy_fast <- bglmer(
    accurate ~
      element * lead_z +
      element * sequence_z +
      issue_cycle +
      valid_hour +
      month +
      (1 | ForecasterID) +
      (1 | forecaster_element) +
      (1 | aerodrome) +
      (1 | weather_block),
    
    data = fast_data,
    
    family = binomial(link = "logit"),
    
    # bglmer defaults:
    # cov.prior   = wishart
    # fixef.prior = NULL
    
    nAGQ = 0,
    
    control = glmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 100000),
      calc.derivs = FALSE
    )
  )
})

summary(fit_accuracy_fast)

fit_accuracy_fast@optinfo$conv$lme4$messages

as.data.frame(VarCorr(fit_accuracy_fast))

continuous_model <- fit_accuracy_fast

# Confirm factor order and reference category
levels(model.frame(continuous_model)$element)
contrasts(model.frame(continuous_model)$element)

extract_element_effect <- function(model,
                                   element_name,
                                   effect_name,
                                   main_term,
                                   interaction_term = NULL) {
  beta <- fixef(model)
  V <- as.matrix(vcov(model))
  
  L <- setNames(rep(0, length(beta)), names(beta))
  
  L[main_term] <- 1
  
  if (!is.null(interaction_term)) {
    if (!interaction_term %in% names(beta)) {
      stop("Coefficient not found: ", interaction_term)
    }
    
    L[interaction_term] <- 1
  }
  
  estimate <- sum(L * beta)
  
  std_error <- sqrt(as.numeric(t(L) %*% V %*% L))
  
  tibble(
    element = element_name,
    effect = effect_name,
    estimate = estimate,
    std_error = std_error,
    odds_ratio = exp(estimate),
    lower_95 = exp(estimate - 1.96 * std_error),
    upper_95 = exp(estimate + 1.96 * std_error)
  )
}

continuous_element_effects <- bind_rows(
  # Wind speed is the reference element
  extract_element_effect(
    continuous_model,
    element_name = "Wind speed",
    effect_name = "Lead time",
    main_term = "lead_z"
  ),
  
  extract_element_effect(
    continuous_model,
    element_name = "Wind speed",
    effect_name = "Forecast vintage",
    main_term = "sequence_z"
  ),
  
  # Wind direction
  extract_element_effect(
    continuous_model,
    element_name = "Wind direction",
    effect_name = "Lead time",
    main_term = "lead_z",
    interaction_term =
      "elementwind_direction:lead_z"
  ),
  
  extract_element_effect(
    continuous_model,
    element_name = "Wind direction",
    effect_name = "Forecast vintage",
    main_term = "sequence_z",
    interaction_term =
      "elementwind_direction:sequence_z"
  ),
  
  # Cloud base
  extract_element_effect(
    continuous_model,
    element_name = "Cloud base",
    effect_name = "Lead time",
    main_term = "lead_z",
    interaction_term =
      "elementcloud_base:lead_z"
  ),
  
  extract_element_effect(
    continuous_model,
    element_name = "Cloud base",
    effect_name = "Forecast vintage",
    main_term = "sequence_z",
    interaction_term =
      "elementcloud_base:sequence_z"
  ),
  
  # Cloud amount
  extract_element_effect(
    continuous_model,
    element_name = "Cloud amount",
    effect_name = "Lead time",
    main_term = "lead_z",
    interaction_term =
      "elementcloud_amount:lead_z"
  ),
  
  extract_element_effect(
    continuous_model,
    element_name = "Cloud amount",
    effect_name = "Forecast vintage",
    main_term = "sequence_z",
    interaction_term =
      "elementcloud_amount:sequence_z"
  ),
  
  # Temperature
  extract_element_effect(
    continuous_model,
    element_name = "Temperature",
    effect_name = "Lead time",
    main_term = "lead_z",
    interaction_term =
      "elementtemperature:lead_z"
  ),
  
  extract_element_effect(
    continuous_model,
    element_name = "Temperature",
    effect_name = "Forecast vintage",
    main_term = "sequence_z",
    interaction_term =
      "elementtemperature:sequence_z"
  )
)

continuous_element_effects |>
  print(n = Inf, width = Inf)

continuous_effects_table <- continuous_element_effects |>
  mutate(`OR (95% interval)` = sprintf("%.3f (%.3f–%.3f)", odds_ratio, lower_95, upper_95)) |>
  select(
    Element = element,
    Effect = effect,
    Estimate = estimate,
    `Standard error` = std_error,
    `OR (95% interval)`
  )

continuous_effects_table |>
  print(n = Inf, width = Inf)

continuous_effects_wide <- continuous_element_effects |>
  mutate(result = sprintf("%.3f (%.3f–%.3f)", odds_ratio, lower_95, upper_95)) |>
  select(element, effect, result) |>
  tidyr::pivot_wider(names_from = effect, values_from = result) |>
  rename(
    Element = element,
    `Lead-time OR (95% interval)` =
      `Lead time`,
    `Vintage OR (95% interval)` =
      `Forecast vintage`
  )

continuous_effects_wide

##################################################
# HIERARCHICAL VARIANCE DECOMPOSITION
##################################################
# Extract random-effect variance components
vc <- as.data.table(as.data.frame(VarCorr(fit_accuracy_fast)))

# Retain the grouping factor, variance and standard deviation
variance_results <- vc[, .(component = grp,
                           variance = vcov,
                           sd = sdcor)]

# Add the conventional latent residual variance
# for a logistic mixed model
variance_results <- rbind(
  variance_results,
  data.table(
    component = "latent_residual",
    variance = pi^2 / 3,
    sd = sqrt(pi^2 / 3)
  ),
  use.names = TRUE
)

# Calculate total latent-scale variance
total_variance <- sum(variance_results$variance, na.rm = TRUE)

# Calculate each component's proportion
variance_results[, proportion := variance / total_variance]

# Add percentages
variance_results[, percentage := 100 * proportion]

# Order from largest to smallest
setorder(variance_results, -proportion)

# Display results
variance_results

variance_results[, label := fifelse(
  component == "ForecasterID",
  "General forecaster effect (Level Noise)",
  fifelse(
    component == "forecaster_element",
    "Forecaster-by-element effect (Stable Pattern Noise)",
    fifelse(
      component == "aerodrome",
      "Aerodrome effect",
      fifelse(
        component == "weather_block",
        "Aerodrome-day weather difficulty",
        "Latent residual variation"
      )
    )
  )
)]

variance_table <- variance_results[, .(
  Component = label,
  Variance = round(variance, 4),
  `Standard deviation` = round(sd, 4),
  `Variance share (%)` = round(percentage, 2)
)]

variance_table

### MEDIAN ODDS RATIO
mor_results <- variance_results[mor_results <- variance != "latent_residual"]

mor_results[, median_odds_ratio := exp(0.6745 * sqrt(2 * variance))]

mor_table <- mor_results[, .(
  Component = label,
  Variance = round(variance, 4),
  `Median odds ratio` = round(median_odds_ratio, 3)
)]

mor_table

##################################################
##### OPERATIONALLY SIGNIFICANT EVENT ANALYSIS
##################################################

event_base <- copy(analysis_data)

# Confirm the visibility threshold encoded by vis_fcst.
# Your extract_block_hazards() currently uses < 5000 m,
# so observed low visibility must use the same threshold.

event_base[, `:=`(
  obs_fg = as.integer(event_fg),
  
  obs_ts = as.integer(event_ts),
  
  obs_low_vis = as.integer(!is.na(vis) &
                             vis < 5000),
  
  # Match the forecast definition:
  # BKN/OVC below 1000 ft
  obs_low_cloud = as.integer(
    !is.na(cloud_base) &
      !is.na(cloud_amount) &
      cloud_base < 1000 &
      cloud_amount >= 3
  ),
  
  fcst_fg = as.integer(fg_fcst),
  
  fcst_ts = as.integer(ts_fcst),
  
  fcst_low_vis = as.integer(vis_fcst),
  
  fcst_low_cloud = as.integer(cld1000_fcst),
  
  valid_date = as.Date(validTime),
  
  weather_block = interaction(aerodrome, as.Date(validTime), drop = TRUE)
)]

event_data <- rbindlist(list(event_base[, .(
  ForecasterID,
  aerodrome,
  Issued,
  validTime,
  taf_id,
  verify_case_id,
  weather_block,
  lead_hours,
  forecast_sequence,
  issue_cycle,
  valid_hour,
  month,
  is_amd,
  forecast_type,
  event = "fog",
  forecast = fcst_fg,
  observed = obs_fg
)], event_base[, .(
  ForecasterID,
  aerodrome,
  Issued,
  validTime,
  taf_id,
  verify_case_id,
  weather_block,
  lead_hours,
  forecast_sequence,
  issue_cycle,
  valid_hour,
  month,
  is_amd,
  forecast_type,
  event = "thunderstorm",
  forecast = fcst_ts,
  observed = obs_ts
)], event_base[, .(
  ForecasterID,
  aerodrome,
  Issued,
  validTime,
  taf_id,
  verify_case_id,
  weather_block,
  lead_hours,
  forecast_sequence,
  issue_cycle,
  valid_hour,
  month,
  is_amd,
  forecast_type,
  event = "low_visibility",
  forecast = fcst_low_vis,
  observed = obs_low_vis
)], event_base[, .(
  ForecasterID,
  aerodrome,
  Issued,
  validTime,
  taf_id,
  verify_case_id,
  weather_block,
  lead_hours,
  forecast_sequence,
  issue_cycle,
  valid_hour,
  month,
  is_amd,
  forecast_type,
  event = "low_cloud",
  forecast = fcst_low_cloud,
  observed = obs_low_cloud
)]),
use.names = TRUE)

# CONVERT TO LONG FORMAT
event_data[, `:=`(
  event = factor(
    event,
    levels = c("fog", "thunderstorm", "low_visibility", "low_cloud")
  ),
  
  ForecasterID = factor(ForecasterID),
  aerodrome = factor(aerodrome),
  weather_block = factor(weather_block),
  issue_cycle = factor(issue_cycle),
  valid_hour = factor(valid_hour),
  month = factor(month),
  is_amd = factor(
    is_amd,
    levels = c(FALSE, TRUE),
    labels = c("original", "amended")
  ),
  
  lead_z = as.numeric(scale(lead_hours)),
  sequence_z = as.numeric(scale(forecast_sequence))
)]

# Create verification outcomes
event_data[, `:=`(
  hit =
    as.integer(forecast == 1 & observed == 1),
  
  miss =
    as.integer(forecast == 0 & observed == 1),
  
  false_alarm =
    as.integer(forecast == 1 & observed == 0),
  
  correct_negative =
    as.integer(forecast == 0 & observed == 0),
  
  correct =
    as.integer(forecast == observed),
  
  signed_error =
    forecast - observed
)]

# Here:
#   signed_error = 1 is a false alarm;
#   signed_error = -1 is a miss;
#   signed_error = 0 is correct.

# Produce the event contingency table
event_contingency <- event_data[, .(
  n = .N,
  
  observed_events =
    sum(observed == 1, na.rm = TRUE),
  
  forecast_events =
    sum(forecast == 1, na.rm = TRUE),
  
  hits =
    sum(hit, na.rm = TRUE),
  
  misses =
    sum(miss, na.rm = TRUE),
  
  false_alarms =
    sum(false_alarm, na.rm = TRUE),
  
  correct_negatives =
    sum(correct_negative, na.rm = TRUE)
  
), by = event]

event_contingency

# Calculate conventional verification statistics
event_scores <- event_contingency[, `:=`(
  observed_rate =
    observed_events / n,
  
  forecast_rate =
    forecast_events / n,
  
  accuracy =
    (hits + correct_negatives) / n,
  
  probability_of_detection =
    fifelse(hits + misses > 0, hits / (hits + misses), NA_real_),
  
  false_alarm_ratio =
    fifelse(
      hits + false_alarms > 0,
      false_alarms / (hits + false_alarms),
      NA_real_
    ),
  
  probability_of_false_detection =
    fifelse(
      false_alarms + correct_negatives > 0,
      false_alarms /
        (false_alarms + correct_negatives),
      NA_real_
    ),
  
  critical_success_index =
    fifelse(
      hits + misses + false_alarms > 0,
      hits /
        (hits + misses + false_alarms),
      NA_real_
    ),
  
  frequency_bias =
    fifelse(hits + misses > 0, (hits + false_alarms) /
              (hits + misses), NA_real_),
  
  equitable_threat_score = {
    random_hits <-
      (hits + misses) *
      (hits + false_alarms) / n
    
    fifelse(
      hits + misses + false_alarms -
        random_hits > 0,
      
      (hits - random_hits) /
        (hits +
           misses +
           false_alarms -
           random_hits),
      
      NA_real_
    )
  },
  
  signed_bias =
    forecast_rate - observed_rate
)]

event_scores

# Calculate scores by forecaster
forecaster_event_scores <- event_data[, {
  h <- sum(hit, na.rm = TRUE)
  m <- sum(miss, na.rm = TRUE)
  f <- sum(false_alarm, na.rm = TRUE)
  c <- sum(correct_negative, na.rm = TRUE)
  n_total <- .N
  
  random_hits <-
    (h + m) * (h + f) / n_total
  
  .(
    n = n_total,
    observed_events = h + m,
    forecast_events = h + f,
    
    probability_of_detection =
      fifelse(h + m > 0, h / (h + m), NA_real_),
    
    false_alarm_ratio =
      fifelse(h + f > 0, f / (h + f), NA_real_),
    
    critical_success_index =
      fifelse(h + m + f > 0, h / (h + m + f), NA_real_),
    
    equitable_threat_score =
      fifelse(
        h + m + f - random_hits > 0,
        (h - random_hits) /
          (h + m + f - random_hits),
        NA_real_
      ),
    
    frequency_bias =
      fifelse(h + m > 0, (h + f) / (h + m), NA_real_)
  )
  
}, by = .(event, ForecasterID)]

forecaster_event_support <- event_data[, .(
  n = .N,
  observed_events = sum(observed),
  forecast_events = sum(forecast),
  n_weather_blocks = uniqueN(weather_block)
), by = .(event, ForecasterID)]

forecaster_event_support[, summary(observed_events), by = event]

# Hierarchical accuracy model for each event
fit_event_accuracy <- function(event_name) {
  d <- droplevels(event_data[event == event_name &
                               complete.cases(
                                 correct,
                                 lead_z,
                                 sequence_z,
                                 issue_cycle,
                                 valid_hour,
                                 month,
                                 is_amd,
                                 ForecasterID,
                                 aerodrome,
                                 weather_block
                               )])
  
  message(
    "Fitting ",
    event_name,
    ": n = ",
    nrow(d),
    ", observed events = ",
    sum(d$observed),
    ", forecast events = ",
    sum(d$forecast)
  )
  
  blme::bglmer(
    correct ~
      lead_z +
      sequence_z +
      is_amd +
      issue_cycle +
      valid_hour +
      month +
      forecast_type +
      
      (1 | ForecasterID) +
      (1 | aerodrome) +
      (1 | weather_block),
    
    data = d,
    
    family = binomial(link = "logit"),
    
    nAGQ = 0,
    
    control = lme4::glmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 100000),
      calc.derivs = FALSE
    )
  )
}

event_names <- levels(event_data$event)

event_accuracy_models <- lapply(event_names, fit_event_accuracy)

names(event_accuracy_models) <- event_names

saveRDS(event_accuracy_models,
        "event_accuracy_hierarchical_models.rds")

##################################################
# EXTRACT FIXED EFFECTS
##################################################
extract_fixed_effects <- function(model, event_name) {
  if (is.null(model)) {
    return(NULL)
  }
  
  x <- as.data.table(coef(summary(model)), keep.rownames = "parameter")
  
  # Rename according to the actual bglmer summary columns
  setnames(
    x,
    old = names(x),
    new = c(
      "parameter",
      "estimate",
      "standard_error",
      "z_value",
      "p_value"
    )
  )
  
  # First create confidence limits
  x[, `:=`(
    event = event_name,
    lower_95 =
      estimate - 1.96 * standard_error,
    upper_95 =
      estimate + 1.96 * standard_error
  )]
  
  # Then transform estimates and limits to odds ratios
  x[, `:=`(
    odds_ratio =
      exp(estimate),
    
    odds_ratio_lower_95 =
      exp(lower_95),
    
    odds_ratio_upper_95 =
      exp(upper_95)
  )]
  
  # Put event first
  setcolorder(
    x,
    c(
      "event",
      "parameter",
      "estimate",
      "standard_error",
      "lower_95",
      "upper_95",
      "odds_ratio",
      "odds_ratio_lower_95",
      "odds_ratio_upper_95",
      "z_value",
      "p_value"
    )
  )
  
  x[]
}

event_accuracy_fixed <- rbindlist(
  Map(
    extract_fixed_effects,
    event_accuracy_models,
    names(event_accuracy_models)
  ),
  use.names = TRUE,
  fill = TRUE
)

head(event_accuracy_fixed)

unique(event_accuracy_fixed$event)

event_accuracy_fixed[, .N, by = event]

event_accuracy_main_effects <- event_accuracy_fixed[parameter %in% c("lead_z", "sequence_z", "is_amdamended"), .(
  event,
  parameter,
  estimate,
  standard_error,
  lower_95,
  upper_95,
  odds_ratio,
  odds_ratio_lower_95,
  odds_ratio_upper_95,
  p_value
)]

event_accuracy_main_effects

event_accuracy_table <- copy(event_accuracy_main_effects)

event_accuracy_table[, effect := fcase(
  parameter == "lead_z",
  "Lead time",
  
  parameter == "sequence_z",
  "Forecast vintage",
  
  parameter == "is_amdamended",
  "Amended TAF",
  
  default = parameter
)]

event_accuracy_table[, `:=`(
  odds_ratio = round(odds_ratio, 3),
  
  odds_ratio_lower_95 = round(odds_ratio_lower_95, 3),
  
  odds_ratio_upper_95 = round(odds_ratio_upper_95, 3),
  
  p_value = signif(p_value, 3)
)]

event_accuracy_table[, .(
  Event = event,
  Effect = effect,
  `Odds ratio` = odds_ratio,
  `Lower 95% limit` =
    odds_ratio_lower_95,
  `Upper 95% limit` =
    odds_ratio_upper_95,
  `p-value` = p_value
)]


event_accuracy_fixed[parameter %in% c("lead_z", "sequence_z", "is_amdamended")]

###
# Extract event-specific forecaster noise
extract_event_variance <- function(model, event_name) {
  vc <- as.data.table(as.data.frame(VarCorr(model)))
  
  out <- vc[, .(
    event = event_name,
    component = grp,
    variance = vcov,
    sd = sdcor
  )]
  
  out <- rbind(
    out,
    data.table(
      event = event_name,
      component = "latent_residual",
      variance = pi^2 / 3,
      sd = sqrt(pi^2 / 3)
    )
  )
  
  out[, `:=`(
    proportion =
      variance / sum(variance),
    
    percentage =
      100 * variance / sum(variance)
  )]
  
  out[component != "latent_residual", median_odds_ratio :=
        exp(0.6745 *
              sqrt(2 * variance))]
  
  out[]
}

event_variance <- rbindlist(Map(
  extract_event_variance,
  event_accuracy_models,
  names(event_accuracy_models)
),
fill = TRUE)

event_variance[component %in% c("ForecasterID", "aerodrome", "weather_block")][order(event, -percentage)]

# Key tables
event_forecaster_noise <- event_variance[component == "ForecasterID", .(
  event,
  forecaster_variance = variance,
  forecaster_sd = sd,
  variance_share_percent = percentage,
  median_odds_ratio
)]

event_forecaster_noise

# Model misses conditional on the event occurring
# more operationally meaningful than overall accuracy

fit_event_miss <- function(event_name) {
  d <- droplevels(event_data[event == event_name &
                               observed == 1 &
                               complete.cases(
                                 miss,
                                 lead_z,
                                 sequence_z,
                                 issue_cycle,
                                 valid_hour,
                                 month,
                                 is_amd,
                                 ForecasterID,
                                 aerodrome,
                                 weather_block
                               )])
  
  message("Miss model: ",
          event_name,
          "; n observed events = ",
          nrow(d),
          "; misses = ",
          sum(d$miss))
  
  if (nrow(d) < 100 || length(unique(d$miss)) < 2)
    return(NULL)
  
  blme::bglmer(
    miss ~
      lead_z +
      sequence_z +
      is_amd +
      issue_cycle +
      valid_hour +
      month +
      forecast_type +
      
      (1 | ForecasterID) +
      (1 | aerodrome) +
      (1 | weather_block),
    
    data = d,
    
    family = binomial(link = "logit"),
    
    nAGQ = 0,
    
    control = glmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 100000),
      calc.derivs = FALSE
    )
  )
}

# Fit
event_miss_models <- lapply(event_names, fit_event_miss)

names(event_miss_models) <- event_names
event_miss_models

miss_model_table <- purrr::imap_dfr(event_miss_models, function(model, event_name) {
  beta <- lme4::fixef(model)
  se   <- sqrt(diag(vcov(model)))
  
  terms <- c(lead_z = "Lead time",
             sequence_z = "Forecast vintage",
             is_amdamended = "Amended TAF")
  
  tibble::tibble(
    event = event_name,
    term = names(terms),
    effect = unname(terms),
    estimate = beta[names(terms)],
    std_error = se[names(terms)],
    odds_ratio = exp(estimate),
    lower_95 = exp(estimate - 1.96 * std_error),
    upper_95 = exp(estimate + 1.96 * std_error)
  )
})

miss_model_table

# Interpretation:
# positive lead_z: higher miss risk at longer lead time;
# negative sequence_z: later forecasts reduce miss risk;
# negative is_amdamended: amendments reduce miss risk.

# Model false alarms conditional on no event occurring
fit_event_false_alarm <- function(event_name) {
  d <- droplevels(event_data[event == event_name &
                               observed == 0 &
                               complete.cases(
                                 false_alarm,
                                 lead_z,
                                 sequence_z,
                                 issue_cycle,
                                 valid_hour,
                                 month,
                                 is_amd,
                                 ForecasterID,
                                 aerodrome,
                                 weather_block
                               )])
  
  message(
    "False-alarm model: ",
    event_name,
    "; n non-events = ",
    nrow(d),
    "; false alarms = ",
    sum(d$false_alarm)
  )
  
  if (nrow(d) < 100 || length(unique(d$false_alarm)) < 2)
    return(NULL)
  
  blme::bglmer(
    false_alarm ~
      lead_z +
      sequence_z +
      is_amd +
      issue_cycle +
      valid_hour +
      month +
      forecast_type +
      
      (1 | ForecasterID) +
      (1 | aerodrome) +
      (1 | weather_block),
    
    data = d,
    
    family = binomial(link = "logit"),
    
    nAGQ = 0,
    
    control = glmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 100000),
      calc.derivs = FALSE
    )
  )
}

# Fit
event_false_alarm_models <- lapply(event_names, fit_event_false_alarm)

names(event_false_alarm_models) <- event_names


library(dplyr)
library(purrr)
library(tibble)
library(lme4)

# Terms to retain in the publication table
focal_terms <- c(lead_z        = "Lead time",
                 sequence_z    = "Forecast vintage",
                 is_amdamended = "Amended TAF")

false_alarm_model_table <- imap_dfr(event_false_alarm_models, function(model, event_name) {
  # Fixed-effect estimates and model-based standard errors
  beta <- fixef(model)
  se   <- sqrt(diag(vcov(model)))
  
  # Confirm that all required coefficients exist
  missing_terms <- setdiff(names(focal_terms), names(beta))
  
  if (length(missing_terms) > 0) {
    stop(
      "Missing coefficient(s) in ",
      event_name,
      " model: ",
      paste(missing_terms, collapse = ", ")
    )
  }
  
  # Model frame used in the fitted model
  model_data <- model.frame(model)
  
  # False-alarm response must be coded 0/1
  if (!"false_alarm" %in% names(model_data)) {
    stop(
      "The response variable 'false_alarm' was not found in the ",
      event_name,
      " model frame."
    )
  }
  
  # Forecaster random-effect standard deviation
  forecaster_sd <- as.data.frame(VarCorr(model)) |>
    filter(grp == "ForecasterID") |>
    pull(sdcor)
  
  if (length(forecaster_sd) != 1) {
    stop(
      "Could not uniquely identify the ForecasterID random-effect SD ",
      "for the ",
      event_name,
      " model."
    )
  }
  
  # Median odds ratio
  forecaster_mor <- exp(qnorm(0.75) * sqrt(2) * forecaster_sd)
  
  tibble(
    event = event_name,
    
    # Counts
    non_event_hours = nobs(model),
    false_alarms = sum(model_data$false_alarm == 1, na.rm = TRUE),
    
    # Lead-time effect
    lead_estimate = unname(beta["lead_z"]),
    lead_std_error = unname(se["lead_z"]),
    lead_odds_ratio = exp(beta["lead_z"]),
    lead_lower_95 = exp(beta["lead_z"] - 1.96 * se["lead_z"]),
    lead_upper_95 = exp(beta["lead_z"] + 1.96 * se["lead_z"]),
    
    # Forecast-vintage effect
    vintage_estimate = unname(beta["sequence_z"]),
    vintage_std_error = unname(se["sequence_z"]),
    vintage_odds_ratio = exp(beta["sequence_z"]),
    vintage_lower_95 = exp(beta["sequence_z"] - 1.96 * se["sequence_z"]),
    vintage_upper_95 = exp(beta["sequence_z"] + 1.96 * se["sequence_z"]),
    
    # Amendment effect
    amendment_estimate = unname(beta["is_amdamended"]),
    amendment_std_error = unname(se["is_amdamended"]),
    amendment_odds_ratio = exp(beta["is_amdamended"]),
    amendment_lower_95 = exp(beta["is_amdamended"] -
                               1.96 * se["is_amdamended"]),
    amendment_upper_95 = exp(beta["is_amdamended"] +
                               1.96 * se["is_amdamended"]),
    
    # Forecaster dispersion
    forecaster_sd = forecaster_sd,
    forecaster_mor = forecaster_mor
  )
})

false_alarm_model_table

table5_output <- false_alarm_model_table |>
  mutate(
    event = recode(
      event,
      fog = "Fog",
      thunderstorm = "Thunderstorm",
      low_visibility = "Low visibility",
      low_cloud = "Low cloud"
    ),
    
    `Lead-time OR (95% interval)` = sprintf(
      "%.3f (%.3f–%.3f)",
      lead_odds_ratio,
      lead_lower_95,
      lead_upper_95
    ),
    
    `Vintage OR (95% interval)` = sprintf(
      "%.3f (%.3f–%.3f)",
      vintage_odds_ratio,
      vintage_lower_95,
      vintage_upper_95
    ),
    
    `Amendment OR (95% interval)` = sprintf(
      "%.3f (%.3f–%.3f)",
      amendment_odds_ratio,
      amendment_lower_95,
      amendment_upper_95
    ),
    
    `Non-event hours` = scales::comma(non_event_hours),
    `False alarms` = scales::comma(false_alarms),
    `Forecaster SD` = sprintf("%.3f", forecaster_sd),
    `Forecaster MOR` = sprintf("%.2f", forecaster_mor)
  ) |>
  select(
    Event = event,
    `Non-event hours`,
    `False alarms`,
    `Lead-time OR (95% interval)`,
    `Vintage OR (95% interval)`,
    `Amendment OR (95% interval)`,
    `Forecaster SD`,
    `Forecaster MOR`
  )

table5_output

false_alarm_model_table |>
  select(
    event,
    amendment_odds_ratio,
    amendment_lower_95,
    amendment_upper_95,
    forecaster_sd,
    forecaster_mor
  ) |>
  print(n = Inf, width = Inf)

# Interpretation:
# positive forecaster random-effect variance means forecasters differ in precautionary tendency after case adjustment;
# negative sequence_z means later vintages reduce false alarms;
# a negative amendment effect means amendments remove unnecessary event signals.

###
# Forecast-propensity model: the cleanest bias model
# models whether the event was forecast, conditional on whether it occurred.
library(parallel)

setDT(event_data)

set.seed(2026)

# Development sample:
# retain up to 10,000 records per event × observed × forecast stratum
event_data_sample <- event_data[complete.cases(
  forecast,
  observed,
  lead_z,
  sequence_z,
  issue_cycle,
  valid_hour,
  month,
  is_amd,
  forecast_type,
  ForecasterID,
  aerodrome,
  weather_block
), .SD[sample.int(.N, min(.N, 10000L))], by = .(event, observed, forecast, aerodrome, is_amd)]

fit_event_bias <- function(event_name, data = event_data_sample) {
  d <- droplevels(data[event == event_name])
  
  message("Fitting ",
          event_name,
          ": ",
          format(nrow(d), big.mark = ","),
          " records")
  
  blme::bglmer(
    forecast ~
      observed * lead_z +
      observed * sequence_z +
      observed * is_amd +
      issue_cycle +
      valid_hour +
      month +
      forecast_type +
      (1 | ForecasterID) +
      (1 | aerodrome) +
      (1 | weather_block),
    
    data = d,
    family = binomial(link = "logit"),
    nAGQ = 0,
    
    control = glmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 30000),
      calc.derivs = FALSE
    )
  )
}

# fit the event models in parallel:
n_workers <- max(1L, min(length(event_names), detectCores() - 1L))

if (.Platform$OS.type == "windows") {
  cl <- makeCluster(n_workers)
  
  clusterExport(cl,
                c("event_names", "event_data_sample", "fit_event_bias"),
                envir = environment())
  
  clusterEvalQ(cl, {
    library(data.table)
    library(blme)
    library(lme4)
  })
  
  event_bias_models <- parLapply(cl, event_names, fit_event_bias, data = event_data_sample)
  
  stopCluster(cl)
  
} else {
  event_bias_models <- mclapply(event_names,
                                fit_event_bias,
                                data = event_data_sample,
                                mc.cores = n_workers)
}

names(event_bias_models) <- event_names

event_bias_models

find_interaction_term <- function(model, term_1, term_2) {
  coefficient_names <- names(fixef(model))
  
  possible_names <- c(paste0(term_1, ":", term_2), paste0(term_2, ":", term_1))
  
  result <- intersect(possible_names, coefficient_names)
  
  if (length(result) != 1) {
    stop("Could not uniquely identify interaction: ",
         term_1,
         " × ",
         term_2)
  }
  
  result
}

extract_propensity_contrast <- function(model,
                                        event_name,
                                        effect_name,
                                        main_term,
                                        interaction_term) {
  beta <- fixef(model)
  V <- as.matrix(vcov(model))
  
  # Non-event effect: main coefficient
  L_non_event <- setNames(rep(0, length(beta)), names(beta))
  
  L_non_event[main_term] <- 1
  
  non_event_estimate <- sum(L_non_event * beta)
  
  non_event_se <- sqrt(as.numeric(t(L_non_event) %*%
                                    V %*%
                                    L_non_event))
  
  # Interaction: difference between observed-event
  # and non-event effects
  L_interaction <- setNames(rep(0, length(beta)), names(beta))
  
  L_interaction[interaction_term] <- 1
  
  interaction_estimate <- sum(L_interaction * beta)
  
  interaction_se <- sqrt(as.numeric(t(L_interaction) %*%
                                      V %*%
                                      L_interaction))
  
  # Observed-event effect: main + interaction
  L_observed <- setNames(rep(0, length(beta)), names(beta))
  
  L_observed[c(main_term, interaction_term)] <- 1
  
  observed_estimate <- sum(L_observed * beta)
  
  observed_se <- sqrt(as.numeric(t(L_observed) %*%
                                   V %*%
                                   L_observed))
  
  tibble(
    event = event_name,
    effect = effect_name,
    
    non_event_estimate = non_event_estimate,
    non_event_std_error = non_event_se,
    non_event_odds_ratio = exp(non_event_estimate),
    non_event_lower_95 = exp(non_event_estimate -
                               1.96 * non_event_se),
    non_event_upper_95 = exp(non_event_estimate +
                               1.96 * non_event_se),
    
    interaction_estimate = interaction_estimate,
    interaction_std_error = interaction_se,
    interaction_odds_ratio = exp(interaction_estimate),
    interaction_lower_95 = exp(interaction_estimate -
                                 1.96 * interaction_se),
    interaction_upper_95 = exp(interaction_estimate +
                                 1.96 * interaction_se),
    
    observed_event_estimate = observed_estimate,
    observed_event_std_error = observed_se,
    observed_event_odds_ratio = exp(observed_estimate),
    observed_event_lower_95 = exp(observed_estimate -
                                    1.96 * observed_se),
    observed_event_upper_95 = exp(observed_estimate +
                                    1.96 * observed_se)
  )
}

table8_effects <- imap_dfr(event_bias_models, function(model, event_name) {
  lead_interaction <- find_interaction_term(model, "observed", "lead_z")
  
  vintage_interaction <- find_interaction_term(model, "observed", "sequence_z")
  
  amendment_interaction <- find_interaction_term(model, "observed", "is_amdamended")
  
  bind_rows(
    extract_propensity_contrast(
      model = model,
      event_name = event_name,
      effect_name = "Lead time",
      main_term = "lead_z",
      interaction_term = lead_interaction
    ),
    
    extract_propensity_contrast(
      model = model,
      event_name = event_name,
      effect_name = "Forecast vintage",
      main_term = "sequence_z",
      interaction_term = vintage_interaction
    ),
    
    extract_propensity_contrast(
      model = model,
      event_name = event_name,
      effect_name = "Amended TAF",
      main_term = "is_amdamended",
      interaction_term = amendment_interaction
    )
  )
})

table8_effects |>
  print(n = Inf, width = Inf)

table8_output <- table8_effects |>
  mutate(
    event = recode(
      event,
      fog = "Fog",
      thunderstorm = "Thunderstorm",
      low_visibility = "Low visibility",
      low_cloud = "Low cloud"
    ),
    
    `Non-event OR (95% interval)` = sprintf(
      "%.3f (%.3f–%.3f)",
      non_event_odds_ratio,
      non_event_lower_95,
      non_event_upper_95
    ),
    
    `Interaction OR (95% interval)` = sprintf(
      "%.3f (%.3f–%.3f)",
      interaction_odds_ratio,
      interaction_lower_95,
      interaction_upper_95
    ),
    
    `Observed-event OR (95% interval)` = sprintf(
      "%.3f (%.3f–%.3f)",
      observed_event_odds_ratio,
      observed_event_lower_95,
      observed_event_upper_95
    )
  ) |>
  select(
    Event = event,
    Effect = effect,
    `Non-event OR (95% interval)`,
    `Interaction OR (95% interval)`,
    `Observed-event OR (95% interval)`
  )

table8_output |>
  print(n = Inf, width = Inf)

table8_random_effects <- event_noise_comparison |>
  filter(model == "forecast_propensity") |>
  mutate(
    event = recode(
      event,
      fog = "Fog",
      thunderstorm = "Thunderstorm",
      low_visibility = "Low visibility",
      low_cloud = "Low cloud"
    )
  ) |>
  select(
    Event = event,
    `Forecaster variance` =
      forecaster_variance,
    `Forecaster SD` =
      forecaster_sd,
    `Forecaster MOR` =
      median_odds_ratio
  )

print(table8_random_effects)
as.data.frame(table8_random_effects) |>
  print(row.names = FALSE)

table8_effects |>
  print(n = Inf, width = Inf)

# This model separates:
# forecasting propensity when no event occurs—false-alarm tendency;
# response to an observed event—detection sensitivity;
# whether sensitivity changes with lead time;
# whether later vintages improve discrimination; and
# persistent forecaster differences in event-forecast propensity.

# Compare forecaster variance across event models
extract_forecaster_component <- function(model, event_name, model_type) {
  if (is.null(model))
    return(NULL)
  
  vc <- as.data.table(as.data.frame(VarCorr(model)))
  
  x <- vc[grp == "ForecasterID"]
  
  if (nrow(x) == 0)
    return(NULL)
  
  data.table(
    event = event_name,
    model = model_type,
    forecaster_variance = x$vcov,
    forecaster_sd = x$sdcor,
    median_odds_ratio =
      exp(0.6745 *
            sqrt(2 * x$vcov))
  )
}

# Combine:
event_noise_comparison <- rbindlist(list(
  rbindlist(
    Map(
      extract_forecaster_component,
      event_accuracy_models,
      names(event_accuracy_models),
      MoreArgs = list(model_type = "overall_accuracy")
    ),
    fill = TRUE
  ),
  
  rbindlist(
    Map(
      extract_forecaster_component,
      event_miss_models,
      names(event_miss_models),
      MoreArgs = list(model_type = "miss_risk")
    ),
    fill = TRUE
  ),
  
  rbindlist(
    Map(
      extract_forecaster_component,
      event_false_alarm_models,
      names(event_false_alarm_models),
      MoreArgs = list(model_type = "false_alarm_risk")
    ),
    fill = TRUE
  ),
  
  rbindlist(
    Map(
      extract_forecaster_component,
      event_bias_models,
      names(event_bias_models),
      MoreArgs = list(model_type = "forecast_propensity")
    ),
    fill = TRUE
  )
), fill = TRUE)

event_noise_comparison[order(event, model)]

# 14. Plot the event-specific forecaster effects without ranking
# Use caterpillar plots ordered only for display and retain anonymity.
extract_forecaster_effects <- function(model, event_name) {
  re <- ranef(model, condVar = TRUE)$ForecasterID
  
  out <- data.table(ForecasterID = rownames(re), estimate = re[, 1])
  
  posterior_var <- attr(re, "postVar")
  
  if (!is.null(posterior_var)) {
    out[, se :=
          sqrt(posterior_var[1, 1, seq_len(dim(posterior_var)[3])])]
    
    out[, `:=`(lower_95 =
                 estimate - 1.96 * se,
               
               upper_95 =
                 estimate + 1.96 * se)]
  }
  
  out[, event := event_name]
  
  out[]
}

event_forecaster_effects <- rbindlist(Map(
  extract_forecaster_effects,
  event_accuracy_models,
  names(event_accuracy_models)
),
fill = TRUE)

#####
names(event_bias_models)        <- event_names
names(event_accuracy_models)    <- event_names
names(event_miss_models)        <- event_names
names(event_false_alarm_models) <- event_names

to_binary <- function(x, variable_name) {
  if (is.logical(x)) {
    return(as.integer(x))
  }
  
  if (is.numeric(x) || is.integer(x)) {
    out <- as.integer(x)
    
  } else {
    x_char <- trimws(tolower(as.character(x)))
    
    out <- dplyr::case_when(
      x_char %in% c(
        "1",
        "yes",
        "true",
        "present",
        "forecast",
        "observed",
        "event"
      ) ~ 1L,
      
      x_char %in% c(
        "0",
        "no",
        "false",
        "absent",
        "not forecast",
        "not observed",
        "non-event",
        "none"
      ) ~ 0L,
      
      TRUE ~ NA_integer_
    )
  }
  
  valid_values <- unique(out[!is.na(out)])
  
  if (!all(valid_values %in% c(0L, 1L))) {
    stop(variable_name, " could not be converted cleanly to 0/1.")
  }
  
  out
}

hazard_sample_flow <- imap_dfr(event_bias_models, function(model, event_name) {
  d <- model.frame(model)
  
  required_variables <- c("forecast",
                          "observed",
                          "ForecasterID",
                          "aerodrome",
                          "weather_block")
  
  missing_variables <- setdiff(required_variables, names(d))
  
  if (length(missing_variables) > 0) {
    stop(
      "Missing variable(s) in the ",
      event_name,
      " model frame: ",
      paste(missing_variables, collapse = ", ")
    )
  }
  
  forecast_binary <- to_binary(d$forecast, "forecast")
  
  observed_binary <- to_binary(d$observed, "observed")
  
  complete_binary <- !is.na(forecast_binary) &
    !is.na(observed_binary)
  
  forecast_binary <- forecast_binary[complete_binary]
  
  observed_binary <- observed_binary[complete_binary]
  
  hits <- sum(forecast_binary == 1L &
                observed_binary == 1L)
  
  misses <- sum(forecast_binary == 0L &
                  observed_binary == 1L)
  
  false_alarms <- sum(forecast_binary == 1L &
                        observed_binary == 0L)
  
  correct_negatives <- sum(forecast_binary == 0L &
                             observed_binary == 0L)
  
  observed_event_hours <- hits + misses
  non_event_hours <- false_alarms +
    correct_negatives
  
  total_records <- observed_event_hours +
    non_event_hours
  
  tibble(
    event = event_name,
    
    total_records = total_records,
    
    observed_event_hours =
      observed_event_hours,
    
    non_event_hours =
      non_event_hours,
    
    hits = hits,
    misses = misses,
    false_alarms = false_alarms,
    
    correct_negatives =
      correct_negatives,
    
    event_prevalence =
      observed_event_hours /
      total_records,
    
    forecast_prevalence =
      (hits + false_alarms) /
      total_records,
    
    probability_of_detection =
      if_else(
        observed_event_hours > 0,
        hits / observed_event_hours,
        NA_real_
      ),
    
    false_alarm_ratio =
      if_else(
        hits + false_alarms > 0,
        false_alarms /
          (hits + false_alarms),
        NA_real_
      ),
    
    critical_success_index =
      if_else(
        hits + misses + false_alarms > 0,
        hits /
          (hits + misses + false_alarms),
        NA_real_
      ),
    
    forecasters =
      n_distinct(d$ForecasterID[complete_binary]),
    
    aerodromes =
      n_distinct(d$aerodrome[complete_binary]),
    
    weather_blocks =
      n_distinct(d$weather_block[complete_binary])
  )
})

hazard_sample_flow |>
  print(n = Inf, width = Inf)

hazard_sample_flow_output <- hazard_sample_flow |>
  mutate(
    event = recode(
      event,
      fog = "Fog",
      thunderstorm = "Thunderstorm",
      low_visibility = "Low visibility",
      low_cloud = "Low cloud"
    ),
    
    across(
      c(
        total_records,
        observed_event_hours,
        non_event_hours,
        hits,
        misses,
        false_alarms,
        correct_negatives
      ),
      ~ format(
        .x,
        big.mark = ",",
        scientific = FALSE,
        trim = TRUE
      )
    ),
    
    event_prevalence = sprintf("%.2f%%", 100 * event_prevalence),
    
    forecast_prevalence = sprintf("%.2f%%", 100 * forecast_prevalence),
    
    probability_of_detection = sprintf("%.3f", probability_of_detection),
    
    false_alarm_ratio = sprintf("%.3f", false_alarm_ratio),
    
    critical_success_index = sprintf("%.3f", critical_success_index)
  ) |>
  select(
    Event = event,
    `Total records` = total_records,
    `Observed-event hours` =
      observed_event_hours,
    `Non-event hours` =
      non_event_hours,
    Hits = hits,
    Misses = misses,
    `False alarms` = false_alarms,
    `Correct negatives` =
      correct_negatives,
    `Event prevalence` =
      event_prevalence,
    `Forecast prevalence` =
      forecast_prevalence,
    POD = probability_of_detection,
    FAR = false_alarm_ratio,
    CSI = critical_success_index,
    Forecasters = forecasters,
    Aerodromes = aerodromes,
    `Weather blocks` = weather_blocks
  )

hazard_sample_flow_output |>
  print(n = Inf, width = Inf)

get_model_n <- function(model_list, event_name) {
  model <- model_list[[event_name]]
  
  if (is.null(model)) {
    return(NA_integer_)
  }
  
  as.integer(nobs(model))
}

model_sample_reconciliation <- tibble(event = event_names) |>
  mutate(
    accuracy_model_n = map_int(event, ~ get_model_n(event_accuracy_models, .x)),
    
    miss_model_n = map_int(event, ~ get_model_n(event_miss_models, .x)),
    
    false_alarm_model_n = map_int(event, ~ get_model_n(event_false_alarm_models, .x)),
    
    propensity_model_n = map_int(event, ~ get_model_n(event_bias_models, .x))
  ) |>
  left_join(
    hazard_sample_flow |>
      select(event, total_records, observed_event_hours, non_event_hours),
    by = "event"
  ) |>
  mutate(
    accuracy_minus_flow =
      accuracy_model_n - total_records,
    
    propensity_minus_flow =
      propensity_model_n - total_records,
    
    miss_minus_observed =
      miss_model_n - observed_event_hours,
    
    false_alarm_minus_non_event =
      false_alarm_model_n - non_event_hours,
    
    accuracy_population_matches =
      accuracy_minus_flow == 0,
    
    propensity_population_matches =
      propensity_minus_flow == 0,
    
    miss_population_matches =
      miss_minus_observed == 0,
    
    false_alarm_population_matches =
      false_alarm_minus_non_event == 0
  )

model_sample_reconciliation |>
  print(n = Inf, width = Inf)

model_reconciliation_output <-
  model_sample_reconciliation |>
  mutate(
    event = recode(
      event,
      fog = "Fog",
      thunderstorm = "Thunderstorm",
      low_visibility = "Low visibility",
      low_cloud = "Low cloud"
    )
  ) |>
  select(
    Event = event,
    
    `Contingency total` =
      total_records,
    
    `Accuracy model` =
      accuracy_model_n,
    
    `Accuracy difference` =
      accuracy_minus_flow,
    
    `Propensity model` =
      propensity_model_n,
    
    `Propensity difference` =
      propensity_minus_flow,
    
    `Observed-event hours` =
      observed_event_hours,
    
    `Miss model` =
      miss_model_n,
    
    `Miss difference` =
      miss_minus_observed,
    
    `Non-event hours` =
      non_event_hours,
    
    `False-alarm model` =
      false_alarm_model_n,
    
    `False-alarm difference` =
      false_alarm_minus_non_event
  )

model_reconciliation_output |>
  print(n = Inf, width = Inf)

names(event_accuracy_models)    <- event_names
names(event_miss_models)        <- event_names
names(event_false_alarm_models) <- event_names
names(event_bias_models)        <- event_names

diagnose_propensity_exclusions <- purrr::imap_dfr(event_accuracy_models, function(model, event_name) {
  d <- model.frame(model)
  
  tibble::tibble(
    event = event_name,
    accuracy_records = nrow(d),
    forecast_type_missing =
      if ("forecast_type" %in% names(d)) {
        sum(is.na(d$forecast_type))
      } else {
        NA_integer_
      },
    forecast_type_missing_observed =
      if (all(c("forecast_type", "observed") %in% names(d))) {
        sum(is.na(d$forecast_type) &
              d$observed == 1, na.rm = TRUE)
      } else {
        NA_integer_
      },
    forecast_type_missing_non_event =
      if (all(c("forecast_type", "observed") %in% names(d))) {
        sum(is.na(d$forecast_type) &
              d$observed == 0, na.rm = TRUE)
      } else {
        NA_integer_
      },
    excluded_from_propensity =
      nrow(d) -
      nobs(event_bias_models[[event_name]])
  )
})

diagnose_propensity_exclusions |>
  print(n = Inf, width = Inf)

purrr::imap_dfr(event_bias_models, function(model, event_name) {
  d <- model.frame(model, na.action = NULL)
  
  tibble::tibble(
    event = event_name,
    records = nrow(d),
    missing_forecast_type =
      sum(is.na(d$forecast_type)),
    missing_forecast =
      sum(is.na(d$forecast)),
    missing_observed =
      sum(is.na(d$observed))
  )
}) |>
  print(n = Inf, width = Inf)

sample_flow_from_accuracy <- purrr::imap_dfr(event_accuracy_models, function(model, event_name) {
  d <- model.frame(model)
  
  tibble::tibble(
    event = event_name,
    total_records = nobs(model),
    observed_event_hours =
      nobs(event_miss_models[[event_name]]),
    non_event_hours =
      nobs(event_false_alarm_models[[event_name]]),
    forecasters =
      dplyr::n_distinct(d$ForecasterID),
    aerodromes =
      dplyr::n_distinct(d$aerodrome),
    weather_blocks =
      dplyr::n_distinct(d$weather_block)
  )
})

sample_flow_from_accuracy |>
  print(n = Inf, width = Inf)

purrr::imap(event_accuracy_models, ~ names(model.frame(.x)))

# Plot
library(ggplot2)

ggplot(event_forecaster_effects,
       aes(
         x = reorder(interaction(event, ForecasterID), estimate),
         y = estimate,
         colour = event
       )) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_pointrange(aes(ymin = lower_95, ymax = upper_95), alpha = 0.7) +
  facet_wrap(~ event, scales = "free_x") +
  coord_flip() +
  labs(x = "Anonymised forecaster", y = "Adjusted forecaster effect on accuracy log-odds", title = "Partially pooled forecaster effects by operational event") +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none"
  )

# Plot 2
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)

# ------------------------------------------------------------
# 1. Reshape the false-alarm results
# ------------------------------------------------------------

false_alarm_plot_data <- false_alarm_model_table |>
  transmute(
    event,
    `Lead time` = lead_odds_ratio,
    `Lead time lower` = lead_lower_95,
    `Lead time upper` = lead_upper_95,
    
    `Forecast vintage` = vintage_odds_ratio,
    `Forecast vintage lower` = vintage_lower_95,
    `Forecast vintage upper` = vintage_upper_95,
    
    `Amended TAF` = amendment_odds_ratio,
    `Amended TAF lower` = amendment_lower_95,
    `Amended TAF upper` = amendment_upper_95
  ) |>
  pivot_longer(
    cols = -event,
    names_to = "variable",
    values_to = "value"
  ) |>
  mutate(
    bound = case_when(
      grepl(" lower$", variable) ~ "lower_95",
      grepl(" upper$", variable) ~ "upper_95",
      TRUE ~ "odds_ratio"
    ),
    effect = variable |>
      sub(" lower$", "", x = _) |>
      sub(" upper$", "", x = _)
  ) |>
  select(
    event,
    effect,
    bound,
    value
  ) |>
  pivot_wider(
    names_from = bound,
    values_from = value
  ) |>
  mutate(
    error_type = "False alarm"
  )

false_alarm_plot_data <- false_alarm_model_table |>
  transmute(
    event,
    lead_odds_ratio,
    lead_lower_95,
    lead_upper_95,
    vintage_odds_ratio,
    vintage_lower_95,
    vintage_upper_95,
    amendment_odds_ratio,
    amendment_lower_95,
    amendment_upper_95
  ) |>
  pivot_longer(
    cols = -event,
    names_to = c("effect", ".value"),
    names_pattern =
      "(lead|vintage|amendment)_(odds_ratio|lower_95|upper_95)"
  ) |>
  mutate(
    effect = recode(
      effect,
      lead = "Lead time",
      vintage = "Forecast vintage",
      amendment = "Amended TAF"
    ),
    error_type = "False alarm"
  )

miss_plot_data <- miss_model_table |>
  transmute(
    event,
    effect,
    odds_ratio,
    lower_95,
    upper_95,
    error_type = "Miss"
  )

paired_error_plot_data <- bind_rows(
  miss_plot_data,
  false_alarm_plot_data
) |>
  mutate(
    event = recode(
      event,
      fog = "Fog",
      thunderstorm = "Thunderstorm",
      low_visibility = "Low visibility",
      low_cloud = "Low cloud"
    ),
    
    event = factor(
      event,
      levels = c(
        "Low cloud",
        "Low visibility",
        "Thunderstorm",
        "Fog"
      )
    ),
    
    effect = factor(
      effect,
      levels = c(
        "Lead time",
        "Forecast vintage",
        "Amended TAF"
      )
    ),
    
    error_type = factor(
      error_type,
      levels = c(
        "Miss",
        "False alarm"
      )
    )
  )

paired_error_forest_plot <- ggplot(
  paired_error_plot_data,
  aes(
    x = odds_ratio,
    y = event,
    colour = error_type,
    shape = error_type
  )
) +
  geom_vline(
    xintercept = 1,
    linewidth = 0.45,
    colour = "grey40",
    linetype = "dashed"
  ) +
  geom_errorbarh(
    aes(
      xmin = lower_95,
      xmax = upper_95
    ),
    height = 0,
    linewidth = 0.65,
    position = position_dodge(
      width = 0.45
    )
  ) +
  geom_point(
    size = 2.5,
    position = position_dodge(
      width = 0.45
    )
  ) +
  facet_wrap(
    ~ effect,
    nrow = 1,
    scales = "free_x"
  ) +
  scale_x_log10() +
  scale_colour_manual(
    values = c(
      "Miss" = "#0072B2",
      "False alarm" = "#D55E00"
    )
  ) +
  scale_shape_manual(
    values = c(
      "Miss" = 16,
      "False alarm" = 17
    )
  ) +
  labs(
    x = "Adjusted odds ratio (log scale)",
    y = NULL,
    colour = "Conditional error",
    shape = "Conditional error"
  ) +
  theme_bw(
    base_size = 11
  ) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    strip.background = element_rect(
      fill = "grey92",
      colour = "grey60"
    ),
    strip.text = element_text(
      face = "bold"
    ),
    axis.text.y = element_text(
      colour = "black"
    ),
    plot.margin = margin(
      6, 8, 6, 6
    )
  )

paired_error_forest_plot

plot_lower <- 0.10
plot_upper <- 4.00

paired_error_plot_display <- paired_error_plot_data |>
  mutate(
    display_or = pmin(
      pmax(odds_ratio, plot_lower),
      plot_upper
    ),
    display_lower = pmax(
      lower_95,
      plot_lower
    ),
    display_upper = pmin(
      upper_95,
      plot_upper
    ),
    lower_truncated =
      lower_95 < plot_lower,
    upper_truncated =
      upper_95 > plot_upper
  )

paired_error_forest_plot <- ggplot(
  paired_error_plot_display,
  aes(
    x = display_or,
    y = event,
    colour = error_type,
    shape = error_type
  )
) +
  geom_vline(
    xintercept = 1,
    linewidth = 0.45,
    colour = "grey40",
    linetype = "dashed"
  ) +
  geom_segment(
    aes(
      x = display_lower,
      xend = display_upper,
      yend = event
    ),
    linewidth = 0.65,
    position = position_dodge(
      width = 0.45
    )
  ) +
  geom_point(
    size = 2.5,
    position = position_dodge(
      width = 0.45
    )
  ) +
  geom_point(
    data = \(x) filter(
      x,
      lower_truncated
    ),
    aes(
      x = plot_lower
    ),
    shape = 60,
    size = 3,
    position = position_dodge(
      width = 0.45
    ),
    show.legend = FALSE
  ) +
  geom_point(
    data = \(x) filter(
      x,
      upper_truncated
    ),
    aes(
      x = plot_upper
    ),
    shape = 62,
    size = 3,
    position = position_dodge(
      width = 0.45
    ),
    show.legend = FALSE
  ) +
  facet_wrap(
    ~ effect,
    nrow = 1
  ) +
  scale_x_log10(
    limits = c(
      plot_lower,
      plot_upper
    ),
    breaks = c(
      0.1,
      0.25,
      0.5,
      1,
      2,
      4
    )
  ) +
  scale_colour_manual(
    values = c(
      "Miss" = "#0072B2",
      "False alarm" = "#D55E00"
    )
  ) +
  scale_shape_manual(
    values = c(
      "Miss" = 16,
      "False alarm" = 17
    )
  ) +
  labs(
    x = "Adjusted odds ratio (log scale)",
    y = NULL,
    colour = NULL,
    shape = NULL
  ) +
  theme_bw(
    base_size = 11
  ) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    strip.background = element_rect(
      fill = "grey92",
      colour = "grey60"
    ),
    strip.text = element_text(
      face = "bold"
    )
  )

paired_error_forest_plot

#### END ########################
