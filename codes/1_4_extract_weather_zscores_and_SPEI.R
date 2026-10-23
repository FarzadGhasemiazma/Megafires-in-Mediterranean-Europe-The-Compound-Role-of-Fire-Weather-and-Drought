#!/usr/bin/env Rscript

## ============================================================
## Extract daily JJAS weather z-scores and monthly SPEI for each fire
## ------------------------------------------------------------
## Fires:
##   /diskonfire/EFFIS_FIRE_DATA/MESOGEOS_DOMAIN/
##      all_fires_mesogeos_JJAS_gt30ha_2008_2022_med_countries.gpkg
##
## Daily JJAS z-score variables (lags -1, 0, +1):
##   lai, smi, lst_day, lst_night, ndvi, rh, t2m, wind_speed, tp
##   Files: <var>_JJAS_zscore_2008_2022.nc in /diskonfire/MESOGEOS
##
## Monthly SPEI variables (time-interpolated, no lag):
##   spei01, spei03, spei06
##   Files: /diskonfire/SPEI_MSWEP_hPET/spei[01|03|06]_mesogeos.nc
##
## For each fire:
##   * Start date d (column `initialdat`)
##   * Daily vars at d-1, d, d+1:
##       - polygon mean; if all NA -> use 5 km, then 10 km, then 20 km buffer
##   * SPEI linearly interpolated in time between previous and current month:
##       - Let d      = fire start date
##       - Let m      = month(d)
##       - Let d_prev = first day of previous month (m-1)
##       - Let d_curr = first day of current month (m)
##       - Let d_next = first day of next month (m+1)
##       - Let SPEI_prev = SPEI at d_prev
##       - Let SPEI_curr = SPEI at d_curr
##       - Let month_len = d_next - d_curr  (number of days in month m)
##       - Let offset    = d - d_curr      (days since start of month m, 0-based)
##         Then:
##           w_curr = offset / month_len
##           w_prev = 1 - w_curr
##           SPEI_interp(d) = w_prev * SPEI_prev + w_curr * SPEI_curr
##
##       - Early in the month (d close to 1)  -> weight on SPEI_prev (month m-1)
##       - Late in the month (d close to 30) -> weight on SPEI_curr (month m)
##
## Output:
##   /diskonfire/MESOGEOS/fire_meteo_JJAS_zscore_SPEI_2008_2022.csv
## ============================================================

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(lubridate)
})

terraOptions(progress = 1)

## -------------------------------
## 1. PATHS AND BASIC SETTINGS
## -------------------------------

fire_gpkg <- "/diskonfire/EFFIS_FIRE_DATA/MESOGEOS_DOMAIN/all_fires_mesogeos_JJAS_gt30ha_2008_2022_med_countries.gpkg"
meteo_dir <- "/diskonfire/MESOGEOS"

## Daily JJAS z-score variables
zvars <- c(
  "lai",
  "smi",
  "lst_day",
  "lst_night",
  "ndvi",
  "rh",
  "t2m",
  "wind_speed",
  "tp"
)

zscore_files <- setNames(
  file.path(meteo_dir, paste0(zvars, "_JJAS_zscore_2008_2022.nc")),
  zvars
)

## Monthly SPEI variables
spei_scales <- c("spei01", "spei03", "spei06")
spei_dir    <- "/diskonfire/SPEI_MSWEP_hPET"
spei_files  <- setNames(
  file.path(spei_dir, paste0(spei_scales, "_mesogeos.nc")),
  spei_scales
)

## Fire attribute columns
ID_COL   <- "id"
DATE_COL <- "initialdat"

## CRS for metric buffers (Europe LAEA)
CRS_LAEA <- "EPSG:3035"

## Daily lags (in days)
lags      <- c(-1, 0, 1)
lag_names <- c("m1", "d0", "p1")

## Output
out_file <- "/diskonfire/MESOGEOS/fire_meteo_JJAS_zscore_SPEI_2008_2022.csv"

## ============================================================
## 2. READ FIRES AND FIX DATE COLUMN
## ============================================================

cat("Reading fires from:\n  ", fire_gpkg, "\n")
fires_raw <- terra::vect(fire_gpkg)
cat("Number of fires read:", nrow(fires_raw), "\n")

if (!DATE_COL %in% names(fires_raw)) {
  stop("DATE_COL '", DATE_COL, "' not found in fire layer.")
}
if (!ID_COL %in% names(fires_raw)) {
  stop("ID_COL '", ID_COL, "' not found in fire layer.")
}

fires_df <- as.data.frame(fires_raw)
date_raw <- fires_df[[DATE_COL]]

## Robust conversion to Date
if (inherits(date_raw, "Date")) {
  fires_df[[DATE_COL]] <- date_raw
} else if (inherits(date_raw, "POSIXt")) {
  fires_df[[DATE_COL]] <- as.Date(date_raw)
} else if (is.numeric(date_raw)) {
  ## Days since 1970-01-01 (e.g. 14037, 19249, ...)
  fires_df[[DATE_COL]] <- as.Date(date_raw, origin = "1970-01-01")
} else {
  fires_df[[DATE_COL]] <- as.Date(date_raw)
}

cat("Class of date column:", class(fires_df[[DATE_COL]]), "\n")
cat(
  "Date range (fires): ",
  format(min(fires_df[[DATE_COL]], na.rm = TRUE)), " -> ",
  format(max(fires_df[[DATE_COL]], na.rm = TRUE)), "\n", sep = ""
)

## Back to SpatVector and filter valid dates
fires <- fires_raw
fires[[DATE_COL]] <- fires_df[[DATE_COL]]

valid   <- !is.na(fires[[DATE_COL]])
fires   <- fires[valid, ]
fires_df <- fires_df[valid, ]

cat("Fires with valid date:", nrow(fires), "\n")

## ============================================================
## 3. REPROJECT FIRES AND BUILD BUFFERS (FOR DAILY METEO & SPEI)
## ============================================================

## Reference CRS from one z-score file (t2m)
ref_file <- zscore_files[["t2m"]]
cat("Using reference CRS from file:\n  ", ref_file, "\n")
r_ref   <- terra::rast(ref_file)
crs_ref <- terra::crs(r_ref)
cat("Reference CRS:\n", crs_ref, "\n")

## Fires to reference CRS (usually lon/lat)
fires_ll <- terra::project(fires, crs_ref)

cat("Building 5, 10 and 20 km buffers...\n")
fires_laea       <- terra::project(fires_ll, CRS_LAEA)
fires_buf0_laea  <- fires_laea
fires_buf5_laea  <- terra::buffer(fires_laea,  5000)   # 5 km
fires_buf10_laea <- terra::buffer(fires_laea, 10000)   # 10 km
fires_buf20_laea <- terra::buffer(fires_laea, 20000)   # 20 km

## Back to reference CRS (for all rasters we will use)
fires_buf0  <- terra::project(fires_buf0_laea,  crs_ref)
fires_buf5  <- terra::project(fires_buf5_laea,  crs_ref)
fires_buf10 <- terra::project(fires_buf10_laea, crs_ref)
fires_buf20 <- terra::project(fires_buf20_laea, crs_ref)

## ============================================================
## 4. LOAD DAILY JJAS Z-SCORE RASTERS
## ============================================================

cat("Loading JJAS z-score SpatRasters:\n")

zrasters <- list()
ztimes   <- list()

for (v in zvars) {
  f <- zscore_files[[v]]
  if (!file.exists(f)) stop("Missing NetCDF for '", v, "': ", f)
  cat("  -", v, ":", f, "\n")
  
  r <- terra::rast(f)
  zrasters[[v]] <- r
  
  tt <- terra::time(r)
  
  ## Robust conversion of time to Date
  if (inherits(tt, "Date")) {
    zt <- tt
  } else if (inherits(tt, "POSIXt")) {
    zt <- as.Date(tt)
  } else if (is.numeric(tt)) {
    ## If terra already applied the origin, as.Date(tt) will be fine.
    ## If not, adjust origin here. (For now we assume terra did the right thing.)
    zt <- as.Date(tt)
  } else {
    stop("Unrecognized time type for variable ", v, ": ", paste(class(tt), collapse = ", "))
  }
  
  ztimes[[v]] <- zt
  
  if (v == "t2m") {
    cat("     class(terra::time) =", paste(class(tt), collapse = ", "), "\n")
    cat("     t2m date range: ",
        format(min(zt)), " -> ", format(max(zt)), "\n", sep = "")
    cat("     first 5 dates: ",
        paste(format(head(zt, 5)), collapse = ", "), "\n", sep = "")
    cat("     last 5 dates : ",
        paste(format(tail(zt, 5)), collapse = ", "), "\n", sep = "")
  }
}

## ============================================================
## 5. LOAD MONTHLY SPEI RASTERS
## ============================================================

cat("Loading monthly SPEI rasters (if available):\n")

spei_rasters <- list()
spei_dates   <- list()
loaded_spei  <- character(0)

for (nm in spei_scales) {
  f <- spei_files[[nm]]
  if (!file.exists(f)) {
    cat("  -", nm, ": file NOT found ->", f, "\n")
    next
  }
  cat("  -", nm, ":", f, "\n")
  r <- terra::rast(f)
  tt <- terra::time(r)
  dtt <- as.Date(tt)  # terra should handle the origin from NetCDF
  
  spei_rasters[[nm]] <- r
  spei_dates[[nm]]   <- dtt
  loaded_spei        <- c(loaded_spei, nm)
  
  cat("      ", nm, "date range: ",
      format(min(dtt)), " -> ", format(max(dtt)), "\n", sep = "")
}

if (length(loaded_spei) == 0) {
  cat("No SPEI rasters loaded. SPEI columns will remain NA.\n")
}
spei_scales <- loaded_spei  # keep only those that actually exist

## ============================================================
## 6. GENERIC EXTRACTION FUNCTION (WITH BUFFERS)
## ============================================================

extract_with_buffers <- function(r_layer,
                                 fires0, fires5, fires10, fires20,
                                 idx_fires) {
  ## r_layer: single-layer SpatRaster
  ## fires*:  SpatVector (same CRS as r_layer)
  ## idx_fires: indices of fires to process (row indices in fires*)
  
  vals <- rep(NA_real_, length(idx_fires))
  
  do_extract <- function(r, v) {
    if (nrow(v) == 0) return(rep(NA_real_, length(idx_fires)))
    e <- terra::extract(r, v, fun = mean, na.rm = TRUE)
    as.numeric(e[[2]])
  }
  
  ## 1) Original polygons
  v0 <- do_extract(r_layer, fires0[idx_fires, ])
  vals[] <- v0
  missing <- which(is.na(vals))
  
  ## 2) 5 km buffer
  if (length(missing) > 0) {
    idx5 <- idx_fires[missing]
    v5   <- do_extract(r_layer, fires5[idx5, ])
    ok   <- which(!is.na(v5))
    if (length(ok) > 0) vals[missing[ok]] <- v5[ok]
    missing <- which(is.na(vals))
  }
  
  ## 3) 10 km buffer
  if (length(missing) > 0) {
    idx10 <- idx_fires[missing]
    v10   <- do_extract(r_layer, fires10[idx10, ])
    ok    <- which(!is.na(v10))
    if (length(ok) > 0) vals[missing[ok]] <- v10[ok]
    missing <- which(is.na(vals))
  }
  
  ## 4) 20 km buffer
  if (length(missing) > 0) {
    idx20 <- idx_fires[missing]
    v20   <- do_extract(r_layer, fires20[idx20, ])
    ok    <- which(!is.na(v20))
    if (length(ok) > 0) vals[missing[ok]] <- v20[ok]
  }
  
  vals
}

## ============================================================
## 7. FUNCTION: SPEI PREVIOUS–CURRENT MONTH & WEIGHTS
## ============================================================

spei_prev_curr_for_fire <- function(d_fire, spei_dates_vec) {
  ## Convert fire date to Date if needed
  if (inherits(d_fire, "Date")) {
    d0 <- d_fire
  } else if (inherits(d_fire, "POSIXt")) {
    d0 <- as.Date(d_fire)
  } else if (is.numeric(d_fire)) {
    d0 <- as.Date(d_fire, origin = "1970-01-01")
  } else {
    d0 <- as.Date(d_fire)
  }
  
  if (is.na(d0)) return(NULL)
  
  y <- as.integer(format(d0, "%Y"))
  m <- as.integer(format(d0, "%m"))
  
  ## Current month start (d_curr)
  d_curr <- as.Date(sprintf("%04d-%02d-01", y, m))
  
  ## Previous month start (d_prev)
  if (m == 1) {
    y_prev <- y - 1L
    m_prev <- 12L
  } else {
    y_prev <- y
    m_prev <- m - 1L
  }
  d_prev <- as.Date(sprintf("%04d-%02d-01", y_prev, m_prev))
  
  ## Next month start (d_next) – just to compute month length
  if (m == 12) {
    y_next <- y + 1L
    m_next <- 1L
  } else {
    y_next <- y
    m_next <- m + 1L
  }
  d_next <- as.Date(sprintf("%04d-%02d-01", y_next, m_next))
  
  ## Both previous and current month must exist in this SPEI series
  if (!(d_prev %in% spei_dates_vec) || !(d_curr %in% spei_dates_vec)) {
    return(NULL)
  }
  
  ## Length of month m (in days) and offset from start of month
  month_len <- as.numeric(d_next - d_curr)
  if (month_len <= 0) return(NULL)
  
  offset <- as.numeric(d0 - d_curr)
  ## Clamp offset to [0, month_len] just in case
  offset <- max(0, min(offset, month_len))
  
  ## Weights: early in month → previous month dominates; late in month → current month dominates
  w_curr <- offset / month_len
  w_prev <- 1 - w_curr
  
  list(
    d_prev = d_prev,
    d_curr = d_curr,
    w_prev = w_prev,
    w_curr = w_curr
  )
}

## ============================================================
## 8. PREPARE OUTPUT DATA FRAME
## ============================================================

n_fire    <- nrow(fires_buf0)
fire_id   <- fires_df[[ID_COL]]
fire_date <- fires_df[[DATE_COL]]

cat("Total number of fires to process:", n_fire, "\n")

out <- data.frame(
  fire_id    = fire_id,
  date_start = fire_date,
  stringsAsFactors = FALSE
)

## Keep coordinates and burned area (if available)
if ("x" %in% names(fires_df))      out$x       <- fires_df[["x"]]
if ("y" %in% names(fires_df))      out$y       <- fires_df[["y"]]
if ("area_ha" %in% names(fires_df)) out$area_ha <- fires_df[["area_ha"]]

## Daily z-score columns
for (v in zvars) {
  for (k in seq_along(lags)) {
    ln  <- lag_names[k]
    col <- paste0(v, "_z_", ln)
    out[[col]] <- NA_real_
  }
}

## SPEI columns (interpolated, no lag)
for (nm in spei_scales) {
  col_spei <- paste0(nm, "_interp")
  out[[col_spei]] <- NA_real_
}

## ============================================================
## 9. MAIN LOOP OVER UNIQUE FIRE START DATES
## ============================================================

cat("Building vector of unique fire start dates...\n")
unique_dates <- sort(unique(as.Date(fire_date)))
cat("Number of unique fire start dates:", length(unique_dates), "\n")

for (d0 in unique_dates) {
  idx_fires <- which(fire_date == d0)
  if (length(idx_fires) == 0) next
  
  cat(
    "\n=== Fire start date:", format(d0),
    " (fires:", length(idx_fires), ") ===\n"
  )
  
  ## ---------------------------------------
  ## 9a. DAILY JJAS Z-SCORES (WITH LAGS)
  ## ---------------------------------------
  for (i_lag in seq_along(lags)) {
    lag      <- lags[i_lag]
    lag_name <- lag_names[i_lag]
    d_target <- d0 + lag
    
    cat("  >> Lag", lag, "-> target date:", format(d_target), "\n")
    
    for (v in zvars) {
      r  <- zrasters[[v]]
      tt <- ztimes[[v]]
      
      ## Robust matching using numeric representation of Dates
      layer_idx <- match(as.numeric(d_target), as.numeric(tt))
      if (is.na(layer_idx)) {
        ## No data for that date (remains NA)
        next
      }
      
      r_layer <- r[[layer_idx]]
      
      vals <- extract_with_buffers(
        r_layer,
        fires_buf0, fires_buf5, fires_buf10, fires_buf20,
        idx_fires
      )
      
      colname <- paste0(v, "_z_", lag_name)
      out[[colname]][idx_fires] <- vals
    }
  }
  
  ## ---------------------------------------
  ## 9b. MONTHLY SPEI (PREVIOUS–CURRENT INTERPOLATION)
  ## ---------------------------------------
  if (length(spei_scales) > 0) {
    for (nm in spei_scales) {
      r_spei <- spei_rasters[[nm]]
      d_spei <- spei_dates[[nm]]
      
      info <- spei_prev_curr_for_fire(d0, d_spei)
      if (is.null(info)) next
      
      i_prev <- match(info$d_prev, d_spei)
      i_curr <- match(info$d_curr, d_spei)
      if (is.na(i_prev) || is.na(i_curr)) next
      
      r_prev <- r_spei[[i_prev]]
      r_curr <- r_spei[[i_curr]]
      
      vals_prev <- extract_with_buffers(
        r_prev,
        fires_buf0, fires_buf5, fires_buf10, fires_buf20,
        idx_fires
      )
      vals_curr <- extract_with_buffers(
        r_curr,
        fires_buf0, fires_buf5, fires_buf10, fires_buf20,
        idx_fires
      )
      
      ## Linear interpolation: early month → previous, late month → current
      vals_interp <- info$w_prev * vals_prev + info$w_curr * vals_curr
      
      col_spei <- paste0(nm, "_interp")
      out[[col_spei]][idx_fires] <- vals_interp
    }
  }
}

## ============================================================
## 10. WRITE OUTPUT
## ============================================================

cat("\nSaving final table to:\n  ", out_file, "\n")
write.csv(out, out_file, row.names = FALSE)
cat("Done.\n")
