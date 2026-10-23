#!/usr/bin/env Rscript

## ============================================================
## Extract daily JJAS LST_night for each fire (lags -2, -1, 0, +1)
## ------------------------------------------------------------
## Fires:
##   /diskonfire/EFFIS_FIRE_DATA/MESOGEOS_DOMAIN/
##      all_fires_mesogeos_JJAS_gt30ha_2008_2022_med_countries.gpkg
##
## Daily JJAS variable:
##   lst_night  (raw, NOT standardized)
##   File: lst_night_JJAS_2008_2022.nc in /diskonfire/MESOGEOS
##
## For each fire:
##   * Start date d (column `initialdat`)
##   * LST_night at d-2, d-1, d, d+1:
##       - polygon mean; if all NA -> use 5 km, then 10 km, then 20 km buffer
##
## Output:
##   /diskonfire/MESOGEOS/fire_LSTnight_JJAS_2008_2022.csv
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

## Single daily JJAS variable (raw, NOT z-score)
var_name   <- "lst_night"
lst_file   <- file.path(meteo_dir, paste0(var_name, "_JJAS_2008_2022.nc"))

## Fire attribute columns
ID_COL   <- "id"
DATE_COL <- "initialdat"

## CRS for metric buffers (Europe LAEA)
CRS_LAEA <- "EPSG:3035"

## Daily lags (in days): -2, -1, 0, +1
lags      <- c(-2, -1, 0, 1)
lag_names <- c("m2", "m1", "d0", "p1")

## Output
out_file <- "/diskonfire/MESOGEOS/fire_LSTnight_JJAS_2008_2022.csv"

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

valid    <- !is.na(fires[[DATE_COL]])
fires    <- fires[valid, ]
fires_df <- fires_df[valid, ]

cat("Fires with valid date:", nrow(fires), "\n")

## ============================================================
## 3. REPROJECT FIRES AND BUILD BUFFERS
## ============================================================

## Reference CRS from LST_night file
if (!file.exists(lst_file)) {
  stop("LST_night NetCDF not found: ", lst_file)
}

cat("Using reference CRS from file:\n  ", lst_file, "\n")
r_ref   <- terra::rast(lst_file)
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
## 4. LOAD DAILY JJAS LST_NIGHT RASTER
## ============================================================

cat("Loading JJAS LST_night SpatRaster:\n")
r_lst <- terra::rast(lst_file)
tt    <- terra::time(r_lst)

## Robust conversion of time to Date
if (inherits(tt, "Date")) {
  lst_dates <- tt
} else if (inherits(tt, "POSIXt")) {
  lst_dates <- as.Date(tt)
} else if (is.numeric(tt)) {
  lst_dates <- as.Date(tt)
} else {
  stop("Unrecognized time type for LST_night: ", paste(class(tt), collapse = ", "))
}

cat("  class(terra::time) =", paste(class(tt), collapse = ", "), "\n")
cat("  LST_night date range: ",
    format(min(lst_dates)), " -> ", format(max(lst_dates)), "\n", sep = "")
cat("  first 5 dates: ",
    paste(format(head(lst_dates, 5)), collapse = ", "), "\n", sep = "")
cat("  last 5 dates : ",
    paste(format(tail(lst_dates, 5)), collapse = ", "), "\n", sep = "")

## ============================================================
## 5. GENERIC EXTRACTION FUNCTION (WITH BUFFERS)
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
## 6. PREPARE OUTPUT DATA FRAME
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
if ("x" %in% names(fires_df))       out$x       <- fires_df[["x"]]
if ("y" %in% names(fires_df))       out$y       <- fires_df[["y"]]
if ("area_ha" %in% names(fires_df)) out$area_ha <- fires_df[["area_ha"]]

## Daily LST_night columns for each lag
for (k in seq_along(lags)) {
  ln  <- lag_names[k]
  col <- paste0(var_name, "_", ln)   # lst_night_m2, lst_night_m1, lst_night_d0, lst_night_p1
  out[[col]] <- NA_real_
}

## ============================================================
## 7. MAIN LOOP OVER UNIQUE FIRE START DATES
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
  ## DAILY LST_NIGHT (WITH LAGS -2, -1, 0, +1)
  ## ---------------------------------------
  for (i_lag in seq_along(lags)) {
    lag      <- lags[i_lag]
    lag_name <- lag_names[i_lag]
    d_target <- d0 + lag
    
    cat("  >> Lag", lag, "-> target date:", format(d_target), "\n")
    
    ## Robust matching using numeric representation of Dates
    layer_idx <- match(as.numeric(d_target), as.numeric(lst_dates))
    if (is.na(layer_idx)) {
      ## No data for that date (remains NA)
      next
    }
    
    r_layer <- r_lst[[layer_idx]]
    
    vals <- extract_with_buffers(
      r_layer,
      fires_buf0, fires_buf5, fires_buf10, fires_buf20,
      idx_fires
    )
    
    colname <- paste0(var_name, "_", lag_name)
    out[[colname]][idx_fires] <- vals
  }
}

## ============================================================
## 8. WRITE OUTPUT
## ============================================================

cat("\nSaving final table to:\n  ", out_file, "\n")
write.csv(out, out_file, row.names = FALSE)
cat("Done.\n")
