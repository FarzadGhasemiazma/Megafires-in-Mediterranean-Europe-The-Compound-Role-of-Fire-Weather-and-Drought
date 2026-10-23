# =============================
# CLEAN START & PACKAGES
# =============================
rm(list = ls()); graphics.off(); cat("\f"); options(scipen = 999)

pkgs <- c("readr","dplyr","tidyr","stringr","purrr","tibble")
to_install <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
if (length(to_install)) install.packages(to_install, quiet = TRUE)
invisible(lapply(pkgs, library, character.only = TRUE))

# =============================
# INPUT
# =============================
csv_path <- "/Users/marco/Library/CloudStorage/Dropbox/estcena/scripts/farzad/data/fire_meteo_JJAS_zscore_SPEI_2008_2022.csv"

# Base variable names whose z-score windows we want to average over (d-1, d0, d+1)
window_var_bases <- c("rh","t2m","tp","wind_speed","smi",
                      "lst_day","lst_night","lai","ndvi")

# SPEI variables (no averaging over days)
spei_vars <- c("spei01_interp","spei03_interp","spei06_interp")

# Day suffixes for the z-score windows (we will add the "_z" inside the function)
day_suffixes <- c("_m1","_d0","_p1")

# Fire size class labels
lab_M  <- "Medium [30–100)"
lab_L  <- "Large [100–1000)"
lab_VL <- "Very Large [1000–10000)"
lab_MG <- "Mega (≥10000)"

# =============================
# READ + BASIC PREP
# =============================
dat <- readr::read_csv(csv_path, show_col_types = FALSE)

# Minimal checks for required columns
stopifnot(all(c("fire_id","date_start","area_ha","x","y") %in% names(dat)))

# Add date components and size class
dat <- dat %>%
  mutate(
    date_start = as.Date(date_start),
    year  = as.integer(format(date_start, "%Y")),
    month = as.integer(format(date_start, "%m")),
    size_class = dplyr::case_when(
      area_ha >= 30    & area_ha < 100    ~ lab_M,
      area_ha >= 100   & area_ha < 1000   ~ lab_L,
      area_ha >= 1000  & area_ha < 10000  ~ lab_VL,
      area_ha >= 10000                    ~ lab_MG,
      TRUE ~ NA_character_
    )
  )

# =============================
# FILTERS (for safety, even if JJAS is already in the file name)
# - years 2008–2022
# - months June–September (6:9)
# - area >= 30 ha
# =============================
dat_filt <- dat %>%
  filter(year >= 2008, year <= 2022,
         month %in% 6:9,
         area_ha >= 30)

# =============================
# HELPER: 3-day mean from z-score windows
# For each base var (e.g. "rh"), we look for:
#   rh_z_m1, rh_z_d0, rh_z_p1
# and compute their row-wise mean (na.rm = TRUE).
# New column name will correspond to "mean3_<base_var>"
# (e.g. mean3_rh, mean3_wind_speed, etc.).
# =============================
decision_log <- list()

compute_mean3_from_zwindow <- function(df, base_var,
                                       day_suffs = day_suffixes) {
  # Name of the 3-day-mean variable we want to create
  target_name <- base_var
  
  # Actual column names in the data, e.g. "rh_z_m1"
  cols_win <- paste0(base_var, "_z", day_suffs)
  
  existing <- cols_win[cols_win %in% names(df)]
  
  if (length(existing) == 0) {
    # No columns found: return all NA
    val <- rep(NA_real_, nrow(df))
    decision_log[[base_var]] <<- "missing_all_z_window_columns"
  } else {
    m <- as.matrix(df[existing])
    val <- rowMeans(m, na.rm = TRUE)
    # If all three were NA, rowMeans gives NaN -> set back to NA
    val[is.nan(val)] <- NA_real_
    decision_log[[base_var]] <<- paste0("computed_from(", paste(existing, collapse = ","), ")")
  }
  
  tibble(.row_id = seq_len(nrow(df)),
         var     = target_name,
         value   = as.numeric(val))
}

# =============================
# COMPUTE:
# - 3-day mean for all z-score window variables
# - direct values for SPEI (no averaging)
# =============================
means_long <- purrr::map_dfr(window_var_bases,
                             ~compute_mean3_from_zwindow(dat_filt, .x))

# SPEI (no averaging)
spei_present <- intersect(spei_vars, names(dat_filt))

if (length(spei_present) == 0) {
  warning("No SPEI variables found in the dataset.")
  spei_long <- tibble(.row_id = integer(), var = character(), value = numeric())
} else {
  spei_long <- dat_filt %>%
    mutate(.row_id = dplyr::row_number()) %>%
    select(.row_id, all_of(spei_present)) %>%
    tidyr::pivot_longer(
      cols      = all_of(spei_present),
      names_to  = "var",
      values_to = "value"
    )
}

# Combined long for NA diagnostics (variables + SPEI)
combined_long <- dplyr::bind_rows(means_long, spei_long)

# =============================
# BUILD WIDE TABLE:
# - mean3_<var> for window variables
# - spei01/03/06_interp (no prefix)
# - plus fire metadata
# =============================

# 3-day means in wide format. Column names: mean3_<base_var>
wide_means <- means_long %>%
  select(.row_id, var, value) %>%
  tidyr::pivot_wider(
    names_from   = var,
    values_from  = value,
    names_prefix = "mean3_"
  )

# SPEI wide format, if present
wide_spei <- dat_filt %>%
  mutate(.row_id = dplyr::row_number()) %>%
  {
    if (length(spei_present) > 0) {
      select(., .row_id, all_of(spei_present))
    } else {
      select(., .row_id)
    }
  }

# Build meta columns (optionally including 'country' if present)
meta_cols_base  <- c(".row_id","fire_id","date_start","year","month")
if ("country" %in% names(dat_filt)) {
  meta_cols_base <- c(meta_cols_base, "country")
}
meta_cols_rest <- c("area_ha","size_class","x","y")
meta_cols      <- c(meta_cols_base, meta_cols_rest)

meta_fire <- dat_filt %>%
  mutate(.row_id = dplyr::row_number()) %>%
  select(all_of(meta_cols))

# Final wide table: means + SPEI + metadata
means_wide_filtered <- wide_means %>%
  left_join(wide_spei, by = ".row_id") %>%
  left_join(meta_fire,  by = ".row_id")

# =============================
# NA DIAGNOSTICS (AFTER 3-DAY MEAN)
# =============================

# NA by variable
na_after_mean_by_var <- combined_long %>%
  group_by(var) %>%
  summarise(
    fires_total    = n(),
    fires_value_na = sum(is.na(value)),
    pct_value_na   = round(100 * fires_value_na / fires_total, 2),
    fires_value_ok = fires_total - fires_value_na,
    .groups        = "drop"
  ) %>%
  arrange(var)

# NA by fire
na_after_mean_by_fire <- combined_long %>%
  group_by(.row_id) %>%
  summarise(
    n_vars            = n(),
    n_vars_value_na   = sum(is.na(value)),
    pct_vars_value_na = round(100 * n_vars_value_na / n_vars, 2),
    .groups           = "drop"
  ) %>%
  left_join(meta_fire, by = ".row_id")

# Add NA summary per fire into the main wide table
means_wide_filtered <- means_wide_filtered %>%
  left_join(
    na_after_mean_by_fire %>%
      select(fire_id, n_vars_value_na, pct_vars_value_na),
    by = "fire_id"
  )

# =============================
# OPTIONAL: COUNTRY COVERAGE/BA SUMMARY
# (only if 'country' column exists)
# =============================
country_coverage <- NULL
country_ba_summary <- NULL

if ("country" %in% names(dat_filt)) {
  # Table 1 — coverage per country
  country_coverage <- dat_filt %>%
    group_by(country) %>%
    summarise(
      first_year_with_data = min(year, na.rm = TRUE),
      last_year_with_data  = max(year, na.rm = TRUE),
      n_fires              = dplyr::n(),
      total_ba_all_years   = sum(area_ha, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(country)
  
  # Table 2 — burned area summary per country
  ba_yearly <- dat_filt %>%
    group_by(country, year) %>%
    summarise(
      ba_year = sum(area_ha, na.rm = TRUE),
      .groups = "drop"
    )
  
  country_ba_summary <- ba_yearly %>%
    group_by(country) %>%
    summarise(
      total_ba_country       = sum(ba_year, na.rm = TRUE),
      first_year_ba_positive = suppressWarnings(min(year[ba_year > 0], na.rm = TRUE)),
      last_year_ba_positive  = suppressWarnings(max(year[ba_year > 0], na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    mutate(
      first_year_ba_positive = ifelse(is.infinite(first_year_ba_positive), NA, first_year_ba_positive),
      last_year_ba_positive  = ifelse(is.infinite(last_year_ba_positive),  NA, last_year_ba_positive)
    ) %>%
    arrange(country)
} else {
  message('Column "country" not found in data: skipping country coverage summaries.')
}

# =============================
# OUTPUTS (keep main names as before)
# =============================
out_dir <- dirname(csv_path)

files_written <- c()

# Main per-fire table
f1 <- file.path(out_dir, "weather_3day_means_per_fire_filtered_anom.csv")
readr::write_csv(means_wide_filtered, f1)
files_written <- c(files_written, f1)

# NA summaries
f2 <- file.path(out_dir, "weather_na_AFTERmean_by_variable_anom.csv")
readr::write_csv(na_after_mean_by_var, f2)
files_written <- c(files_written, f2)

f3 <- file.path(out_dir, "weather_na_AFTERmean_by_fire_anom.csv")
readr::write_csv(na_after_mean_by_fire, f3)
files_written <- c(files_written, f3)

# Optional country tables
if (!is.null(country_coverage)) {
  f4 <- file.path(out_dir, "country_coverage_first_last_year_anom.csv")
  readr::write_csv(country_coverage, f4)
  files_written <- c(files_written, f4)
}
if (!is.null(country_ba_summary)) {
  f5 <- file.path(out_dir, "country_ba_summary_total_firstlast_nonzero_anom.csv")
  readr::write_csv(country_ba_summary, f5)
  files_written <- c(files_written, f5)
}

cat("\nFiles written:\n - ",
    paste(files_written, collapse = "\n - "),
    "\n", sep = "")

# =============================
# DIAGNOSTIC: which path used for each base variable?
# =============================
diag_tbl <- tibble::tibble(
  var_base  = names(decision_log),
  path_used = unlist(decision_log, use.names = FALSE)
) %>%
  arrange(var_base)

print(diag_tbl, n = Inf)
