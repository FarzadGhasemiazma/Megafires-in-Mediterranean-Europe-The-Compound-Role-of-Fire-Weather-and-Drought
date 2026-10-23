# =============================
# CLEAN START & PACKAGES
# =============================
rm(list = ls()); graphics.off(); cat("\f"); options(scipen = 999)

pkgs <- c("sf", "dplyr", "readr", "tidyr", "stringr")
to_install <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
if (length(to_install)) install.packages(to_install, quiet = TRUE)
invisible(lapply(pkgs, library, character.only = TRUE))

# =============================
# PATHS
# =============================
base_dir  <- "/Users/marco/Library/CloudStorage/Dropbox/estcena/scripts/farzad/data"

gpkg_file <- file.path(
  base_dir,
  "all_fires_mesogeos_JJAS_gt30ha_2008_2022_med_countries.gpkg"
)

out_dir   <- file.path(base_dir, "tables_size_classes")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

f_global_counts   <- file.path(out_dir, "fire_size_global_counts_pct_BA.csv")
f_country_counts  <- file.path(out_dir, "fire_size_by_country_counts_pct_BA.csv")

# =============================
# READ GPKG
# =============================
fires_sf <- sf::st_read(gpkg_file, quiet = TRUE)

stopifnot("country" %in% names(fires_sf))
stopifnot("area_ha" %in% names(fires_sf))

fires <- fires_sf |> sf::st_drop_geometry()

# =============================
# NORMALIZE COUNTRY NAMES
# =============================
fires <- fires |>
  dplyr::mutate(
    country_norm = dplyr::case_when(
      country %in% c("Greece", "Ελλάδα")      ~ "Greece",
      country %in% c("Bulgaria", "България")  ~ "Bulgaria",
      country %in% c("Croatia", "Hrvatska")   ~ "Croatia",
      TRUE                                    ~ country
    )
  )

# =============================
# DEFINE SIZE CLASSES
# =============================
class_order <- c("Medium", "Large", "Very Large", "Mega")

fires <- fires |>
  dplyr::mutate(
    size_class = dplyr::case_when(
      area_ha >=   30 & area_ha <    100 ~ "Medium",
      area_ha >=  100 & area_ha <   1000 ~ "Large",
      area_ha >= 1000 & area_ha <  10000 ~ "Very Large",
      area_ha >= 10000                   ~ "Mega",
      TRUE ~ NA_character_
    ),
    size_class = factor(size_class, levels = class_order)
  ) |>
  dplyr::filter(!is.na(size_class))

# =============================
# 1) GLOBAL SUMMARY: frequency + burned area
# =============================
global_summary <- fires |>
  dplyr::group_by(size_class) |>
  dplyr::summarise(
    n_fires      = dplyr::n(),
    BA_total_ha  = sum(area_ha, na.rm = TRUE),
    BA_mean_ha   = mean(area_ha, na.rm = TRUE),
    BA_median_ha = stats::median(area_ha, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::arrange(size_class)

total_fires_global <- sum(global_summary$n_fires)
total_BA_global    <- sum(global_summary$BA_total_ha)

global_summary <- global_summary |>
  dplyr::mutate(
    pct_fires = 100 * n_fires / total_fires_global,
    pct_BA    = 100 * BA_total_ha / total_BA_global
  )

readr::write_csv(global_summary, f_global_counts)
cat("Global (freq + BA) summary saved to:\n  ", f_global_counts, "\n\n")

# =============================
# 2) COUNTRY SUMMARY: frequency + burned area
# =============================
country_class_summary <- fires |>
  dplyr::group_by(country_norm, size_class) |>
  dplyr::summarise(
    n_fires     = dplyr::n(),
    BA_total_ha = sum(area_ha, na.rm = TRUE),
    BA_mean_ha  = mean(area_ha, na.rm = TRUE),
    BA_median_ha = stats::median(area_ha, na.rm = TRUE),
    .groups = "drop"
  )

# add totals per country to compute within-country %
country_class_summary <- country_class_summary |>
  dplyr::group_by(country_norm) |>
  dplyr::mutate(
    total_fires_country = sum(n_fires),
    total_BA_country    = sum(BA_total_ha),
    pct_fires_country   = dplyr::if_else(total_fires_country > 0,
                                         100 * n_fires / total_fires_country, 0),
    pct_BA_country      = dplyr::if_else(total_BA_country > 0,
                                         100 * BA_total_ha / total_BA_country, 0)
  ) |>
  dplyr::ungroup()

# (Optional) wide format, if you want one row per country
country_class_wide <- country_class_summary |>
  dplyr::mutate(
    size_class_clean = stringr::str_replace_all(size_class, " ", "_")
  ) |>
  tidyr::pivot_wider(
    id_cols = c(country_norm,
                total_fires_country, total_BA_country),
    names_from = size_class_clean,
    values_from = c(n_fires, BA_total_ha,
                    BA_mean_ha, BA_median_ha,
                    pct_fires_country, pct_BA_country),
    names_glue = "{size_class_clean}_{.value}",
    values_fill = 0
  ) |>
  dplyr::arrange(country_norm)

readr::write_csv(country_class_wide, f_country_counts)
cat("Country (freq + BA) summary saved to:\n  ", f_country_counts, "\n")

# quick look in console
cat("\n=== GLOBAL SUMMARY (freq + BA) ===\n")
print(global_summary)

cat("\n=== COUNTRY SUMMARY (head; wide) ===\n")
print(head(country_class_wide))
