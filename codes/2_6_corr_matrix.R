# ============================
# Correlation matrix (pooled) — ANOM (z-scores + SPEI_interp)
# ============================
rm(list = setdiff(ls(), "df"))  # keep df if already in memory

req <- c("PerformanceAnalytics","dplyr","readr")
to_install <- req[!req %in% installed.packages()[,"Package"]]
if (length(to_install)) install.packages(to_install, quiet = TRUE)
invisible(lapply(req, library, character.only = TRUE))

# --- INPUT ---
base_dir  <- "/Users/marco/Library/CloudStorage/Dropbox/estcena/scripts/farzad"
csv_means <- file.path(base_dir, "data", "weather_3day_means_per_fire_filtered_anom.csv")

# If df doesn't exist OR is not a data frame, (re)load it
if (!exists("df") || !is.data.frame(df)) {
  df <- readr::read_csv(csv_means, show_col_types = FALSE)
} else {
  # Convert tibbles/data.tables etc. to plain data.frame if needed
  df <- as.data.frame(df)
}

# --- helper: base assignment only (no mutate) ---
add_if_exists <- function(.data, new, src, fun) {
  if (!is.data.frame(.data)) stop("Input is not a data.frame.")
  if (src %in% names(.data)) {
    .data[[new]] <- fun(.data[[src]])
  } else {
    .data[[new]] <- NA_real_
  }
  .data
}

# Aquí los mean3_* ya son z-scores, así que usamos identidad
corr_df <- df
corr_df <- add_if_exists(corr_df, "AT_z",    "mean3_t2m",        identity)
corr_df <- add_if_exists(corr_df, "LSTd_z",  "mean3_lst_day",    identity)
corr_df <- add_if_exists(corr_df, "LSTn_z",  "mean3_lst_night",  identity)
corr_df <- add_if_exists(corr_df, "RH_z",    "mean3_rh",         identity)
corr_df <- add_if_exists(corr_df, "WS_z",    "mean3_wind_speed", identity)
corr_df <- add_if_exists(corr_df, "PR_z",    "mean3_tp",         identity)
corr_df <- add_if_exists(corr_df, "SMI_z",   "mean3_smi",        identity)
corr_df <- add_if_exists(corr_df, "NDVI_z",  "mean3_ndvi",       identity)
corr_df <- add_if_exists(corr_df, "LAI_z",   "mean3_lai",        identity)

# SPEI interpolados
corr_df <- add_if_exists(corr_df, "SPEI_1",  "spei01_interp",    identity)
corr_df <- add_if_exists(corr_df, "SPEI_3",  "spei03_interp",    identity)
corr_df <- add_if_exists(corr_df, "SPEI_6",  "spei06_interp",    identity)

# Select only the desired numeric columns that actually exist
keep_vars <- c("AT_z","LSTd_z","LSTn_z","RH_z","WS_z","PR_z",
               "SMI_z","NDVI_z","LAI_z","SPEI_1","SPEI_3","SPEI_6")
corr_df <- corr_df[, intersect(keep_vars, names(corr_df)), drop = FALSE]
corr_df <- corr_df[, vapply(corr_df, is.numeric, logical(1)), drop = FALSE]

# Drop all-NA columns
all_na_cols <- names(corr_df)[colSums(!is.na(corr_df)) == 0]
if (length(all_na_cols)) {
  message("Dropping all-NA columns: ", paste(all_na_cols, collapse = ", "))
  corr_df <- corr_df[, setdiff(names(corr_df), all_na_cols), drop = FALSE]
}

# Keep only rows with finite values across remaining columns
stopifnot(ncol(corr_df) > 0)
finite_mask <- Reduce("&", lapply(corr_df, function(x) is.finite(x)))
corr_df <- corr_df[finite_mask, , drop = FALSE]

print(colnames(corr_df))
cat("N used for correlations:", nrow(corr_df), "\n")

# --- Plot & Save ---
out_pdf <- file.path(base_dir, "FIGS", "Correlation_matrix_pooled_anom.pdf")
dir.create(dirname(out_pdf), showWarnings = FALSE, recursive = TRUE)

PerformanceAnalytics::chart.Correlation(corr_df, histogram = TRUE, pch = 19)

pdf(out_pdf, width = 8.2, height = 8.2)
PerformanceAnalytics::chart.Correlation(corr_df, histogram = TRUE, pch = 19)
dev.off()
cat("Saved: ", out_pdf, "\n", sep = "")

# Also export numeric matrix for reproducibility
cor_mat <- stats::cor(corr_df, use = "pairwise.complete.obs", method = "pearson")
readr::write_csv(as.data.frame(cor_mat),
                 file.path(base_dir, "FIGS", "Correlation_matrix_pooled_anom_values.csv"))
