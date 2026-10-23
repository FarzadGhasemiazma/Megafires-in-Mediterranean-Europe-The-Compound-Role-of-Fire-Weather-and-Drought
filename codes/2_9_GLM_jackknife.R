# ============================================================
# LOO JACKKNIFE GLM STABILITY ANALYSIS FOR MEGA-FIRES
# Best models only, 5-fold CV × 50 repeats = up to 250 CV fits
# per Mega-fire removal (LOO)
# ============================================================

user_lib <- "~/R/libs"
if (!dir.exists(user_lib)) dir.create(user_lib, recursive = TRUE)
.libPaths(c(user_lib, .libPaths()))

rm(list = ls()); graphics.off(); cat("\f"); options(scipen = 999)

install.packages(c("broom", "car"))

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(caret)
  library(broom)
  library(purrr)
  library(car)
  library(yardstick)
  library(ggplot2)
})

# ------------------------------------------------------------
# SETTINGS
# ------------------------------------------------------------
SEED_MAIN     <- 2025
set.seed(SEED_MAIN)

N_REPEATS     <- 1   # outer repetitions per jackknife
CV_FOLDS      <- 5
CV_REPEATS    <- 50  # repeated CV
THR           <- 0.5
VIF_MAX       <- 2   # relaxed VIF threshold

# Best models for the three transitions:
BEST_MODELS <- list(
  M_to_MG  = c("mean3_lst_night", "mean3_wind_speed"),
  L_to_MG  = c("mean3_lst_night", "mean3_wind_speed"),
  VL_to_MG = c("mean3_lst_night")
)

# ------------------------------------------------------------
# PATHS
# ------------------------------------------------------------
in_csv  <- "/diskonfire/MESOGEOS/weather_3day_means_per_fire_filtered_anom.csv"
out_dir <- "/home/farzad/1_farzad_project/1_Codes/JACKKNIFE_250"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------
# LOAD & PREPARE DATA
# ------------------------------------------------------------
df0 <- readr::read_csv(in_csv, show_col_types = FALSE)

lab_M  <- "Medium [30–100)"
lab_L  <- "Large [100–1000)"
lab_VL <- "Very Large [1000–10000)"
lab_MG <- "Mega (≥10000)"

df <- df0 %>%
  mutate(size_class = case_when(
    area_ha >= 30 & area_ha < 100     ~ lab_M,
    area_ha >= 100 & area_ha < 1000   ~ lab_L,
    area_ha >= 1000 & area_ha < 10000 ~ lab_VL,
    area_ha >= 10000                  ~ lab_MG
  )) %>%
  filter(!is.na(size_class))

TRANSITIONS <- tibble::tribble(
  ~lower, ~upper, ~tag,
  lab_M,  lab_MG, "M_to_MG",
  lab_L,  lab_MG, "L_to_MG",
  lab_VL, lab_MG, "VL_to_MG"
)

# ------------------------------------------------------------
# FUNCTION: FIT GLM WITH CV + DOWNSAMPLING (ONE REPETITION)
# Returns:
#   $metrics : one row per CV run (up to 250 per LOO)
#   $betas   : coefficient estimates per run
# ------------------------------------------------------------
fit_glm_one_rep <- function(df_sub, preds){

  # Pre-VIF check (only if >1 predictor)
  if (length(preds) > 1) {
    gv <- try(
      glm(as.integer(target == "L1") ~ .,
          data   = df_sub[, c(preds, "target")],
          family = binomial),
      silent = TRUE
    )
    if (!inherits(gv, "try-error")) {
      vf <- try(car::vif(gv), silent = TRUE)
      if (!inherits(vf, "try-error") && any(vf > VIF_MAX)) return(NULL)
    }
  }

  folds <- caret::createMultiFolds(df_sub$target, k = CV_FOLDS, times = CV_REPEATS)

  MET <- list()
  BET <- list()
  ii  <- 0

  for (rn in names(folds)) {

    idx_tr <- folds[[rn]]
    idx_te <- setdiff(seq_len(nrow(df_sub)), idx_tr)

    dtr <- df_sub[idx_tr, ]
    dte <- df_sub[idx_te, ]

    if (length(unique(dtr$target)) < 2) next

    # Downsample training to balance L0/L1
    dtr_ds <- caret::downSample(
      x     = dtr[, preds, drop = FALSE],
      y     = dtr$target,
      yname = "target"
    )
    if (length(unique(dtr_ds$target)) < 2) next

    # Fit GLM
    gfit <- glm(as.integer(target == "L1") ~ .,
                data   = dtr_ds[, c(preds, "target")],
                family = binomial)

    # Coefficients (no intercept)
    bet <- broom::tidy(gfit) %>%
      filter(term != "(Intercept)") %>%
      transmute(term, beta_raw = estimate)
    BET[[length(BET) + 1]] <- bet

    # Predict
    prob  <- predict(gfit, newdata = dte[, preds], type = "response")
    truth <- factor(ifelse(dte$target == "L1", "L1", "L0"), levels = c("L1", "L0"))
    pred  <- factor(ifelse(prob >= THR, "L1", "L0"), levels = c("L1", "L0"))

    # Metrics
    AUCv <- yardstick::roc_auc_vec(truth, prob, event_level = "first")
    RECv <- yardstick::sens_vec(truth, pred, event_level = "first")
    FPRv <- 1 - yardstick::spec_vec(truth, pred, event_level = "first")

    ii <- ii + 1
    MET[[ii]] <- tibble(
      cv_run = rn,   # ID of the CV split (fold × repeat)
      AUC    = AUCv,
      Recall = RECv,
      FPR    = FPRv
    )
  }

  if (ii == 0) return(NULL)

  list(
    metrics = bind_rows(MET),
    betas   = bind_rows(BET)
  )
}

# ------------------------------------------------------------
# MAIN JACKKNIFE LOOP
# ------------------------------------------------------------
RESULTS <- list()

for (ti in seq_len(nrow(TRANSITIONS))) {

  lower <- TRANSITIONS$lower[ti]
  upper <- TRANSITIONS$upper[ti]
  tag   <- TRANSITIONS$tag[ti]
  preds <- BEST_MODELS[[tag]]

  message("\n==== Transition: ", tag, " ====")

  df_t <- df %>%
    filter(size_class %in% c(lower, upper)) %>%
    mutate(
      target = factor(
        ifelse(size_class == upper, "L1", "L0"),
        levels = c("L0", "L1")
      )
    ) %>%
    select(all_of(preds), target)

  df_L1 <- df_t %>% filter(target == "L1")
  N_L1  <- nrow(df_L1)
  message("Mega-fires in this transition: ", N_L1)

  # Jackknife: remove one Mega-fire at a time
  for (i in seq_len(N_L1)) {
    message("  Leaving out Mega-fire ", i, "/", N_L1)

    df_sub_jk <- bind_rows(
      df_t %>% filter(target == "L0"),
      df_L1[-i, ]
    )

    # Repeat outer rep if needed (currently 1)
    for (rep in 1:N_REPEATS) {

      out <- fit_glm_one_rep(df_sub_jk, preds)
      key <- paste(tag, "LOO", i, "R", rep, sep = "::")

      RESULTS[[key]] <- list(
        transition = tag,
        loo_id     = i,
        repetition = rep,
        preds      = paste(preds, collapse = "+"),
        out        = out
      )
    }
  }
}

# ------------------------------------------------------------
# SAVE RAW RESULTS OBJECT (OPTIONAL)
# ------------------------------------------------------------
saveRDS(RESULTS, file.path(out_dir, "JACKKNIFE_250_raw.rds"))

# ------------------------------------------------------------
# RAW PERFORMANCE METRICS:
# ONE ROW PER CV RUN (~250 PER LOO)
# ------------------------------------------------------------
METRICS_RAW <- purrr::map_df(RESULTS, function(x){
  if (is.null(x$out)) return(NULL)
  x$out$metrics %>%
    mutate(
      transition = x$transition,
      loo_id     = x$loo_id,
      repetition = x$repetition,
      model      = x$preds
    )
})

# This is the file you care about: 250 rows per LOO (if all folds used)
write_csv(METRICS_RAW, file.path(out_dir, "JACKKNIFE_250_metrics_raw.csv"))

# ------------------------------------------------------------
# (OPTIONAL) MEDIAN METRICS PER LOO
# ------------------------------------------------------------
STAB <- METRICS_RAW %>%
  group_by(transition, loo_id, repetition, model) %>%
  summarise(
    AUC_med    = median(AUC),
    Recall_med = median(Recall),
    FPR_med    = median(FPR),
    .groups    = "drop"
  )

write_csv(STAB, file.path(out_dir, "JACKKNIFE_250_metrics_summary.csv"))

# ------------------------------------------------------------
# β STABILITY (RAW)
# ------------------------------------------------------------
BETA_RAW <- purrr::map_df(RESULTS, function(x){
  if (is.null(x$out)) return(NULL)
  x$out$betas %>%
    mutate(
      transition = x$transition,
      loo_id     = x$loo_id,
      repetition = x$repetition,
      model      = x$preds
    )
})

write_csv(BETA_RAW, file.path(out_dir, "JACKKNIFE_250_beta_raw.csv"))

# ------------------------------------------------------------
# β STABILITY (SUMMARY) - OPTIONAL
# ------------------------------------------------------------
BETA_SUMM <- BETA_RAW %>%
  group_by(transition, model, term) %>%
  summarise(
    beta_med = median(beta_raw),
    beta_L   = quantile(beta_raw, 0.025),
    beta_U   = quantile(beta_raw, 0.975),
    sign_consistency = max(
      mean(beta_raw > 0),
      mean(beta_raw < 0)
    ),
    .groups = "drop"
  )

write_csv(BETA_SUMM, file.path(out_dir, "JACKKNIFE_250_beta_summary.csv"))

message("\nDONE. Results saved in: ", out_dir)
