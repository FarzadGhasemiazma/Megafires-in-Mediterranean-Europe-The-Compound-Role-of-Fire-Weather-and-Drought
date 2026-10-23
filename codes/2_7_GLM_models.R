# ============================================================
# GLM (logit) — 5×50 CV with β CIs from resamples
# - Input predictors already standardized (z-scores / anomalies)
# - 4 predictors: mean3_wind_speed, mean3_lst_night, mean3_ndvi, spei03_interp
# - All 1..4 combinations per transition
# - 5×50 CV:
#     * caret::createMultiFolds (k=5, times=50)
#     * in each resample: downSample(train), GLM on train, predict on test
# - Metrics:
#     * AUC, Recall (thr=0.5), FPR from test per resample
#     * ROC ribbons from resample ROC curves
# - β:
#     * For every resample, save coefficients (no standardization)
#     * Per transition–model–term: median β and 95% empirical CI (2.5–97.5%)
# - VIF filter:
#     * Pre-CV GLM on full data (no downsample), VIF <= 2 if nvars > 1
# - Best models:
#     * By AUC_med (then Recall_med, then nvars)
#     * Filtered best: only models where all β CIs do NOT cross 0
# Outputs:
#   .../OUT_TABLES/GLM_transitions_all_models.csv
#   .../OUT_TABLES/GLM_coefficients_summary.csv   (β CV-based)
#   .../OUT_TABLES/GLM_best_models_by_transition.csv
#   .../OUT_TABLES/GLM_best_models_by_transition_filtered.csv
#   .../OUT_TABLES/VIF_posthoc_log.csv
#   .../ROC_PER_RESAMPLE/ROC_points_<transition>_<model>.csv
#   .../ROC_PER_RESAMPLE/ROC_CIgrid_<transition>_<model>.csv
# ============================================================

rm(list=ls()); graphics.off(); cat("\f"); options(scipen = 999)

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(purrr); library(stringr); library(tibble)
  library(caret); library(pROC); library(yardstick); library(car); library(broom)
})

# ---------- SETTINGS ----------
SEED       <- 1000
N_FOLDS    <- 5
N_REPEATS  <- 50
THRESH     <- 0.5

# Predictores (ya estandarizados)
PRED_VARS  <- c("mean3_wind_speed",   # WS
                "mean3_lst_night",    # LSTn
                "mean3_ndvi",         # NDVI
                "spei03_interp")      # SPEI-3

VIF_MAX    <- 2
FPR_GRID   <- seq(0,1,by=0.01)

# ---------- PATHS ----------
in_base   <- "/Users/marco/Library/CloudStorage/Dropbox/estcena/scripts/farzad"
in_csv    <- file.path(in_base,"data","weather_3day_means_per_fire_filtered_anom.csv")

out_root  <- "/Users/marco/Documents/output/farzad"
out_dir   <- file.path(out_root, "OUT_TABLES");       dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
roc_dir   <- file.path(out_root, "ROC_PER_RESAMPLE"); dir.create(roc_dir, recursive=TRUE, showWarnings=FALSE)

# ---------- DATA & CLASSES ----------
df0 <- readr::read_csv(in_csv, show_col_types = FALSE)

lab_M  <- "Medium [30–100)";         lab_L  <- "Large [100–1000)"
lab_VL <- "Very Large [1000–10000)"; lab_MG <- "Mega (≥10000)"

df <- df0 %>%
  mutate(size_class = case_when(
    area_ha >= 30   & area_ha < 100    ~ lab_M,
    area_ha >= 100  & area_ha < 1000   ~ lab_L,
    area_ha >= 1000 & area_ha < 10000  ~ lab_VL,
    area_ha >= 10000                   ~ lab_MG,
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(size_class))

TRANSITIONS <- tibble::tribble(
  ~lower, ~upper, ~tag,       ~label,
  lab_M,  lab_L,  "M_to_L",   "Medium \u2192 Large",
  lab_M,  lab_VL, "M_to_VL",  "Medium \u2192 Very Large",
  lab_M,  lab_MG, "M_to_MG",  "Medium \u2192 Mega",
  lab_L,  lab_VL, "L_to_VL",  "Large \u2192 Very Large",
  lab_L,  lab_MG, "L_to_MG",  "Large \u2192 Mega",
  lab_VL, lab_MG, "VL_to_MG", "Very Large \u2192 Mega"
)

# ---------- HELPERS ----------
combos_of <- function(vars){
  k <- seq_along(vars)
  unlist(lapply(k, function(m) combn(vars, m, simplify=FALSE)), recursive=FALSE)
}
safe_slug <- function(x){
  x <- gsub("[^A-Za-z0-9_\\-]+","_",x)
  x <- gsub("_+","_",x)
  gsub("^_+|_+$","",x)
}

# ROC points per resample
roc_points_df <- function(truth_char, prob){
  rr <- pROC::roc(response = truth_char,
                  predictor = prob,
                  levels = c("L0","L1"), quiet = TRUE)
  as_tibble(pROC::coords(rr, "all", ret=c("specificity","sensitivity","threshold"))) %>%
    transmute(FPR = 1-specificity, TPR = sensitivity, threshold)
}

# Aggrega a griglia FPR (mediana + CI95)
roc_ci_grid <- function(all_df, fpr_grid=seq(0,1,by=0.01)){
  all_dfs <- split(all_df, all_df$resample_id)
  interp <- lapply(all_dfs, function(d){
    d2 <- d[order(d$FPR),]
    approx(x=d2$FPR, y=d2$TPR, xout=fpr_grid, ties="ordered",
           yleft=0, yright=1)$y
  })
  M <- do.call(cbind, interp)
  tibble(FPR=fpr_grid,
         TPR_med = apply(M,1,median,na.rm=TRUE),
         TPR_L   = apply(M,1,quantile,0.025,na.rm=TRUE),
         TPR_U   = apply(M,1,quantile,0.975,na.rm=TRUE))
}

# Metricas per resample
auc_vec <- function(truth, prob){
  yardstick::roc_auc_vec(truth=truth, estimate=prob, event_level="first")
}
recall_at <- function(truth, prob, thr=0.5){
  cls <- factor(ifelse(prob >= thr, levels(truth)[1], levels(truth)[2]),
                levels=levels(truth))
  yardstick::sens_vec(truth=truth, estimate=cls, event_level="first")
}
fpr_at <- function(truth, prob, thr=0.5){
  cls <- factor(ifelse(prob >= thr, levels(truth)[1], levels(truth)[2]),
                levels=levels(truth))
  1 - yardstick::spec_vec(truth=truth, estimate=cls, event_level="first")
}

# ---------- MAIN ----------
ALL_MODELS   <- list()
ALL_BETA_RES <- list()  # beta per resample
VIF_LOG      <- list()

model_list <- combos_of(PRED_VARS)

for (ti in seq_len(nrow(TRANSITIONS))) {
  lower <- TRANSITIONS$lower[ti]; upper <- TRANSITIONS$upper[ti]
  ttag  <- TRANSITIONS$tag[ti];    tlab  <- TRANSITIONS$label[ti]
  
  message("\n==============================")
  message("Transition: ", tlab, " [", ttag, "]")
  message("==============================")
  
  base_df <- df %>%
    filter(size_class %in% c(lower, upper)) %>%
    select(all_of(PRED_VARS), size_class) %>%
    mutate(
      target = ifelse(size_class == upper, "L1", "L0"),
      target = factor(target, levels = c("L0","L1"))
    ) %>%
    select(-size_class)
  
  if (nrow(base_df) < 20 || length(unique(base_df$target)) < 2) {
    warning("Skipping transition ", ttag, " — insufficient data.")
    next
  }
  
  for (preds in model_list) {
    mod_id <- paste(preds, collapse = "+")
    message("  → Model: ", mod_id)
    
    sub_df <- base_df %>%
      select(all_of(preds), target) %>%
      filter(if_all(all_of(preds), ~ is.finite(.) & !is.na(.)))
    
    if (nrow(sub_df) < 50 || length(unique(sub_df$target)) < 2) {
      message("    Skipped (n<50 or single class).")
      next
    }
    
    # ----- VIF pre-filter (single GLM on full data, no downsample) -----
    vif_stat <- NA_real_
    if (length(preds) > 1L) {
      y_bin <- as.integer(sub_df$target == "L1")
      df_vif <- data.frame(y_bin = y_bin, sub_df[, preds, drop=FALSE])
      glm_vif <- try(glm(y_bin ~ ., data=df_vif, family=binomial), silent=TRUE)
      if (!inherits(glm_vif, "try-error")) {
        vf <- try(car::vif(glm_vif), silent=TRUE)
        if (!inherits(vf, "try-error")) {
          vif_stat <- max(vf, na.rm=TRUE)
        }
      }
    }
    VIF_LOG[[length(VIF_LOG)+1L]] <- tibble(
      transition = ttag,
      model      = mod_id,
      nvars      = length(preds),
      vif_max    = as.numeric(vif_stat)
    )
    if (length(preds) > 1L && is.finite(vif_stat) && vif_stat > VIF_MAX) {
      message("    Skipped (VIF_max > ", VIF_MAX, ").")
      next
    }
    
    # ----- 5×50 CV with downsample + β per resample -----
    set.seed(SEED)
    folds <- caret::createMultiFolds(sub_df$target, k = N_FOLDS, times = N_REPEATS)
    
    MET_R   <- list()
    ROC_all <- list()
    BETA_r  <- list()
    
    r_index <- 0L
    for (rn in names(folds)) {
      idx_tr <- folds[[rn]]
      idx_te <- setdiff(seq_len(nrow(sub_df)), idx_tr)
      
      dtr <- sub_df[idx_tr, , drop=FALSE]
      dte <- sub_df[idx_te, , drop=FALSE]
      
      if (length(unique(dtr$target)) < 2L || length(unique(dte$target)) < 2L) next
      
      # Downsample training set
      dtr_ds <- caret::downSample(
        x = dtr[, preds, drop=FALSE],
        y = dtr$target,
        yname = "target"
      )
      if (length(unique(dtr_ds$target)) < 2L) next
      
      # GLM on standardized predictors (already z-scores)
      y_bin_tr <- as.integer(dtr_ds$target == "L1")
      df_glm   <- data.frame(y_bin = y_bin_tr, dtr_ds[, preds, drop=FALSE])
      
      fit_glm <- try(glm(y_bin ~ ., data=df_glm, family=binomial), silent=TRUE)
      if (inherits(fit_glm, "try-error")) next
      
      # Coeff per resample (no intercept)
      bt <- try(broom::tidy(fit_glm), silent=TRUE)
      if (!inherits(bt, "try-error") && nrow(bt)) {
        bt2 <- bt %>%
          filter(term != "(Intercept)") %>%
          transmute(
            transition = ttag,
            model      = mod_id,
            resample_id= rn,
            term       = term,
            beta_raw   = estimate
          )
        BETA_r[[length(BETA_r)+1L]] <- bt2
      }
      
      # Predictions on test
      prob <- try(predict(fit_glm, newdata=dte[, preds, drop=FALSE], type="response"),
                  silent=TRUE)
      if (inherits(prob, "try-error")) next
      
      truth <- factor(ifelse(dte$target == "L1","L1","L0"), levels=c("L1","L0"))
      if (length(unique(truth)) < 2L) next
      
      AUC   <- auc_vec(truth, prob)
      RECALL<- recall_at(truth, prob, THRESH)
      FPR   <- fpr_at(truth, prob, THRESH)
      
      r_index <- r_index + 1L
      MET_R[[r_index]] <- tibble(
        Resample = rn,
        AUC      = AUC,
        Recall   = RECALL,
        FPR      = FPR
      )
      
      # ROC points for this resample
      rp <- roc_points_df(truth_char = as.character(truth), prob = prob)
      rp$resample_id <- rn
      ROC_all[[length(ROC_all)+1L]] <- rp
    }
    
    if (!length(MET_R)) {
      message("    No valid resamples, skipping.")
      next
    }
    
    MET_R_df <- bind_rows(MET_R)
    
    # Aggregate metrics across resamples
    AGG <- MET_R_df %>%
      summarise(
        AIC_med   = NA_real_, AIC_L = NA_real_, AIC_U = NA_real_,  # placeholder
        AUC_med   = median(AUC,   na.rm=TRUE),
        AUC_L     = quantile(AUC,   0.025, na.rm=TRUE),
        AUC_U     = quantile(AUC,   0.975, na.rm=TRUE),
        RECALL_med= median(Recall, na.rm=TRUE),
        RECALL_L  = quantile(Recall, 0.025, na.rm=TRUE),
        RECALL_U  = quantile(Recall, 0.975, na.rm=TRUE),
        FPR_med   = median(FPR,    na.rm=TRUE)
      ) %>%
      mutate(
        transition = ttag,
        label      = tlab,
        model      = mod_id,
        nvars      = length(preds),
        .before = 1
      )
    
    ALL_MODELS[[paste(ttag,mod_id,sep="::")]] <- AGG
    
    # β per resample (for CI)
    if (length(BETA_r)) {
      ALL_BETA_RES[[paste(ttag,mod_id,sep="::")]] <- bind_rows(BETA_r)
    }
    
    # ROC save
    if (length(ROC_all)) {
      ROC_df <- bind_rows(ROC_all)
      ROC_df$transition <- ttag
      ROC_df$model      <- mod_id
      
      readr::write_csv(ROC_df,
                       file.path(roc_dir, paste0("ROC_points_", ttag, "_", safe_slug(mod_id), ".csv")))
      
      ci_grid <- roc_ci_grid(ROC_df %>% select(resample_id, FPR, TPR))
      ci_grid$transition <- ttag
      ci_grid$model      <- mod_id
      
      readr::write_csv(ci_grid,
                       file.path(roc_dir, paste0("ROC_CIgrid_", ttag, "_", safe_slug(mod_id), ".csv")))
    }
  }
}

# ---------- SAVE TABLES ----------
TAB_ALL <- bind_rows(ALL_MODELS)
readr::write_csv(TAB_ALL, file.path(out_dir, "GLM_transitions_all_models.csv"))

VIF_TAB <- bind_rows(VIF_LOG)
readr::write_csv(VIF_TAB, file.path(out_dir, "VIF_posthoc_log.csv"))

# β distributions over 5×50 resamples
if (length(ALL_BETA_RES)) {
  BETA_RES_ALL <- bind_rows(ALL_BETA_RES)
} else {
  BETA_RES_ALL <- tibble()
}

# Summaries: median + 95% empirical interval for β (no extra standardization)
if (nrow(BETA_RES_ALL)) {
  COEF_SUMM <- BETA_RES_ALL %>%
    group_by(transition, model, term) %>%
    summarise(
      beta_med = median(beta_raw, na.rm=TRUE),
      beta_L   = quantile(beta_raw, 0.025, na.rm=TRUE),
      beta_U   = quantile(beta_raw, 0.975, na.rm=TRUE),
      .groups  = "drop"
    )
} else {
  COEF_SUMM <- tibble()
}

readr::write_csv(COEF_SUMM, file.path(out_dir, "GLM_coefficients_summary.csv"))

# ---------- BEST MODELS (unfiltered) ----------
BEST <- TAB_ALL %>%
  group_by(transition) %>%
  arrange(desc(AUC_med), desc(RECALL_med), nvars, .by_group=TRUE) %>%
  slice(1) %>% ungroup()
readr::write_csv(BEST, file.path(out_dir, "GLM_best_models_by_transition.csv"))

# ---------- FILTERED BEST: all β CIs NOT crossing 0 ----------
if (nrow(COEF_SUMM)) {
  ok_models <- COEF_SUMM %>%
    group_by(transition, model) %>%
    summarise(all_sig = all(!(beta_L < 0 & beta_U > 0)), .groups="drop") %>%
    filter(all_sig)
  
  TAB_F <- TAB_ALL %>% inner_join(ok_models, by=c("transition","model"))
  
  BEST_F <- TAB_F %>%
    group_by(transition) %>%
    arrange(desc(AUC_med), desc(RECALL_med), nvars, .by_group=TRUE) %>%
    slice(1) %>% ungroup()
  
  readr::write_csv(BEST_F, file.path(out_dir, "GLM_best_models_by_transition_filtered.csv"))
} else {
  BEST_F <- tibble()
  readr::write_csv(BEST_F, file.path(out_dir, "GLM_best_models_by_transition_filtered.csv"))
}

cat("Saved in:\n - ", file.path(out_dir, "GLM_transitions_all_models.csv"),
    "\n - ", file.path(out_dir, "GLM_coefficients_summary.csv"),
    "\n - ", file.path(out_dir, "GLM_best_models_by_transition.csv"),
    "\n - ", file.path(out_dir, "GLM_best_models_by_transition_filtered.csv"),
    "\n - ", file.path(out_dir, "VIF_posthoc_log.csv"),
    "\nROC data in:\n - ", roc_dir, "\n\nDone.\n")
