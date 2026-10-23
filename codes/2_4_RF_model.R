# ============================
# PFI-only pipeline (fast vs slow colors, NEW anom data)
# ============================
rm(list=ls()); graphics.off(); cat("\f"); options(scipen = 999)

pkgs <- c("readr","dplyr","tibble","stringr","tidyr","ggplot2","patchwork",
          "caret","ranger","vip","yardstick","tidytext")
to_install <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
if(length(to_install)) install.packages(to_install, quiet=TRUE)
invisible(lapply(pkgs, library, character.only=TRUE))

# ---- SETTINGS
EXCLUDE_SPEI6  <- FALSE  # set TRUE if you want to drop spei06_interp entirely
B_PFI          <- 250    # permutations
NTREES         <- 500
TR_FRAC        <- 0.80
SEED           <- 2025
DO_DOWNSAMPLE  <- TRUE

# ---- PATHS
base_dir  <- "/Users/marco/Library/CloudStorage/Dropbox/estcena/scripts/farzad"
csv_means <- file.path(base_dir, "data", "weather_3day_means_per_fire_filtered_anom.csv")
out_dir   <- file.path(base_dir, "FIGS"); dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

# ---- DATA
df <- readr::read_csv(csv_means, show_col_types = FALSE)

# ---- SIZE CLASSES
lab_M  <- "Medium [30–100)";     lab_L  <- "Large [100–1000)"
lab_VL <- "Very Large [1000–10000)"; lab_MG <- "Mega (≥10000)"

df <- df |>
  dplyr::mutate(size_class = dplyr::case_when(
    area_ha >= 30   & area_ha < 100    ~ lab_M,
    area_ha >= 100  & area_ha < 1000   ~ lab_L,
    area_ha >= 1000 & area_ha < 10000  ~ lab_VL,
    area_ha >= 10000                   ~ lab_MG,
    TRUE ~ NA_character_
  )) |>
  dplyr::filter(!is.na(size_class))

# ---- PREDICTORS: mean3_* + spei01/03/06_interp (optional drop spei06_interp)
mean3_all <- grep("^mean3_", names(df), value = TRUE)
spei_all  <- intersect(c("spei01_interp","spei03_interp","spei06_interp"), names(df))
if (EXCLUDE_SPEI6) spei_all <- setdiff(spei_all, "spei06_interp")

all_pred <- c(mean3_all, spei_all)

# drop VPD/D2M/DEW just in case
drop_patterns <- "(^|_)vpd($|_)|(^|_)d2m($|_)|dew"
pred_vars <- all_pred[!grepl(drop_patterns, all_pred, ignore.case = TRUE)]

# ---- NA AUDIT + FILTER (keep only rows with finite values across chosen predictors)
na_audit <- tibble::tibble(
  variable = pred_vars,
  n_na = vapply(df[pred_vars], function(x) sum(!is.finite(x) | is.na(x)), integer(1)),
  n_total = nrow(df)
) |>
  dplyr::mutate(na_pct = 100*n_na/pmax(1,n_total)) |>
  dplyr::arrange(dplyr::desc(na_pct))
readr::write_csv(na_audit, file.path(out_dir, "NA_audit_predictors.csv"))

df <- df |>
  dplyr::filter(dplyr::if_all(dplyr::all_of(pred_vars), ~ is.finite(.) & !is.na(.)))

# ---- LABELING & GROUPING
pretty_var <- function(x){
  y <- x
  # quitar prefijos/sufijos técnicos
  y <- gsub("^mean3_", "", y)
  y <- gsub("_interp$", "", y, ignore.case = TRUE)
  
  # renombrar variables
  y <- gsub("(^|_)tp($|_)",       "Precipitation",             y, ignore.case = TRUE)
  y <- gsub("lst[_ ]?day",        "LST Day Temperature",       y, ignore.case = TRUE)
  y <- gsub("lst[_ ]?night",      "LST Night Temperature",     y, ignore.case = TRUE)
  y <- gsub("\\b(tas|t2m)\\b",    "Air Temperature",           y, ignore.case = TRUE)
  y <- gsub("\\brh\\b",           "Relative Humidity",         y, ignore.case = TRUE)
  y <- gsub("\\bsmi\\b",          "Soil Moisture (SMI)",       y, ignore.case = TRUE)
  y <- gsub("\\blai\\b",          "Leaf Area Index (LAI)",     y, ignore.case = TRUE)
  y <- gsub("\\bndvi\\b",         "NDVI",                      y, ignore.case = TRUE)
  y <- gsub("\\bspei01\\b",       "SPEI-1",                    y, ignore.case = TRUE)
  y <- gsub("\\bspei03\\b",       "SPEI-3",                    y, ignore.case = TRUE)
  y <- gsub("\\bspei06\\b",       "SPEI-6",                    y, ignore.case = TRUE)
  
  gsub("_"," ", y)
}

# Two-color speed class:
# Slow = mean3_smi, mean3_ndvi, mean3_lai, spei*
# Fast = all other predictors
speed_group_of <- function(v){
  vlow <- tolower(v)
  is_slow <- grepl("^spei", vlow) |
    grepl("(^|_)smi($|_)|(^|_)ndvi($|_)|(^|_)lai($|_)", vlow)
  ifelse(is_slow, "Slow", "Fast")
}

# colors for speed group
pal_speed <- c(Fast = "#8B0000", Slow = "#555555")  # dark red, dark grey

# ---- METRICS & HELPERS
auc_metric <- function(truth, estimate){
  yardstick::roc_auc_vec(truth = truth, estimate = estimate, event_level = "first")
}
pred_fun <- function(object, newdata){
  pr <- predict(object, data = newdata)$predictions
  pr[,"L1"]
}
safe_slug <- function(x){
  x <- gsub("[^A-Za-z0-9_\\-]+", "_", x)
  x <- gsub("_+", "_", x)
  gsub("^_+|_+$","",x)
}
save_csv_safe <- function(df_, prefix, tag){
  readr::write_csv(df_, file.path(out_dir, paste0(prefix, "_", safe_slug(tag), ".csv")))
}
save_pdf_safe <- function(plot, prefix, tag, w=26, h=22, units="cm"){
  ggplot2::ggsave(file.path(out_dir, paste0(prefix, "_", safe_slug(tag), ".pdf")),
                  plot, width=w, height=h, units=units)
}

# ---- PROGRESS BAR
PROG <- new.env(parent=emptyenv())
init_progress <- function(n_transitions, B_pfi){
  PROG$total_ops <- as.integer(n_transitions * B_pfi)
  PROG$done_ops  <- 0L
  PROG$pb <- utils::txtProgressBar(min=0, max=PROG$total_ops, style=3)
}
tick_progress <- function(n=1L, msg=NULL){
  PROG$done_ops <- PROG$done_ops + as.integer(n)
  if (!is.null(msg)) message(sprintf("[%,d/%,d] %s", PROG$done_ops, PROG$total_ops, msg))
  utils::setTxtProgressBar(PROG$pb, PROG$done_ops); flush.console()
}
finish_progress <- function(){
  close(PROG$pb)
  message(sprintf("✅ Finished %d / %d.", PROG$done_ops, PROG$total_ops))
}

# ---- CORE: PFI (B splits)
perm_importance_with_ci <- function(data, lower_label, upper_label, tag, pretty_title,
                                    B=250, train_frac=0.8, ntree=400, seed=2025,
                                    DO_DOWNSAMPLE=FALSE,
                                    K_FOLDS=5, N_REPEATS=50, MTRY=1){
  
  sub <- data |>
    dplyr::filter(size_class %in% c(lower_label, upper_label)) |>
    dplyr::select(dplyr::all_of(pred_vars), size_class) |>
    dplyr::mutate(target = factor(ifelse(size_class==upper_label,"L1","L0")),
                  .keep="unused")
  sub$target <- stats::relevel(sub$target, ref="L1")
  
  if (length(unique(sub$target)) < 2L || ncol(sub) < 3L) {
    warning("Transition ", tag, ": insufficient data.")
    empty_auc <- tibble::tibble(transition=tag,label=pretty_title,method="pfi_cv",AUC=NA_real_)
    empty_ci  <- tibble::tibble(method="pfi_cv",mean=NA,sd=NA,p2.5=NA,p97.5=NA,
                                transition=tag,label=pretty_title)
    save_csv_safe(empty_auc,"AUC_dist",paste0(tag,"_pfi"))
    save_csv_safe(empty_ci,"AUC_CI",paste0(tag,"_pfi"))
    save_csv_safe(tibble::tibble(), "PFI_CI", tag)
    return(list(auc_df=empty_auc, auc_summary=empty_ci, pfi=tibble::tibble()))
  }
  
  # ---- 5×50 repeated CV: caret returns TRAIN indices
  set.seed(seed)
  folds <- caret::createMultiFolds(sub$target, k = K_FOLDS, times = N_REPEATS)
  B_eff <- length(folds)  # = 250
  
  all_vi   <- vector("list", B_eff)
  auc_base <- rep(NA_real_, B_eff)
  
  for (b in seq_len(B_eff)) {
    
    idx_tr <- folds[[b]]
    dtr <- sub[idx_tr, , drop = FALSE]
    dte <- sub[-idx_tr, , drop = FALSE]
    if (length(unique(dte$target)) < 2L) { tick_progress(1); next }
    
    feats <- setdiff(colnames(dtr), "target")
    if (!length(feats)) { tick_progress(1); next }
    
    # Zero-variance guard (pre)
    zv <- vapply(dtr[, feats, drop = FALSE], function(x) length(unique(x)) <= 1L, logical(1))
    if (any(zv)) {
      feats <- feats[!zv]
      if (!length(feats)) { tick_progress(1); next }
      dtr <- dtr[, c("target", feats), drop = FALSE]
      dte <- dte[, c("target", feats), drop = FALSE]
    }
    
    # Down-sampling SOLO nel training fold
    if (DO_DOWNSAMPLE) {
      dtr <- caret::downSample(x = dtr[, feats, drop = FALSE], y = dtr$target, yname = "target")
      dtr$target <- stats::relevel(dtr$target, ref="L1")
      
      # Zero-variance guard (post)
      zv2 <- vapply(dtr[, feats, drop = FALSE], function(x) length(unique(x)) <= 1L, logical(1))
      if (any(zv2)) {
        feats <- feats[!zv2]
        if (!length(feats)) { tick_progress(1); next }
        dtr <- dtr[, c("target", feats), drop = FALSE]
        dte <- dte[, c("target", feats), drop = FALSE]
      }
    }
    
    # Class weights (opzionale; con downsampling diventano ~1)
    tab <- table(dtr$target)
    w <- as.numeric(sum(tab)/(2*tab)); names(w) <- names(tab)
    
    fit <- ranger::ranger(
      target ~ ., data = dtr,
      mtry = MTRY, num.trees = ntree,
      probability = TRUE, classification = TRUE,
      case.weights = w[as.character(dtr$target)],
      respect.unordered.factors = TRUE,
      seed = seed + b
    )
    
    # Baseline AUC on held-out fold
    pr <- predict(fit, data = dte)$predictions[, "L1"]
    auc_base[b] <- auc_metric(dte$target, pr)
    
    # PFI on held-out fold (AUC drop)
    vi_b <- try(vip::vi_permute(
      object = fit,
      train  = dte[, feats, drop = FALSE],
      target = dte$target,
      metric = auc_metric,
      pred_wrapper = function(object, newdata)
        predict(object, data = newdata)$predictions[, "L1"],
      nsim = 1, parallel = FALSE, smaller_is_better = FALSE
    ), silent = TRUE)
    
    if (!inherits(vi_b, "try-error") && !is.null(vi_b) && nrow(vi_b)) {
      vi_b$run <- b
      all_vi[[b]] <- vi_b
    }
    
    tick_progress(1)
  }
  
  # ---- AUC dist & CI
  auc_df <- tibble::tibble(transition=tag, label=pretty_title, method="pfi_cv", AUC=auc_base)
  save_csv_safe(auc_df,"AUC_dist",paste0(tag,"_pfi"))
  auc_ci <- auc_df |>
    dplyr::summarise(
      mean=mean(AUC,na.rm=TRUE),
      sd=sd(AUC,na.rm=TRUE),
      p2.5=quantile(AUC,0.025,na.rm=TRUE),
      p97.5=quantile(AUC,0.975,na.rm=TRUE)
    ) |>
    dplyr::mutate(transition=tag,label=pretty_title,method="pfi_cv",.before=1)
  save_csv_safe(auc_ci,"AUC_CI",paste0(tag,"_pfi"))
  
  # ---- PFI summary
  vi_list <- Filter(Negate(is.null), all_vi)
  pfi_summ <- if (!length(vi_list)) tibble::tibble() else
    dplyr::bind_rows(vi_list) |>
    dplyr::rename(importance_drop=Importance, variable=Variable) |>
    dplyr::mutate(variable=as.character(variable)) |>
    dplyr::group_by(variable) |>
    dplyr::summarise(
      mean_imp = mean(importance_drop, na.rm=TRUE),
      p2.5     = quantile(importance_drop, 0.025, na.rm=TRUE),
      p97.5    = quantile(importance_drop, 0.975, na.rm=TRUE),
      .groups  = "drop"
    ) |>
    dplyr::mutate(
      speed_group     = speed_group_of(variable),
      variable_pretty = pretty_var(variable),
      transition      = tag,
      label           = pretty_title
    ) |>
    dplyr::arrange(dplyr::desc(mean_imp))
  
  save_csv_safe(pfi_summ, "PFI_CI", tag)
  
  return(list(auc_df=auc_df, auc_summary=auc_ci, pfi=pfi_summ))
}

# ---- RUN ALL TRANSITIONS
TRANSITIONS <- tibble::tribble(
  ~lower, ~upper, ~tag,                         ~label,
  lab_M,  lab_L,  "Medium_to_Large",            "Medium → Large",
  lab_M,  lab_VL, "Medium_to_VeryLarge",        "Medium → Very Large",
  lab_M,  lab_MG, "Medium_to_Mega",             "Medium → Mega",
  lab_L,  lab_VL, "Large_to_VeryLarge",         "Large → Very Large",
  lab_L,  lab_MG, "Large_to_Mega",              "Large → Mega",
  lab_VL, lab_MG, "VeryLarge_to_Mega",          "Very Large → Mega"
)

all_auc_dists <- list(); all_auc_CIs <- list(); all_pfi <- list()
init_progress(nrow(TRANSITIONS), B_PFI)

for(i in seq_len(nrow(TRANSITIONS))){
  lo <- TRANSITIONS$lower[i]
  up <- TRANSITIONS$upper[i]
  tg <- TRANSITIONS$tag[i]
  lb <- TRANSITIONS$label[i]
  message("\n==============================\n→ Transition: ", lb, "  [", tg, "]\n==============================")
  pf <- perm_importance_with_ci(df, lo, up, paste0(tg, "_PFI"), lb,
                                ntree=NTREES, seed=SEED,
                                DO_DOWNSAMPLE=TRUE,
                                K_FOLDS=5, N_REPEATS=50, MTRY=1)
  
  all_auc_dists[[tg]] <- pf$auc_df
  all_auc_CIs[[tg]]   <- pf$auc_summary
  all_pfi[[tg]]       <- pf$pfi
}
finish_progress()

AUC_DISTS <- dplyr::bind_rows(all_auc_dists)
AUC_CIs   <- dplyr::bind_rows(all_auc_CIs)
PFI_ALL   <- dplyr::bind_rows(all_pfi)

readr::write_csv(AUC_DISTS, file.path(out_dir, "AUC_dist_ALL_ALL_PFI.csv"))
readr::write_csv(AUC_CIs,   file.path(out_dir, "AUC_CI_ALL_ALL_PFI.csv"))
readr::write_csv(PFI_ALL,   file.path(out_dir, "PFI_CI_ALL_ALL_PFI.csv"))

# ---- PLOT rankings (PFI ± CI) with 2 colors (Fast=dark red; Slow=dark grey)
plot_pfi <- function(summ, top_k=12){
  if (is.null(summ) || !nrow(summ))
    return(ggplot2::ggplot() + ggplot2::theme_void())
  
  sub <- summ |> dplyr::slice_head(n = top_k)
  
  ggplot2::ggplot(
    sub,
    ggplot2::aes(x = mean_imp,
                 y = tidytext::reorder_within(variable_pretty, mean_imp, speed_group),
                 fill = speed_group)
  ) +
    ggplot2::geom_col(width = 0.8) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = p2.5, xmax = p97.5), height = 0) +
    tidytext::scale_y_reordered() +
    ggplot2::scale_fill_manual(values = pal_speed) +
    ggplot2::labs(
      title = unique(sub$label),
      x = "Mean AUC decrease on permutation (higher = more important)",
      y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.title    = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor   = ggplot2::element_blank(),
      plot.margin    = ggplot2::margin(6, 8, 6, 6),
      plot.title     = ggplot2::element_text(size=11)
    )
}

blank_panel <- function(){ ggplot2::ggplot() + ggplot2::theme_void() }

g11 <- plot_pfi(all_pfi[["Medium_to_Large"]])
g12 <- plot_pfi(all_pfi[["Medium_to_VeryLarge"]])
g13 <- plot_pfi(all_pfi[["Medium_to_Mega"]])

g21 <- blank_panel()
g22 <- plot_pfi(all_pfi[["Large_to_VeryLarge"]])
g23 <- plot_pfi(all_pfi[["Large_to_Mega"]])

g31 <- blank_panel()
g32 <- blank_panel()
g33 <- plot_pfi(all_pfi[["VeryLarge_to_Mega"]])

g_3x3 <- (g11 + g12 + g13) /
  (g21 + g22 + g23) /
  (g31 + g32 + g33) +
  patchwork::plot_layout(guides = "collect") &
  ggplot2::theme(legend.position = "bottom")

ggplot2::ggsave(file.path(out_dir, "PFI_CI_PFI_rankings_3x3_onlyPFI_anomRF.pdf"),
                g_3x3, width=26, height=22, units="cm")

message("✅ Done. CSVs and PDFs in: ", out_dir)
