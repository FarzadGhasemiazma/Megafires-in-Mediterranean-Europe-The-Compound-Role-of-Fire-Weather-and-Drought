# ============================================================
# PFI (Random Forest) — 2×3 rankings with CI (2 colors: Fast vs Slow)
#  - Per-panel ranking preserved (each panel has its own y-order)
#  - Short y labels (acronyms: AT, LSTd, LSTn, RH, WS, PR, SMI, NDVI, LAI, SPEI-1/3/6)
#  - Font size 12 everywhere
#  - 2 rows × 3 columns
#  - Legend collected at bottom (single legend)
#  - Saves: FIGS/PFI_rankings_2x3_RF_SPEED2COLORS_anom.[pdf|png]
# ============================================================

rm(list=ls()); graphics.off(); cat("\f"); options(scipen=999)

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(stringr); library(tidyr)
  library(ggplot2); library(patchwork); library(tibble)
  library(tidytext)
  # fallback deps only if DO_FALLBACK <- TRUE
  library(caret); library(ranger); library(yardstick); library(vip)
})

# ---------- PARAMS ----------
DO_FALLBACK <- FALSE     # TRUE to recompute quick PFI if CSVs missing
TOP_K       <- 12        # bars per panel
DROP_SPEI6  <- FALSE     # TRUE to hide SPEI-6 in plots

# ---------- PATHS ----------
base_dir  <- "/Users/marco/Library/CloudStorage/Dropbox/estcena/scripts/farzad"
csv_means <- file.path(base_dir, "data", "weather_3day_means_per_fire_filtered_anom.csv")
out_dir   <- file.path(base_dir, "FIGS"); dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

# ---------- TRANSITIONS (panel order) ----------
ORDER_TAGS <- c("Medium_to_Large",
                "Medium_to_VeryLarge",
                "Medium_to_Mega",
                "Large_to_VeryLarge",
                "Large_to_Mega",
                "VeryLarge_to_Mega")

PANEL_TITLES <- c(
  "Medium_to_Large"      = "Medium \u2192 Large",
  "Medium_to_VeryLarge"  = "Medium \u2192 Very Large",
  "Medium_to_Mega"       = "Medium \u2192 Mega",
  "Large_to_VeryLarge"   = "Large \u2192 Very Large",
  "Large_to_Mega"        = "Large \u2192 Mega",
  "VeryLarge_to_Mega"    = "Very Large \u2192 Mega"
)

# ---------- HELPERS ----------
# Short acronyms (same as Fig. 2)
short_var <- function(x){
  y <- x
  y <- gsub("^mean3_", "", y, ignore.case = TRUE)
  y <- gsub("_interp$", "", y, ignore.case = TRUE)
  y_low <- tolower(y)
  
  dplyr::case_when(
    y_low %in% c("t2m","tas") ~ "AT",
    y_low == "lst_day"        ~ "LSTd",
    y_low == "lst_night"      ~ "LSTn",
    y_low == "rh"             ~ "RH",
    y_low == "wind_speed"     ~ "WS",
    y_low == "tp"             ~ "PR",
    y_low == "smi"            ~ "SMI",
    y_low == "ndvi"           ~ "NDVI",
    y_low == "lai"            ~ "LAI",
    y_low == "spei01"         ~ "SPEI-1",
    y_low == "spei03"         ~ "SPEI-3",
    y_low == "spei06"         ~ "SPEI-6",
    TRUE ~ toupper(y)
  )
}

# Slow = SMI, NDVI, LAI, all SPEIs; Fast = everything else
speed_group_of <- function(v){
  vlow <- tolower(v)
  is_slow <- grepl("^spei", vlow) |
    grepl("(^|_)smi($|_)|(^|_)ndvi($|_)|(^|_)lai($|_)", vlow)
  ifelse(is_slow, "Slow", "Fast")
}

pal_speed <- c(Fast = "#8B0000", Slow = "#555555")  # dark red, dark grey

# ---------- LOAD SAVED PFI ----------
load_saved_pfi <- function(){
  comb <- file.path(out_dir, "PFI_CI_ALL_ALL_PFI.csv")
  if (file.exists(comb)) {
    df <- readr::read_csv(comb, show_col_types = FALSE)
    if (nrow(df)) return(df)
  }
  pieces <- list()
  cand <- list.files(out_dir, pattern="^PFI_CI_.*\\.csv$", full.names=TRUE)
  if (length(cand)) {
    for (p in cand) {
      dd <- try(readr::read_csv(p, show_col_types = FALSE), silent=TRUE)
      if (!inherits(dd, "try-error") && nrow(dd)) pieces[[length(pieces)+1]] <- dd
    }
  }
  if (length(pieces)) dplyr::bind_rows(pieces) else NULL
}

# ---------- FALLBACK (optional quick recomputation) ----------
compute_pfi_fast <- function(){
  stopifnot(file.exists(csv_means))
  df <- readr::read_csv(csv_means, show_col_types = FALSE)
  
  mean3_all <- names(df)[grepl("^mean3_", names(df))]
  spei_all  <- names(df)[grepl("^spei0[136]_interp$", names(df))]
  all_pred  <- c(mean3_all, spei_all)
  
  drop_patterns <- "(^|_)vpd($|_)|(^|_)d2m($|_)|dew"
  pred_vars <- all_pred[!grepl(drop_patterns, all_pred, ignore.case = TRUE)]
  
  # quick size classes
  df <- df |>
    dplyr::mutate(size_class = dplyr::case_when(
      area_ha >= 30   & area_ha < 100    ~ "Medium",
      area_ha >= 100  & area_ha < 1000   ~ "Large",
      area_ha >= 1000 & area_ha < 10000  ~ "Very Large",
      area_ha >= 10000                   ~ "Mega",
      TRUE ~ NA_character_
    )) |>
    dplyr::filter(!is.na(size_class)) |>
    dplyr::filter(dplyr::if_all(dplyr::all_of(pred_vars), ~ is.finite(.) & !is.na(.)))
  
  auc_metric <- function(truth, estimate){
    yardstick::roc_auc_vec(truth = truth, estimate = estimate, event_level = "first")
  }
  pred_fun <- function(object, newdata){
    pr <- predict(object, data = newdata)$predictions
    pr[,"L1"]
  }
  
  set.seed(2025); B <- 100; ntree <- 300; tr_frac <- 0.8
  
  transitions <- tibble::tribble(
    ~lower, ~upper, ~tag,
    "Medium", "Large",      "Medium_to_Large",
    "Medium", "Very Large", "Medium_to_VeryLarge",
    "Medium", "Mega",       "Medium_to_Mega",
    "Large",  "Very Large", "Large_to_VeryLarge",
    "Large",  "Mega",       "Large_to_Mega",
    "Very Large", "Mega",   "VeryLarge_to_Mega"
  )
  
  out_list <- list()
  
  for (i in seq_len(nrow(transitions))){
    lo <- transitions$lower[i]; up <- transitions$upper[i]
    tg <- transitions$tag[i]
    
    sub <- df |>
      dplyr::filter(size_class %in% c(lo, up)) |>
      dplyr::select(dplyr::all_of(pred_vars), size_class) |>
      dplyr::mutate(target = factor(ifelse(size_class==up,"L1","L0")),
                    .keep="unused")
    sub$target <- stats::relevel(sub$target, ref="L1")
    if (length(unique(sub$target)) < 2L || ncol(sub) < 3L) next
    
    all_vi <- vector("list", B)
    for (b in seq_len(B)) {
      idx_tr <- caret::createDataPartition(sub$target, p=tr_frac, list=FALSE)
      dtr <- sub[idx_tr,,drop=FALSE]; dte <- sub[-idx_tr,,drop=FALSE]
      if (length(unique(dte$target)) < 2L) next
      
      feats <- setdiff(colnames(dtr), "target"); if (!length(feats)) next
      zv <- vapply(dtr[,feats,drop=FALSE], function(x) length(unique(x))<=1L, logical(1))
      if (any(zv)) {
        feats <- feats[!zv]; if (!length(feats)) next
        dtr <- dtr[,c("target",feats),drop=FALSE]; dte <- dte[,c("target",feats),drop=FALSE]
      }
      
      tab <- table(dtr$target); w <- as.numeric(sum(tab)/(2*tab)); names(w) <- names(tab)
      
      fit <- ranger::ranger(
        target ~ ., data=dtr, mtry=1, num.trees=ntree,
        probability=TRUE, classification=TRUE,
        case.weights=w[as.character(dtr$target)],
        respect.unordered.factors=TRUE,
        seed=2025+b
      )
      
      vi_b <- try(vip::vi_permute(
        object=fit,
        train=dte[,feats,drop=FALSE],
        target=dte$target,
        metric=auc_metric,
        pred_wrapper=pred_fun,
        nsim=1, parallel=FALSE, smaller_is_better=FALSE
      ), silent=TRUE)
      
      if (!inherits(vi_b,"try-error") && !is.null(vi_b) && nrow(vi_b)) {
        vi_b$run <- b; all_vi[[b]] <- vi_b
      }
    }
    
    vi_list <- Filter(Negate(is.null), all_vi)
    if (!length(vi_list)) next
    
    pfi_summ <- dplyr::bind_rows(vi_list) |>
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
        speed_group    = speed_group_of(variable),
        variable_short = short_var(variable),
        transition     = tg
      ) |>
      dplyr::arrange(dplyr::desc(mean_imp))
    
    out_list[[tg]] <- pfi_summ
  }
  
  dplyr::bind_rows(out_list)
}

# ---------- LOAD OR COMPUTE ----------
pfi_raw <- load_saved_pfi()
if (is.null(pfi_raw) || !nrow(pfi_raw)) {
  if (!DO_FALLBACK) stop("❌ No PFI CSVs found. Set DO_FALLBACK <- TRUE to recompute quickly.")
  message("⚠️ PFI CSVs not found. Running quick fallback...")
  pfi_raw <- compute_pfi_fast()
  if (is.null(pfi_raw) || !nrow(pfi_raw)) stop("❌ No PFI data available.")
}

# ---------- NORMALIZE COLUMNS ----------
if (!"variable" %in% names(pfi_raw) && "variable_pretty" %in% names(pfi_raw)) {
  pfi_raw$variable <- pfi_raw$variable_pretty
}
if (!"speed_group" %in% names(pfi_raw)) pfi_raw$speed_group <- speed_group_of(pfi_raw$variable)
if (!"variable_short" %in% names(pfi_raw)) pfi_raw$variable_short <- short_var(pfi_raw$variable)

# Drop SPEI-6 if requested
if (DROP_SPEI6) {
  pfi_raw <- pfi_raw %>%
    dplyr::filter(
      !grepl("\\bSPEI-6\\b", variable_short, ignore.case = TRUE),
      !grepl("spei06", variable, ignore.case = TRUE)
    )
}

# Normalize transition keys
pfi_raw <- pfi_raw %>%
  dplyr::mutate(
    transition = as.character(transition),
    transition = gsub("\u2192", "->", transition, fixed = TRUE),
    transition = gsub("[[:space:]]+", "_", transition),
    transition = gsub("_PFI$", "", transition)
  )

# Split by transition
pfi_list <- pfi_raw %>%
  dplyr::arrange(dplyr::desc(mean_imp)) %>%
  { split(., .$transition) }

# ---------- PANEL PLOTTER (ranked per panel; short y labels; font 12; y labels ALWAYS shown) ----------
plot_pfi_panel <- function(summ_df, panel_title){
  if (is.null(summ_df) || !nrow(summ_df)) return(ggplot2::ggplot() + ggplot2::theme_void())
  
  sub <- summ_df |> dplyr::slice_head(n = TOP_K)
  
  ggplot2::ggplot(
    sub,
    ggplot2::aes(
      x = mean_imp,
      y = tidytext::reorder_within(variable_short, mean_imp, transition),
      fill = speed_group
    )
  ) +
    ggplot2::geom_col(width = 0.8) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = p2.5, xmax = p97.5), height = 0) +
    tidytext::scale_y_reordered() +
    ggplot2::scale_fill_manual(values = pal_speed) +
    ggplot2::labs(
      title = panel_title,
      x = "Permutation importance (AUC decrease)",
      y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.title = ggplot2::element_blank(),
      legend.text  = ggplot2::element_text(size = 12),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor   = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(6, 8, 6, 6),
      plot.title   = ggplot2::element_text(size = 12),
      axis.title.x = ggplot2::element_text(size = 12),
      axis.text.x  = ggplot2::element_text(size = 12),
      axis.text.y  = ggplot2::element_text(size = 12)
    )
}

# ---------- BUILD PANELS (2×3) ----------
plots <- vector("list", length(ORDER_TAGS))
letters_idx <- 0

for (i in seq_along(ORDER_TAGS)) {
  tag  <- ORDER_TAGS[i]
  df_i <- pfi_list[[tag]]
  
  if (!is.null(df_i) && nrow(df_i)) {
    letters_idx <- letters_idx + 1
    title_i <- paste0(letters[letters_idx], ") ", PANEL_TITLES[[tag]])
    plots[[i]] <- plot_pfi_panel(df_i, title_i)
  } else {
    plots[[i]] <- ggplot2::ggplot() + ggplot2::theme_void()
  }
}

# ---------- LAYOUT 2×3 ----------
g_2x3 <- (plots[[1]] + plots[[2]] + plots[[3]]) /
  (plots[[4]] + plots[[5]] + plots[[6]]) +
  patchwork::plot_layout(guides = "collect") &
  ggplot2::theme(legend.position = "bottom")

# ---------- SAVE ----------
ggplot2::ggsave(file.path(out_dir, "PFI_rankings_2x3_RF_SPEED2COLORS_anom.pdf"),
                g_2x3, width = 14, height = 8)
ggplot2::ggsave(file.path(out_dir, "PFI_rankings_2x3_RF_SPEED2COLORS_anom.png"),
                g_2x3, width = 14, height = 8, dpi = 300)

message("✅ Plots saved in: ", out_dir,
        "\n - PFI_rankings_2x3_RF_SPEED2COLORS_anom.[pdf|png]")
