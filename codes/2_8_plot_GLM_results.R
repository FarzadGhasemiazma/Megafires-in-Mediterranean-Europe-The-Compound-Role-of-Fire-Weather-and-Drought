# ============================================================
# ROC GRID — Best per transizione con filtro β* (no-cross-zero)
# Layout:
#   [a] Med→Lrg      [b] Med→VLrg     [c] Med→Mega
#   [ ] vuoto        [d] Lrg→VLrg     [e] Lrg→Mega
#   [ ] vuoto        [ ] vuoto        [f] VLrg→Mega
# Output:
#   /Users/marco/Documents/output/farzad/FIGS/ROC_GRID_best_filtered.[pdf|png]
# ============================================================

rm(list=ls()); graphics.off(); cat("\f"); options(scipen=999)

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(ggplot2); library(patchwork); library(stringr)
})

out_root  <- "/Users/marco/Documents/output/farzad"
tab_dir   <- file.path(out_root,"OUT_TABLES")
roc_dir   <- file.path(out_root,"ROC_PER_RESAMPLE")
fig_dir   <- file.path(out_root,"FIGS"); dir.create(fig_dir, recursive=TRUE, showWarnings=FALSE)

panel_titles <- c(
  M_to_L   = "Medium \u2192 Large",
  M_to_VL  = "Medium \u2192 Very Large",
  M_to_MG  = "Medium \u2192 Mega",
  L_to_VL  = "Large \u2192 Very Large",
  L_to_MG  = "Large \u2192 Mega",
  VL_to_MG = "Very Large \u2192 Mega"
)

order_keys <- c("M_to_L","M_to_VL","M_to_MG","L_to_VL","L_to_MG","VL_to_MG")

# Abbrevia nomi variabili (adaptado a WS, LSTn, NDVI, SPEI3, etc.)
abbr_var <- function(v){
  v <- trimws(gsub("^mean3_", "", v, ignore.case = TRUE))
  
  # LST Day / Night
  v <- gsub("lst_day",   "LSTd", v, ignore.case = TRUE)
  v <- gsub("lst_night", "LSTn", v, ignore.case = TRUE)
  
  # Viento
  v <- gsub("wind_speed", "WS",   v, ignore.case = TRUE)
  
  # Vegetación / humedad suelo
  v <- gsub("\\bndvi\\b", "NDVI", v, ignore.case = TRUE)
  v <- gsub("\\bsmi\\b",  "SMI",  v, ignore.case = TRUE)
  v <- gsub("\\blai\\b",  "LAI",  v, ignore.case = TRUE)
  
  # RH (por si alguna vez aparece)
  v <- gsub("\\brh\\b",   "RH",   v, ignore.case = TRUE)
  
  # SPEI interpolados o no (por si mezclas)
  v <- gsub("spei0?1(_interp)?", "SPEI1", v, ignore.case = TRUE)
  v <- gsub("spei0?3(_interp)?", "SPEI3", v, ignore.case = TRUE)
  v <- gsub("spei0?6(_interp)?", "SPEI6", v, ignore.case = TRUE)
  
  v
}

# ---------- Carica best filtrati (fallback al best "grezzo") ----------
path_best_f <- file.path(tab_dir,"GLM_best_models_by_transition_filtered.csv")
path_best   <- file.path(tab_dir,"GLM_best_models_by_transition.csv")

if (file.exists(path_best_f) && file.info(path_best_f)$size > 0){
  best_raw <- readr::read_csv(path_best_f, show_col_types = FALSE)
  use_filtered <- TRUE
} else {
  best_raw <- readr::read_csv(path_best, show_col_types = FALSE)
  use_filtered <- FALSE
}

stopifnot(all(c("transition","model","AUC_med","RECALL_med") %in% names(best_raw)))

best <- best_raw %>%
  transmute(
    transition = transition,
    model      = model,
    title      = panel_titles[transition],
    AUC        = AUC_med,
    Recall     = RECALL_med
  ) %>%
  arrange(factor(transition, levels=order_keys))

# ---------- Carica betas (β*) e prepara testo ----------
coef_all <- readr::read_csv(file.path(tab_dir,"GLM_coefficients_summary.csv"), show_col_types = FALSE)

beta_text <- function(tr, mod, max_terms = 6){
  cf <- coef_all %>%
    filter(transition == tr, model == mod) %>%
    filter(term != "(Intercept)") %>%           # por si acaso
    mutate(
      term_abbr = abbr_var(term),
      med_beta  = beta_med,                    # nuevos nombres
      lwr_beta  = beta_L,
      upr_beta  = beta_U
    ) %>%
    arrange(desc(abs(med_beta))) %>%
    slice_head(n = max_terms)
  
  if (!nrow(cf)) return("")
  
  paste0(
    sprintf("%+.2f (%+.2f,%+.2f) · %s",
            cf$med_beta, cf$lwr_beta, cf$upr_beta, cf$term_abbr),
    collapse = "\n"
  )
}

# ---------- Carica ROC CI grid ----------
load_ci <- function(tr, mod){
  tag <- gsub("[^A-Za-z0-9]+","_",mod)
  path <- file.path(roc_dir, paste0("ROC_CIgrid_", tr, "_", tag, ".csv"))
  if (!file.exists(path)) return(NULL)
  readr::read_csv(path, show_col_types = FALSE)
}

# ---------- Un pannello ----------
make_panel <- function(row, letter){
  if (is.null(row) || !nrow(row)) return(patchwork::plot_spacer())
  tr  <- row$transition[1];  mod <- row$model[1]
  AUC <- row$AUC[1];         REC <- row$Recall[1]
  ttl <- row$title[1]
  
  ci <- load_ci(tr, mod)
  if (is.null(ci) || !nrow(ci)) return(patchwork::plot_spacer())
  bet <- beta_text(tr, mod)
  
  ggplot(ci, aes(FPR, TPR_med)) +
    geom_ribbon(aes(ymin = TPR_L, ymax = TPR_U), fill = "grey80", alpha = 0.6) +
    geom_line(linewidth = 1.1, colour = "#b2182b") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    scale_x_continuous(limits = c(0,1)) +
    scale_y_continuous(limits = c(0,1)) +
    labs(
      title = paste0(letter, ") ", ttl),
      x = "False Positive Rate (1 - Specificity)",
      y = "True Positive Rate (Sensitivity)"
    ) +
    annotate("text", x=0.99, y=0.25, hjust=1, vjust=0,
             label=sprintf("AUC = %.2f\nRecall = %.0f%%", AUC, 100*REC),
             size=4.0, color="#b2182b") +
    annotate("text", x=0.99, y=0.02, hjust=1, vjust=0, lineheight=1.05,
             label=bet, size=3.3) +
    coord_equal() +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          plot.title = element_text(face = "bold"))
}

# ---------- Costruisci layout 2×3 ----------
rows <- split(best, factor(best$transition, levels=order_keys))
p1 <- make_panel(rows[["M_to_L"]],   "a")
p2 <- make_panel(rows[["M_to_VL"]],  "b")
p3 <- make_panel(rows[["M_to_MG"]],  "c")
p4 <- make_panel(rows[["L_to_VL"]],  "d")
p5 <- make_panel(rows[["L_to_MG"]],  "e")
p6 <- make_panel(rows[["VL_to_MG"]], "f")

g <- (p1 + p2 + p3) /
  (p4 + p5 + p6) +
  patchwork::plot_layout(guides = "collect")

out_pdf <- file.path(fig_dir,"ROC_GRID_best_filtered.pdf")
out_png <- file.path(fig_dir,"ROC_GRID_best_filtered.png")
ggsave(out_pdf, g, width = 13, height = 12)
ggsave(out_png, g, width = 13, height = 12, dpi = 300)

message("Saved → ", out_pdf, " and .png")
