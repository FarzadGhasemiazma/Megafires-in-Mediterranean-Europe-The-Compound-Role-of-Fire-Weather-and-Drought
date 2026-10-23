# ============================================================
# JACKKNIFE STABILITY PLOTS (AUC + β)
# Using EXISTING CSV (AUC_med + beta_raw)
# ============================================================

user_lib <- "~/R/libs"
if (!dir.exists(user_lib)) dir.create(user_lib, recursive = TRUE)
.libPaths(c(user_lib, .libPaths()))

rm(list=ls()); graphics.off(); cat("\f"); options(scipen=999)

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

# ------------------------------------------------------------
# PATHS
# ------------------------------------------------------------
root_dir <- "/home/farzad/1_farzad_project/1_Codes/JACKKNIFE_250"
metrics_file <- file.path(root_dir, "JACKKNIFE_250_metrics.csv")
beta_file    <- file.path(root_dir, "JACKKNIFE_250_beta_raw.csv")

fig_dir      <- file.path(root_dir, "FIGS_JACKKNIFE")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------
# LOAD DATA
# ------------------------------------------------------------
MET  <- read_csv(metrics_file, show_col_types = FALSE)
BETA <- read_csv(beta_file, show_col_types = FALSE)

# ------------------------------------------------------------
# PANEL ORDER
# ------------------------------------------------------------
panel_order <- c("M_to_MG", "L_to_MG", "VL_to_MG")

panel_titles <- c(
  M_to_MG  = "Medium → Mega",
  L_to_MG  = "Large → Mega",
  VL_to_MG = "Very Large → Mega"
)

panel_letter <- c("a", "b", "c")

# ------------------------------------------------------------
# COLORS
# ------------------------------------------------------------
col_line   <- "#8B0000"
col_ribbon <- "#e6b3b3"

# ------------------------------------------------------------
# AUC SUMMARY (CI)
# ------------------------------------------------------------
AUC_SUMM <- MET %>%
  group_by(transition, loo_id) %>%
  summarise(
    AUC_mid = median(AUC_med),
    AUC_L   = quantile(AUC_med, 0.025),
    AUC_U   = quantile(AUC_med, 0.975),
    .groups = "drop"
  )

# ------------------------------------------------------------
# BETA SUMMARY (CI)
# ------------------------------------------------------------
BETA_SUMM <- BETA %>%
  group_by(transition, loo_id, term) %>%
  summarise(
    beta_mid = median(beta_raw),
    beta_L   = quantile(beta_raw, 0.025),
    beta_U   = quantile(beta_raw, 0.975),
    .groups = "drop"
  )

# ------------------------------------------------------------
# PLOT FUNCTIONS
# ------------------------------------------------------------
plot_auc_panel <- function(df, title, letter){
  ggplot(df, aes(x = loo_id, y = AUC_mid)) +
    geom_ribbon(aes(ymin = AUC_L, ymax = AUC_U),
                fill = col_ribbon, alpha = 0.55) +
    geom_line(color = col_line, linewidth = 1) +
    labs(title = paste0(letter, ") ", title, "\nAUC Jackknife Stability"),
         x = "Leave-One-Out Mega-fire",
         y = "AUC") +
    theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold")
    )
}

plot_beta_panel <- function(df){
  ggplot(df, aes(x = loo_id, y = beta_mid, color = term)) +
    geom_ribbon(aes(ymin = beta_L, ymax = beta_U, fill = term),
                alpha = 0.30, color = NA) +
    geom_line(linewidth = 1) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(title = "β Coefficient Stability",
         x = "Leave-One-Out Mega-fire",
         y = "β estimate") +
    theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      legend.position = "bottom"
    )
}

# ------------------------------------------------------------
# BUILD PANELS
# ------------------------------------------------------------
panels <- list()

for (i in seq_along(panel_order)) {

  tr   <- panel_order[i]
  ttl  <- panel_titles[[tr]]
  let  <- panel_letter[i]

  df_auc  <- AUC_SUMM  %>% filter(transition == tr)
  df_beta <- BETA_SUMM %>% filter(transition == tr)

  gAUC  <- plot_auc_panel(df_auc, ttl, let)
  gBETA <- plot_beta_panel(df_beta)

  panels[[i]] <- gAUC / gBETA
}

# ------------------------------------------------------------
# COMBINE INTO 1×3 LAYOUT
# ------------------------------------------------------------
g_final <- panels[[1]] | panels[[2]] | panels[[3]] +
  plot_layout(guides = "collect")

# ------------------------------------------------------------
# SAVE
# ------------------------------------------------------------
pdf_out <- file.path(fig_dir, "JACKKNIFE_combined_panels.pdf")
png_out <- file.path(fig_dir, "JACKKNIFE_combined_panels.png")

ggsave(pdf_out, g_final, width = 18, height = 6)
ggsave(png_out, g_final, width = 18, height = 6, dpi = 300)

message("Saved combined jackknife figure:")
message("  ", pdf_out)
message("  ", png_out)
