# =========================================
# DISTRIBUTIONS (overlay) — 12 paneles (3x4)
# Anomalías / z-scores + SPEI
# Ejes X fijos: [-2.5, 2.5]
# =========================================
rm(list=ls()); graphics.off(); cat("\f"); options(scipen=999)

pkgs <- c("dplyr","tidyr","ggplot2","readr","stringr","purrr","scales","patchwork","ggtext")
to_install <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
if(length(to_install)) install.packages(to_install, quiet=TRUE)
invisible(lapply(pkgs, library, character.only=TRUE))

# --- PATHS ---
base_dir  <- "/Users/marco/Library/CloudStorage/Dropbox/estcena/scripts/farzad"
csv_means <- file.path(base_dir, "data", "weather_3day_means_per_fire_filtered_anom.csv")
out_dir   <- file.path(base_dir, "FIGS"); dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

# --- LOAD ---
df <- readr::read_csv(csv_means, show_col_types = FALSE)

# --- SIZE CLASSES (4 niveles) ---
df <- df |>
  dplyr::mutate(
    size_class = dplyr::case_when(
      area_ha >=   30 & area_ha <    100   ~ "Medium",
      area_ha >=  100 & area_ha <   1000   ~ "Large",
      area_ha >= 1000 & area_ha <  10000   ~ "Very Large",
      area_ha >= 10000                     ~ "Mega",
      TRUE ~ NA_character_
    )
  ) |>
  dplyr::filter(!is.na(size_class))

# ======================================================
# VARIABLES + ACRÓNIMOS + ETIQUETAS EJE X
# (mean3_* son z-scores, conv = identidad)
# ======================================================
VAR_DEF <- list(
  # 1ª fila — temperaturas (z-score)
  AT    = list(col="mean3_t2m",        conv=function(x) x, xlab="AT — Air Temperature (z-score)",        type="linear"),
  LSTd  = list(col="mean3_lst_day",    conv=function(x) x, xlab="LSTd — LST Day (z-score)",              type="linear"),
  LSTn  = list(col="mean3_lst_night",  conv=function(x) x, xlab="LSTn — LST Night (z-score)",            type="linear"),
  # 2ª fila — RH, viento, precipitación (z-score)
  RH    = list(col="mean3_rh",         conv=function(x) x, xlab="RH — Relative Humidity (z-score)",      type="linear"),
  WS    = list(col="mean3_wind_speed", conv=function(x) x, xlab="WS — Wind Speed (z-score)",             type="linear"),
  PR    = list(col="mean3_tp",         conv=function(x) x, xlab="PR — Precipitation (z-score)",          type="linear"),
  # 3ª fila — suelo/vegetación (z-score)
  SMI   = list(col="mean3_smi",        conv=function(x) x, xlab="SMI — Soil Moisture Index (z-score)",   type="linear"),
  NDVI  = list(col="mean3_ndvi",       conv=function(x) x, xlab="NDVI — Vegetation Index (z-score)",     type="linear"),
  LAI   = list(col="mean3_lai",        conv=function(x) x, xlab="LAI — Leaf Area Index (z-score)",       type="linear"),
  # 4ª fila — sequía SPEI (índices, no z-score)
  `SPEI-1` = list(col="spei01_interp", conv=function(x) x, xlab="SPEI-1 — 1-month SPEI (–)",             type="linear"),
  `SPEI-3` = list(col="spei03_interp", conv=function(x) x, xlab="SPEI-3 — 3-month SPEI (–)",             type="linear"),
  `SPEI-6` = list(col="spei06_interp", conv=function(x) x, xlab="SPEI-6 — 6-month SPEI (–)",             type="linear")
)

# deja solo las variables cuyas columnas existen
VAR_DEF <- VAR_DEF[vapply(VAR_DEF, function(v) v$col %in% names(df), logical(1))]
stopifnot(length(VAR_DEF) > 0)

# --- paleta clases ---
pal_sc <- c("Medium"="#9FB3C8","Large"="#F4B266","Very Large"="#E67E22","Mega"="#C0392B")
class_levels <- c("Medium","Large","Very Large","Mega")

# --- helpers ---
finite_vec <- function(x) x[is.finite(x)]

build_long <- function(key){
  v <- VAR_DEF[[key]]
  vals <- v$conv(df[[v$col]])
  out <- tibble::tibble(
    size_class = factor(df$size_class, levels = class_levels),
    value      = vals
  )
  out <- dplyr::filter(out, is.finite(value))
  out
}

pstars <- function(p) ifelse(is.na(p),"", ifelse(p<.001,"***", ifelse(p<.01,"**", ifelse(p<.05,"*",""))))

test_adjacent_lines <- function(dat){
  if (nrow(dat) == 0) return("")
  
  pairs <- list(
    "Med→Lrg"    = c("Medium","Large"),
    "Lrg→V–Lrg"  = c("Large","Very Large"),
    "V–Lrg→Mega" = c("Very Large","Mega")
  )
  
  out <- lapply(names(pairs), function(lbl){
    lv <- pairs[[lbl]][1]
    uv <- pairs[[lbl]][2]
    
    d <- dat |>
      dplyr::filter(size_class %in% c(lv, uv)) |>
      dplyr::mutate(cls = ifelse(size_class == uv, "Upper", "Lower"))
    
    if(nrow(d) < 2 || length(unique(d$cls)) < 2) return(NA_character_)
    
    p <- tryCatch(stats::wilcox.test(value ~ cls, data=d, exact=FALSE)$p.value,
                  error=function(e) NA_real_)
    
    stars <- pstars(p)
    if(stars == "") return(NA_character_)
    
    # color of the higher fire-size class
    col <- pal_sc[[uv]]
    
    # return colored html
    paste0("<span style='color:", col, ";'>", lbl, stars, "</span>")
  })
  
  out <- unlist(out)
  out <- out[!is.na(out) & out != ""]
  
  if(length(out) == 0) "" else paste(out, collapse = "<br>")
}

# --- escala x fija: [-2.5, 2.5] para TODAS las variables ---
scale_x_for <- function(key, data_vec){
  ggplot2::scale_x_continuous(
    limits = c(-4, 4),
    expand = ggplot2::expansion(mult = c(0, 0.02))
  )
}

theme_panel <- ggplot2::theme_minimal(base_size=12) +
  ggplot2::theme(
    plot.title = ggtext::element_markdown(
      size = 12,
      hjust = 0,
      margin = ggplot2::margin(t = 10, r = 0, b = 8, l = 0),
      lineheight = 1.1
    ),
    axis.text    = ggplot2::element_text(size = 12),
    axis.title.x = ggplot2::element_text(size = 12, margin=ggplot2::margin(t=7)),
    axis.title.y = ggplot2::element_text(size = 12, margin=ggplot2::margin(r=7)),
    plot.margin  = ggplot2::margin(t=14, r=8, b=10, l=12)
  )

make_panel <- function(key, tag_letter){
  v   <- VAR_DEF[[key]]
  dat <- build_long(key)
  if (nrow(dat) == 0) {
    # panel placeholder si no hay datos finitos
    return(
      ggplot2::ggplot() +
        ggplot2::annotate("text", x=0.5, y=0.5,
                          label=paste0(tag_letter, ") ", key, "\nNo finite data"),
                          size=3.5, hjust=0.5, vjust=0.5) +
        ggplot2::theme_void() +
        ggplot2::theme(plot.margin = ggplot2::margin(t=14, r=8, b=10, l=12))
    )
  }
  tests_txt <- test_adjacent_lines(dat)
  title_txt <- if(nchar(tests_txt)==0) {
    paste0("**", tag_letter, ")**")
  } else {
    paste0("**", tag_letter, ")**<br>", tests_txt)
  }
  
  ggplot2::ggplot(dat, ggplot2::aes(x=value, color=size_class, fill=size_class)) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.5) +
    ggplot2::stat_density(ggplot2::aes(y=after_stat(density)),
                          adjust=1.15, position="identity",
                          alpha=0.22, linewidth=0) +
    ggplot2::stat_density(ggplot2::aes(y=after_stat(density)),
                          adjust=1.15, position="identity",
                          alpha=1, linewidth=0.6, fill=NA) +
    scale_x_for(key, dat$value) +
    ggplot2::scale_y_continuous(
      limits = c(0, NA),
      breaks = scales::pretty_breaks(n = 3),   # or 4 if you prefer
      expand = ggplot2::expansion(mult = c(0, 0.08))
    ) +
    ggplot2::scale_color_manual(values=pal_sc, name=NULL) +
    ggplot2::scale_fill_manual(values=pal_sc,  name=NULL) +
    ggplot2::guides(
      color = ggplot2::guide_legend(override.aes = list(fill = NA, linewidth = 2),
                                    nrow = 1, byrow = TRUE),
      fill  = "none"
    ) +
    ggplot2::labs(title = title_txt, x = v$xlab, y = "Density") +
    theme_panel
}

# Orden de 12 claves en 4 filas × 3 columnas
order_keys <- c("AT","LSTd","LSTn",
                "RH","WS","PR",
                "SMI","NDVI","LAI",
                "SPEI-1","SPEI-3","SPEI-6")
order_keys <- order_keys[order_keys %in% names(VAR_DEF)]  # por si faltara alguna col
panel_tags <- letters[seq_along(order_keys)]

plots <- purrr::imap(order_keys, ~ make_panel(.x, panel_tags[.y]))

# --- patchwork 3x4, leyenda arriba ---
ncol_out <- 3
g_panels <- patchwork::wrap_plots(plots, ncol = ncol_out)

g_all <- (patchwork::guide_area() / g_panels) +
  patchwork::plot_layout(
    heights = c(0.08, 1),   # first row reserved for legend
    guides  = "collect"
  ) &
  ggplot2::theme(
    legend.position = "top",
    legend.title    = ggplot2::element_blank(),
    legend.text     = ggplot2::element_text(size = 12),
    plot.margin     = ggplot2::margin(t=16, r=10, b=10, l=10)
  )

# --- export ---
f_pdf <- file.path(out_dir, "Distributions_4Classes_12Vars_overlay_PATCHWORK_anom_fixedx.pdf")
f_png <- file.path(out_dir, "Distributions_4Classes_12Vars_overlay_PATCHWORK_anom_fixedx.png")
ggplot2::ggsave(f_pdf, g_all, width=24, height=28, units="cm")
ggplot2::ggsave(f_png, g_all, width=24, height=28, units="cm", dpi=600)

message("✅ Guardados:\n - ", f_pdf, "\n - ", f_png)
