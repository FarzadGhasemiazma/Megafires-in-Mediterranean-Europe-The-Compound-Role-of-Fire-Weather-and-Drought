# =============================
# CLEAN START & PACKAGES
# =============================
rm(list = ls()); graphics.off(); cat("\f"); options(scipen = 999)

pkgs <- c("readr","dplyr","lubridate","sf","ggplot2",
          "rnaturalearth","rnaturalearthdata","tidyr")
to_install <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
if (length(to_install)) install.packages(to_install, quiet = TRUE)
invisible(lapply(pkgs, library, character.only = TRUE))

# =============================
# PATHS
# =============================
base_dir  <- "/Users/marco/Library/CloudStorage/Dropbox/estcena/scripts/farzad"

csv_means <- file.path(base_dir, "data",
                       "weather_3day_means_per_fire_filtered_anom.csv")

out_dir   <- file.path(base_dir, "FIGS")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# =============================
# READ
# =============================
df <- readr::read_csv(csv_means, show_col_types = FALSE)

# controlli minimi con i nuovi nomi
stopifnot(all(c("fire_id","date_start","area_ha","x","y") %in% names(df)))

# se esiste la vecchia size_class (tipo "Medium [30–100)") la salvo per riferimento
if ("size_class" %in% names(df)) {
  df <- df |>
    dplyr::rename(size_class_label = size_class)
} else {
  df$size_class_label <- NA_character_
}

# =============================
# CLASSI + PALETTE
# =============================
pal_sc <- c("Medium"     = "#9FB3C8",
            "Large"      = "#F4B266",
            "Very Large" = "#E67E22",
            "Mega"       = "#C0392B")

class_order <- c("Medium","Large","Very Large","Mega")

label_base <- c("Medium"     = "Medium",
                "Large"      = "Large",
                "Very Large" = "Very Large",
                "Mega"       = "Mega")

df <- df |>
  dplyr::mutate(
    date = suppressWarnings(lubridate::as_date(date_start)),
    size_class = dplyr::case_when(
      area_ha >=   30 & area_ha <    100 ~ "Medium",
      area_ha >=  100 & area_ha <   1000 ~ "Large",
      area_ha >= 1000 & area_ha <  10000 ~ "Very Large",
      area_ha >= 10000                   ~ "Mega",
      TRUE ~ NA_character_
    ),
    size_class = factor(size_class, levels = class_order)
  ) |>
  dplyr::filter(!is.na(size_class), !is.na(x), !is.na(y))

# --- (i conteggi li puoi lasciare, ma NON li usiamo per la legenda) ---
class_counts <- df |>
  dplyr::count(size_class, name = "n") |>
  tidyr::complete(size_class = factor(class_order, levels = class_order),
                  fill = list(n = 0)) |>
  dplyr::arrange(size_class)

# legenda SENZA numeri
legend_labels_simple <- label_base[class_order]

# =============================
# SF OBJECT + MAP DATA
# =============================
fires_sf <- sf::st_as_sf(df, coords = c("x","y"), crs = 4326)
world    <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

# Bounding box Euro–Med (un po' più grande)
xlim <- c(-10, 37)
ylim <- c(34, 48)


# Dimensioni dei cerchi — più contrastate
size_map <- c("Medium"     = 0.7,
              "Large"      = 1.8,
              "Very Large" = 3.2,
              "Mega"       = 6.0)

# Split per disegnare in ordine
fires_medium    <- dplyr::filter(fires_sf, size_class == "Medium")
fires_large     <- dplyr::filter(fires_sf, size_class == "Large")
fires_verylarge <- dplyr::filter(fires_sf, size_class == "Very Large")
fires_mega      <- dplyr::filter(fires_sf, size_class == "Mega")

# =============================
# PLOT
# =============================
p <- ggplot() +
  geom_sf(data = world, fill = "grey95", color = "grey80", linewidth = 0.2) +
  
  # solo contorni colorati (no riempimento): uso colour=, shape=1
  geom_sf(data = fires_medium,
          aes(colour = size_class),
          shape = 1,
          alpha = 0.5,
          size  = size_map["Medium"],
          show.legend = TRUE) +
  geom_sf(data = fires_large,
          aes(colour = size_class),
          shape = 1,
          alpha = 0.6,
          size  = size_map["Large"],
          show.legend = TRUE) +
  geom_sf(data = fires_verylarge,
          aes(colour = size_class),
          shape = 1,
          alpha = 0.7,
          size  = size_map["Very Large"],
          show.legend = TRUE) +
  geom_sf(data = fires_mega,
          aes(colour = size_class),
          shape = 1,
          alpha = 0.9,
          size  = size_map["Mega"],
          show.legend = TRUE) +
  
  scale_colour_manual(
    values = pal_sc,
    breaks = class_order,
    labels = legend_labels_simple,
    name   = "Fire class"
  ) +
  guides(
    colour = guide_legend(
      override.aes = list(
        size  = as.numeric(size_map[class_order]),
        alpha = 1,
        shape = 1
      )
    )
  ) +
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
  labs(x = "Longitude", y = "Latitude") +
  theme_minimal(base_size = 12) +
  theme(
    # legenda piccola, dentro la mappa in alto a destra
    plot.margin = margin(0, 0, 0, 0, "pt"),
    legend.position      = c(0.98, 0.98),
    legend.justification = c("right", "top"),
    legend.background    = element_rect(fill = "white", colour = "grey80"),
    legend.title         = element_text(size = 12),
    legend.text          = element_text(size = 12),
    panel.grid.minor     = element_blank(),
    panel.grid.major     = element_line(color = "grey90"),
    plot.title           = element_blank(),
    plot.subtitle        = element_blank()
  )

# =============================
# SAVE (un po' più grande)
# =============================
f_png <- file.path(out_dir, "Fires_EuroMed_map_byClass_FILTERED_ANOM_COUNTS_outlineOnly.png")
f_pdf <- file.path(out_dir, "Fires_EuroMed_map_byClass_FILTERED_ANOM_COUNTS_outlineOnly.pdf")
ggsave(f_png, p, width = 9, height = 6, dpi = 300)
ggsave(f_pdf, p, width = 9, height = 6)

cat("Saved:\n", f_png, "\n", f_pdf, "\n", sep = "")

