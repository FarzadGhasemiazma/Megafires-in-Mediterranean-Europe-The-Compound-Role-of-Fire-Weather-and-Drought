library(sf)
library(dplyr)
library(stringr)

# Ruta al fichero
base_dir  <- "/Users/marco/Documents/dati/effis/MESOGEOS_DOMAIN"
gpkg_file <- file.path(base_dir, "all_fires_mesogeos.gpkg")

# Leer todos los incendios
fires <- st_read(gpkg_file, quiet = TRUE)

# (opción 1) Filtrar por nombre de país tal y como viene en el atributo `country`
countries_keep <- c(
  "Albania",
  "Bosnia and Herzegovina",
  "Bulgaria",
  "Croatia",
  "España",
  "France",
  "Greece",
  "Hrvatska",
  "Italia",
  "Kosovo",
  "Macedonia",
  "Montenegro",
  "Portugal",
  "Serbia",
  "Turkey",
  "Ελλάδα",
  "България"
)

fires_sel <- fires %>%
  mutate(
    country = str_trim(country),           # por si hay espacios raros
    year    = as.integer(year),
    month   = as.integer(month),
    area_ha = as.numeric(area_ha)
  ) %>%
  filter(
    area_ha > 30,               # > 30 ha
    month %in% 6:9,             # JJAS (Jun–Sep)
    year >= 2008, year <= 2022,
    country %in% countries_keep
  )

# (opcional) echar un ojo rápido
table(fires_sel$country)
summary(fires_sel$area_ha)

# Guardar el resultado si quieres
out_file <- file.path(base_dir,
                      "all_fires_mesogeos_JJAS_gt30ha_2008_2022_med_countries.gpkg")

st_write(fires_sel, out_file, delete_dsn = TRUE)
