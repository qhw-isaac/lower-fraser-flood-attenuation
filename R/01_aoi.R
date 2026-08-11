# ==============================================================================
# 01_aoi.R: downstream AOI + regional-district reference layers
# ------------------------------------------------------------------------------
# downstream_aoi = where downstream flood risk is evaluated (MVRD ∪ FVRD).
#
# Outputs (data/processed/):
#   01_aoi_downstream.gpkg
#   01_mvrd.gpkg, 01_fvrd.gpkg
# ==============================================================================

source(here::here("R", "00_setup.R"))

rds <- bcmaps::regional_districts(ask = FALSE) |>
  sf::st_transform(PROJECT_CRS)

mvrd <- dplyr::filter(rds, ADMIN_AREA_NAME == "Metro Vancouver Regional District")
fvrd <- dplyr::filter(rds, ADMIN_AREA_NAME == "Fraser Valley Regional District")

# the two districts together are the downstream demand area
lower_mainland <- sf::st_union(
  c(sf::st_geometry(mvrd), sf::st_geometry(fvrd))
)

aoi_downstream <- sf::st_sf(role = "downstream", geometry = lower_mainland,
                            crs = PROJECT_CRS)

p <- paths()$processed

# write one polygon layer, reporting its area as a sanity check
file_save <- function(x, name) {
  f <- file.path(p, paste0(name, ".gpkg"))
  sf::st_write(x, f, delete_dsn = TRUE, quiet = TRUE)
  message("  ✓ ", basename(f),
          "  (area = ", round(sum(as.numeric(sf::st_area(x))) / 1e6, 1), " km²)")
}

file_save(mvrd, "01_mvrd")
file_save(fvrd, "01_fvrd")
file_save(aoi_downstream, "01_aoi_downstream")