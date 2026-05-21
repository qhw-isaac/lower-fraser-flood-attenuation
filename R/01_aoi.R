# ==============================================================================
# 01_aoi.R — outcome + downstream AOIs and regional-district reference layers
# ------------------------------------------------------------------------------
# Three sub-AOIs are used downstream:
#
#   outcome_aoi    — final deliverable footprint (MVRD ∪ FVRD)
#   downstream_aoi — where downstream risk is evaluated (MVRD ∪ FVRD)
#   upstream_aoi   — built by R/02_subbasins.R from HydroBASINS, NOT here
#                    (a placeholder upstream AOI based on regional districts
#                    misses Fraser/Thompson basins that drain into the Lower
#                    Mainland from outside MVRD ∪ FVRD ∪ SLRD)
#
# Outputs (data/processed/):
#   01_aoi_outcome.gpkg, 01_aoi_downstream.gpkg
#   01_mvrd.gpkg, 01_fvrd.gpkg, 01_slrd.gpkg   (kept for admin joins)
# ==============================================================================

source(here::here("R", "00_setup.R"))

# fetch BC regional districts and reproject to project CRS
rds <- bcmaps::regional_districts() |>
  sf::st_transform(PROJECT_CRS)

mvrd <- dplyr::filter(rds, ADMIN_AREA_NAME == "Metro Vancouver Regional District")
fvrd <- dplyr::filter(rds, ADMIN_AREA_NAME == "Fraser Valley Regional District")
slrd <- dplyr::filter(rds, ADMIN_AREA_NAME == "Squamish-Lillooet Regional District")

lower_mainland <- sf::st_union(
  c(sf::st_geometry(mvrd), sf::st_geometry(fvrd))
)

aoi_outcome    <- sf::st_sf(role = "outcome",    geometry = lower_mainland, crs = PROJECT_CRS)
aoi_downstream <- sf::st_sf(role = "downstream", geometry = lower_mainland, crs = PROJECT_CRS)

p <- paths()$processed

file_save <- function(x, name) {
  f <- file.path(p, paste0(name, ".gpkg"))
  sf::st_write(x, f, delete_dsn = TRUE, quiet = TRUE)
  message("  ✓ ", basename(f),
          "  (area = ", round(sum(as.numeric(sf::st_area(x))) / 1e6, 1), " km²)")
}

file_save(mvrd, "01_mvrd")
file_save(fvrd, "01_fvrd")
file_save(slrd, "01_slrd")
file_save(aoi_outcome,    "01_aoi_outcome")
file_save(aoi_downstream, "01_aoi_downstream")
