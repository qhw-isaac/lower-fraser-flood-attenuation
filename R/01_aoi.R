# ==============================================================================
# 01_aoi.R — create project AOIs (Area of Interest)
# ------------------------------------------------------------------------------
# This project involves three sub_AOIs:
#
#   outcome_aoi    — concerns final deliverables. spans MVRD ∪ FVRD
#                    (approx. the Hope-to-coast Lower Mainland strip)
#   downstream_aoi — polygon used to calculate downstream risks
#                    spans the same area as outcome_aoi
#   upstream_aoi   — where run-off retention is calculated
#
# Inputs:
#   bcmaps::regional_districts(); MVRD + FVRD + SLRD polygons
#
# Outputs (data/processed/):
#   01_aoi_outcome.gpkg                MVRD ∪ FVRD
#   01_aoi_downstream.gpkg             MVRD ∪ FVRD
#   01_aoi_upstream.gpkg               MVRD ∪ FVRD ∪ SLRD (buffer)
#                                      (placeholder; refined by R/07_subbasins.R)
#
#   individual districts are kept for accessibility
#    - 01_mvrd.gpkg, 01_fvrd.gpkg, 01_slrd.gpkg
# ==============================================================================

# run setup script with libraries
source(here::here("R", "00_setup.R"))

# downloads bcmaps polygons and transforms it to project CRS (EPSG:3005)
rds <- bcmaps::regional_districts() |>
  sf::st_transform(PROJECT_CRS)

# pulls three regional districts of interest
mvrd <- dplyr::filter(rds, ADMIN_AREA_NAME == "Metro Vancouver Regional District")
fvrd <- dplyr::filter(rds, ADMIN_AREA_NAME == "Fraser Valley Regional District")
slrd <- dplyr::filter(rds, ADMIN_AREA_NAME == "Squamish-Lillooet Regional District")

# merge mvrd and fvrd into lower_mainland
lower_mainland <- sf::st_union(
  c(sf::st_geometry(mvrd), sf::st_geometry(fvrd))
)

# set both aoi_outcome / aoi_downstream
aoi_outcome    <- sf::st_sf(role = "outcome",    geometry = lower_mainland, crs = PROJECT_CRS)
aoi_downstream <- sf::st_sf(role = "downstream", geometry = lower_mainland, crs = PROJECT_CRS)

# placeholder upstream aoi with buffer to address future alignment issues with borders
# added slrd as it drains into the lower mainland
aoi_upstream <- sf::st_union(
  c(sf::st_geometry(mvrd), sf::st_geometry(fvrd), sf::st_geometry(slrd))
) |>
  sf::st_buffer(EDGE_BUFFER_M) |>
  sf::st_sf(role = "upstream_placeholder", crs = PROJECT_CRS)

# location to output processed data files
p <- paths()$processed

# take x and save as gpkg file
file_save <- function(x, name) {
  f <- file.path(p, paste0(name, ".gpkg"))
  sf::st_write(x, f, delete_dsn = TRUE, quiet = TRUE)
  message("  ✓ ", basename(f),
          "  (area = ", round(sum(as.numeric(sf::st_area(x))) / 1e6, 1), " km²)")
}

file_save(mvrd, "01_mvrd")
file_save(fvrd, "01_fvrd")
file_save(slrd, "01_slrd")
file_save(aoi_outcome, "01_aoi_outcome")
file_save(aoi_downstream, "01_aoi_downstream")
file_save(aoi_upstream, "01_aoi_upstream")