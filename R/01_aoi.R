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

# ---- coast-clipped copy, for drawing only ------------------------------------
# The regional-district boundaries are administrative: MVRD reaches well out
# into the Strait of Georgia and FVRD's western limb runs past Bowen Island.
# Every map is framed on land data, so that marine ground falls outside the
# frame and the demand border leaves the panel. Clip to the coastline instead.
# The census DA boundary file is already coast-clipped, so its dissolved
# footprint is the land base (07_floodplains.R uses it the same way).
#
# Display only. Analysis keeps using 01_aoi_downstream.gpkg, since the demand
# AOI is an administrative unit and clipping it would change what the model
# counts.
aoi_land <- tryCatch({
  da_file <- data_path("da_boundary_2021")
  bbox_wkt <- sf::st_as_text(sf::st_as_sfc(sf::st_bbox(
    sf::st_transform(aoi_downstream, sf::st_crs(sf::st_read(
      da_file, quiet = TRUE,
      query = paste0("SELECT * FROM ",
                     tools::file_path_sans_ext(basename(da_file)), " LIMIT 1")))))))
  land <- sf::st_read(da_file, quiet = TRUE, wkt_filter = bbox_wkt) |>
    sf::st_transform(PROJECT_CRS) |>
    sf::st_make_valid() |>
    sf::st_union()
  sf::st_sf(role = "downstream_land",
            geometry = sf::st_make_valid(
              sf::st_intersection(sf::st_geometry(aoi_downstream), land)),
            crs = PROJECT_CRS)
}, error = function(e) {
  message("  ! coast clip skipped (", conditionMessage(e),
          "); maps will draw the administrative boundary")
  NULL
})

if (!is.null(aoi_land)) file_save(aoi_land, "01_aoi_downstream_land")