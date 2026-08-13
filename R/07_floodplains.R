# ==============================================================================
# 07_floodplains.R: floodplain raster used to mask "demand" pixels
# ------------------------------------------------------------------------------
# Rasterizes the NHC Lower Fraser 2D flood extent onto the project grid. The
# 1894 flood of record is the natural floodplain footprint the dikes defend; the
# Historic 1% AEP extent (residual risk with dikes in place) is registered as
# nhc_floodplain_1aep for a later sensitivity run.
#
# Scope: the NHC model covers the Fraser freshet corridor, Hope to the Strait of
# Georgia. Floodplains outside it (Serpentine/Nicomekl, coastal storm surge,
# small tributaries) carry no demand here. The retired national source (Mohanty
# & Simonovic CMIP6 100-yr) is coarser and is not a fallback.
#
# Inputs (data/processed/):
#   03_lulc_values.tif
#   data_path("nhc_floodplain_1894")   flood-extent polygons (EPSG:26910)
#
# Outputs (data/processed/):
#   07_floodplain.tif
# ==============================================================================

source(here::here("R", "00_setup.R"))

lulc_grid <- terra::rast(file.path(paths()$processed, "03_lulc_values.tif"))

fp_poly <- sf::st_read(data_path("nhc_floodplain_1894"), quiet = TRUE) |>
  sf::st_transform(PROJECT_CRS)

fp <- terra::rasterize(terra::vect(fp_poly), lulc_grid, field = 1, background = 0)

# punch out mapped WATER BODIES (FWA lakes + double-line rivers), not NALCMS
# water pixels: punching pixels turned every ditch and slough into a hole,
# 438 km² of speckle. In-floodplain ponds stay floodplain, since they flood
# along with the land around them.
downstream_aoi <- read_aoi("downstream")
waterbodies <- tryCatch({
  dom <- sf::st_as_sfc(sf::st_bbox(downstream_aoi))
  lk <- bcdc_query_geodata("WHSE_BASEMAPPING.FWA_LAKES_POLY") |>
    filter(INTERSECTS(dom)) |> collect()
  rv <- bcdc_query_geodata("WHSE_BASEMAPPING.FWA_RIVERS_POLY") |>
    filter(INTERSECTS(dom)) |> collect()
  rbind(lk["geometry"], rv["geometry"]) |> sf::st_transform(PROJECT_CRS)
}, error = function(e) {
  message("  ! FWA water-body pull failed (", conditionMessage(e),
          "). Falling back to masking all NALCMS water")
  NULL
})
if (!is.null(waterbodies)) {
  fp <- terra::mask(fp, terra::rasterize(terra::vect(waterbodies), fp, field = 1),
                    inverse = TRUE)
} else {
  fp <- terra::mask(fp, lulc_grid)
}

# clip to the coastline. NALCMS maps the Sturgeon/Roberts Bank intertidal marsh
# as shrubland, so the LULC mask alone keeps foreshore that reads as ocean. The
# census DA boundary file is already coast-clipped, so it gives the land base.
da_file <- data_path("da_boundary_2021")
da_crs <- sf::st_crs(sf::st_read(da_file, quiet = TRUE,
  query = paste0("SELECT * FROM ", tools::file_path_sans_ext(basename(da_file)),
                 " LIMIT 1")))
bbox_wkt <- sf::st_as_text(sf::st_as_sfc(sf::st_bbox(
  sf::st_transform(downstream_aoi, da_crs))))
land <- sf::st_read(da_file, quiet = TRUE, wkt_filter = bbox_wkt) |>
  sf::st_transform(PROJECT_CRS)
fp <- terra::mask(fp, terra::rasterize(terra::vect(land), fp, field = 1))

# keep only floodplain inside the downstream AOI (where demand is evaluated)
fp <- terra::mask(fp, terra::rasterize(terra::vect(downstream_aoi), fp, field = 1))

safe_writeRaster(fp, file.path(paths()$processed, "07_floodplain.tif"))

# ---- QA preview --------------------------------------------------------------
# the floodplain extent that carries all downstream demand, inside the demand AOI
qa_png("07_floodplain.png", function() {
  # trim() frames on the floodplain, which stops at the coast, so open the frame
  # wide enough to hold the coast-clipped demand boundary too
  aoi_disp <- aoi_display()
  fp_trim <- terra::trim(fp)
  view <- span_bbox(fp_trim, aoi_disp)
  terra::plot(fp_trim, type = "classes", levels = c("dry", "floodplain"),
              col = c("grey90", "#3182bd"), axes = TRUE,
              xlim = view$xlim, ylim = view$ylim,
              main = "Modelled floodplain inside the demand area")
  plot(sf::st_geometry(aoi_disp), add = TRUE, border = "black", lwd = 1)
})

fp_km2 <- terra::global(fp, "sum", na.rm = TRUE)[1, 1] *
  terra::xres(fp) * terra::yres(fp) * 1e-6
message(sprintf(
  "✓ 07_floodplains.R: wrote floodplain raster (NHC 1894 flood of record, %.0f km² in demand AOI)",
  fp_km2))
