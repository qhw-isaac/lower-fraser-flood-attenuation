# ==============================================================================
# 07_floodplains.R: floodplain raster used to mask "demand" pixels
# ------------------------------------------------------------------------------
# Creates the floodplain raster from the JRC/CEMS global river flood hazard map 
# (100-year return period)
#
# Inputs (data/processed/):
#   03_lulc_values.tif
#   data_path("jrc_floodmap_rp100")  RP100 water depth, m (EPSG:3005, 90 m)
#
# Outputs (data/processed/):
#   07_floodplain.tif
# ==============================================================================

source(here::here("R", "00_setup.R"))

lulc_grid <- terra::rast(file.path(paths()$processed, "03_lulc_values.tif"))

# depth -> 0/1, then resample to the 30 m grid
depth <- terra::rast(data_path("jrc_floodmap_rp100"))
fp <- terra::ifel(is.na(depth) | depth <= FLOOD_DEPTH_MIN_M, 0, 1)
fp <- align_to_grid(fp, method = "near", template = lulc_grid)

# punch out mapped WATER BODIES (FWA lakes + double-line rivers)
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

# clip to the coastline
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

# minimum mapping unit
mmu_m2 <- terra::xres(depth) * terra::yres(depth)
clumps <- terra::patches(terra::ifel(fp == 1, 1, NA), directions = 8,
                         zeroAsNA = TRUE)
sz <- terra::freq(clumps)
tiny <- sz$value[sz$count * terra::xres(fp) * terra::yres(fp) < mmu_m2]
if (length(tiny)) {
  keep <- terra::subst(clumps, from = tiny, to = NA)
  fp <- terra::ifel(!is.na(fp) & is.na(keep), 0, fp)
}
message(sprintf("  · sieved %d clumps below the %.2f ha minimum mapping unit",
                length(tiny), mmu_m2 / 1e4))

safe_writeRaster(fp, file.path(paths()$processed, "07_floodplain.tif"))

# ---- QA preview --------------------------------------------------------------
# the floodplain extent that intersects the demand AOI
qa_png("07_floodplain.png", function() {
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
  "✓ 07_floodplains.R: wrote floodplain raster (JRC/CEMS RP100 river flood hazard, %.0f km² in demand AOI)",
  fp_km2))
