# ==============================================================================
# 05_dem_slope.R — DEM and slope raster
# ------------------------------------------------------------------------------
# Replace HydroSHEDS with Copernicus GLO-30 (30m native DSM, TanDEM-X SAR)
#
# Inputs (data/processed/):
#   03_lulc_values.tif
#
# Outputs (data/processed/):
#   05_dem.tif
#   05_slope_deg.tif
# ==============================================================================

source(here::here("R", "00_setup.R"))

dem_path <- data_path("dem_glo30", must_exist = TRUE)
message("loading DEM 'dem_glo30' from ", dem_path)
dem <- terra::rast(dem_path) |> align_to_grid(method = "bilinear")

slope_real <- terra::terrain(dem, v = "slope", neighbors = 8, unit = "degrees")
lulc <- terra::rast(file.path(paths()$processed, "03_lulc_values.tif"))
slope_deg <- terra::ifel(is.na(slope_real) & !is.na(lulc), 10, slope_real)
aoi_mask <- terra::rasterize(terra::vect(read_aoi("upstream")), slope_deg, field = 1)
slope_deg <- terra::mask(slope_deg, aoi_mask)

safe_writeRaster(dem, file.path(paths()$processed, "05_dem.tif"))
safe_writeRaster(slope_deg, file.path(paths()$processed, "05_slope_deg.tif"))

# coverage diagnostic
land_mask <- !is.na(lulc)
total_px <- terra::global(land_mask, "sum", na.rm = TRUE)[1, 1]
fill_mask <- is.na(slope_real) & land_mask
fill_px <- terra::global(fill_mask, "sum", na.rm = TRUE)[1, 1]
dem_px <- total_px - fill_px

message(sprintf("GLO-30-derived slope: %10d px  (%5.1f%%)",
                as.integer(dem_px), 100 * dem_px / total_px))

message(sprintf("10° Duarte fallback (Huang ≈ off): %10d px  (%5.1f%%)",
                as.integer(fill_px), 100 * fill_px / total_px))

if (fill_px / total_px > 0.05) {
  message("  ! >5% of slope pixels are still on the 10° fallback. The ",
          "GLO-30 clip probably doesn't fully cover the upstream AOI — ",
          "re-export it from Earth Engine for the full provider footprint.")
}

# ---- QA preview --------------------------------------------------------------
qa_png("05_dem_slope.png", function() {
  op <- graphics::par(mar = c(2, 2, 3, 7))
  on.exit(graphics::par(op), add = TRUE)
  terra::plot(slope_deg, main = "Slope (degrees)",
              col = grDevices::hcl.colors(50, palette = "YlOrRd", rev = FALSE),
              range = c(0, 60))
})

message("✓ 05_dem_slope.R — wrote DEM + slope (deg) from GLO-30")