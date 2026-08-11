# ==============================================================================
# 05_dem_slope.R: DEM and slope raster
# ------------------------------------------------------------------------------
# Elevation from Copernicus GLO-30 (30 m native DSM, TanDEM-X SAR), which
# replaced HydroSHEDS. Where GLO-30 has no data but land cover does, slope falls
# back to 10°, below the 0.14 threshold at which 08 applies Huang.
#
# Inputs (data/processed/):
#   03_lulc_values.tif
#
# Outputs (data/processed/):
#   05_dem.tif
#   05_slope_deg.tif
# ==============================================================================

source(here::here("R", "00_setup.R"))

dem <- terra::rast(data_path("dem_glo30")) |> align_to_grid(method = "bilinear")

slope_real <- terra::terrain(dem, v = "slope", neighbors = 8, unit = "degrees")
lulc <- terra::rast(file.path(paths()$processed, "03_lulc_values.tif"))
slope_deg <- terra::ifel(is.na(slope_real) & !is.na(lulc), 10, slope_real)
aoi_mask <- terra::rasterize(terra::vect(read_aoi("upstream")), slope_deg, field = 1)
slope_deg <- terra::mask(slope_deg, aoi_mask)
dem <- terra::mask(dem, aoi_mask)

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

message(sprintf("10° flat-slope fallback (Huang ≈ off): %10d px  (%5.1f%%)",
                as.integer(fill_px), 100 * fill_px / total_px))


# ---- QA preview --------------------------------------------------------------
# slope should follow the terrain (i.e. steep in the mountains and flat otherwise)
qa_png("05_dem_slope.png", function() {
  terra::plot(slope_deg, col = grDevices::hcl.colors(50, "YlOrRd"),
              range = c(0, 60), axes = TRUE, main = "")
})

message("✓ 05_dem_slope.R: wrote DEM + slope (deg) from GLO-30")