# ==============================================================================
# 04_dem_slope.R — DEM and slope raster
# ------------------------------------------------------------------------------
# Mirrors Duarte's `scripts_OSF/01_harmonize/02_slope.R` but with Copernicus
# GLO-30 in place of HydroSHEDS:
#
#   DEM:   Copernicus GLO-30 (30 m native DSM, TanDEM-X SAR)
#   fill:  Duarte 10° fallback over LULC pixels where the DEM is NA
#
# HRDEM (NRCan LiDAR, 1–2 m) is parked in `R/extensions/prep_hrdem_mosaic.R`
# as a Phase-2 upgrade
#
# Outputs (data/processed/):
#   04_dem.tif          — DEM on the working grid (GLO-30)
#   04_slope_deg.tif    — slope in *degrees*
#
# Notes:
#   - Duarte's slope unit (degrees) is preserved because 08_curve_numbers.R
#     applies the Huang TR-55 correction with the convention `α = degrees/100`.
# ==============================================================================

source(here::here("R", "00_setup.R"))

dem_path <- data_path("dem_glo30", must_exist = TRUE)
message("  · loading DEM 'dem_glo30' from ", dem_path)
dem <- terra::rast(dem_path) |> align_to_grid(method = "bilinear")

slope_real <- terra::terrain(dem, v = "slope", neighbors = 8, unit = "degrees")
lulc <- terra::rast(file.path(paths()$processed, "02_lulc_values.tif"))
slope_deg <- terra::ifel(is.na(slope_real) & !is.na(lulc), 10, slope_real)
aoi_mask <- terra::rasterize(terra::vect(read_aoi("upstream")), slope_deg, field = 1)
slope_deg <- terra::mask(slope_deg, aoi_mask)

safe_writeRaster(dem, file.path(paths()$processed, "04_dem.tif"))
safe_writeRaster(slope_deg, file.path(paths()$processed, "04_slope_deg.tif"))

# coverage diagnostic
land_mask <- !is.na(lulc)
total_px  <- terra::global(land_mask, "sum", na.rm = TRUE)[1, 1]
fill_mask <- is.na(slope_real) & land_mask
fill_px   <- terra::global(fill_mask, "sum", na.rm = TRUE)[1, 1]
dem_px    <- total_px - fill_px

message(sprintf("  · GLO-30-derived slope:               %10d px  (%5.1f%%)",
                as.integer(dem_px), 100 * dem_px / total_px))

message(sprintf("  · 10° Duarte fallback (Huang ≈ off):  %10d px  (%5.1f%%)",
                as.integer(fill_px), 100 * fill_px / total_px))

if (fill_px / total_px > 0.05) {
  message("  ! >5% of slope pixels are still on the 10° fallback. The ",
          "GLO-30 clip probably doesn't fully cover the upstream AOI — ",
          "re-export it from Earth Engine for the full provider footprint.")
}

# ---- QA preview --------------------------------------------------------------
qa_png("04_dem_slope.png", function() {
  op <- graphics::par(mar = c(2, 2, 3, 7))
  on.exit(graphics::par(op), add = TRUE)
  terra::plot(slope_deg, main = "Slope (degrees)",
              col = grDevices::hcl.colors(50, palette = "YlOrRd", rev = FALSE),
              range = c(0, 60))
})

message("✓ 04_dem_slope.R — wrote DEM + slope (deg) from GLO-30")