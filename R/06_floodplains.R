# ==============================================================================
# 06_floodplains.R — floodplain raster used to mask "demand" pixels
# ------------------------------------------------------------------------------
# Mirrors Duarte's `scripts_OSF/01_harmonize/02_floodplains.R`
#
# Outputs (data/processed/):
#   06_floodplain.tif    — 0/1 raster on the working grid
# ==============================================================================

source(here::here("R", "00_setup.R"))

fp_path <- data_path("mohanty_floodplain_100yr")
ext <- tools::file_ext(fp_path)

if (ext %in% c("tif", "tiff")) {
  fp <- terra::rast(fp_path) |> align_to_grid(method = "near")
  fp <- terra::ifel(fp > 0, 1, 0)
} else {
  v <- sf::st_read(fp_path, quiet = TRUE) |> sf::st_transform(PROJECT_CRS)
  template <- align_to_grid(terra::rast(file.path(paths()$processed, "02_lulc_values.tif")))
  fp <- terra::rasterize(terra::vect(v), template, field = 1, background = 0)
}

outcome_aoi <- read_aoi("outcome")
fp <- terra::mask(fp, terra::rasterize(terra::vect(outcome_aoi), fp, field = 1))

safe_writeRaster(fp, file.path(paths()$processed, "06_floodplain.tif"))

# ---- QA preview --------------------------------------------------------------
qa_png("06_floodplain.png", function() {
  op <- graphics::par(mar = c(2, 2, 3, 6))
  on.exit(graphics::par(op), add = TRUE)
  terra::plot(fp, main = "Floodplain mask (Mohanty CMIP6 100-yr; 1 = inundated)",
              type = "classes", levels = c("dry", "floodplain"),
              col = c("grey90", "#3182bd"))
})

message("✓ 06_floodplains.R — wrote floodplain raster (Mohanty CMIP6 100-yr; swap to FBC when handoff wired)")