# ==============================================================================
# 07_floodplains.R — floodplain raster used to mask "demand" pixels
# ------------------------------------------------------------------------------
# Creates the floodplain used in subsequent analysis.
#
# Inputs (data/processed/):
#   03_lulc_values.tif
#
# Outputs (data/processed/):
#   07_floodplain.tif
# ==============================================================================

source(here::here("R", "00_setup.R"))

lulc_path <- file.path(paths()$processed, "03_lulc_values.tif")
lulc_grid <- terra::rast(lulc_path)

fp_path <- data_path("mohanty_floodplain_100yr")
ext <- tools::file_ext(fp_path)

if (ext %in% c("tif", "tiff")) {
  fp <- terra::rast(fp_path) |> terra::project(lulc_grid, method = "near")
  fp <- terra::ifel(fp > 0, 1, 0)
} else {
  v <- sf::st_read(fp_path, quiet = TRUE) |> sf::st_transform(PROJECT_CRS)
  fp <- terra::rasterize(terra::vect(v), lulc_grid, field = 1, background = 0)
}

# keep only floodplain inside the downstream AOI (where demand is evaluated)
downstream_aoi <- read_aoi("downstream")
fp <- terra::mask(fp, terra::rasterize(terra::vect(downstream_aoi), fp, field = 1))

safe_writeRaster(fp, file.path(paths()$processed, "07_floodplain.tif"))

# ---- QA preview --------------------------------------------------------------
qa_png("07_floodplain.png", function() {
  op <- graphics::par(mar = c(2, 2, 3, 6))
  on.exit(graphics::par(op), add = TRUE)
  terra::plot(terra::trim(fp),
              main = "Floodplain mask (Mohanty CMIP6 100-yr; 1 = inundated)",
              type = "classes", levels = c("dry", "floodplain"),
              col = c("grey90", "#3182bd"))
})

message("✓ 07_floodplains.R — wrote floodplain raster (Mohanty CMIP6 100-yr; swap to FBC when handoff wired)")