# ==============================================================================
# 05_precipitation.R — Wettest-month precipitation raster (P, in mm)
# ------------------------------------------------------------------------------
# PCIC PRISM 800 m monthly normals — pixel-wise max across 12 months, matching
# the CHELSA bio13 (precipitation of wettest month) convention used in
# Duarte et al. (2024).
#
# This is the sole precipitation input to the CN model. The wettest-month
# normal is used rather than an event-based scenario because the goal is to
# characterize the *relative spatial pattern* of flood attenuation capacity
# across the landscape — which ecosystem areas retain the most runoff under
# climatologically high precipitation. Absolute PRR magnitudes are
# scenario-dependent but the spatial ranking is stable across plausible
# precipitation inputs (Duarte et al. 2024).
#
# The November 2021 atmospheric river event motivates the policy relevance of
# this analysis but is not modelled explicitly; the wettest-month normal
# captures the seasonal precipitation regime that makes the Lower Fraser
# corridor flood-prone.
#
# Inputs:
#   data_path("prism_pcic")   PCIC PRISM 800 m monthly normals (12 bands)
#
# Outputs (data/processed/precip/):
#   05_p_wettest_month.tif    wettest-month P (mm), aligned to working grid
# ==============================================================================
source(here::here("R", "00_setup.R"))

precip_dir <- file.path(paths()$processed, "precip")

dir.create(precip_dir, showWarnings = FALSE, recursive = TRUE)

prism <- terra::rast(data_path("prism_pcic")) |> align_to_grid(method = "bilinear")

p <- terra::app(prism, fun = max, na.rm = TRUE)

lulc <- terra::rast(file.path(paths()$processed, "02_lulc_values.tif"))

p <- terra::mask(p, lulc)

safe_writeRaster(p, file.path(precip_dir, "05_p_wettest_month.tif"))

# ---- QA preview --------------------------------------------------------------
qa_png("05_p_wettest_month.png", function() {
  op <- graphics::par(mar = c(2, 2, 3, 6))
  on.exit(graphics::par(op), add = TRUE)
  terra::plot(p, main = "Wettest-month precipitation (mm; PRISM 800 m, pixel-wise max)",
              col = grDevices::hcl.colors(50, palette = "Blues 3", rev = TRUE))
})

message("✓ 05_precipitation.R — wrote wettest-month P raster (mm)")