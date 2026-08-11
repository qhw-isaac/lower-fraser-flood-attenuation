# ==============================================================================
# 04_soils_hysog.R: Hydrologic Soil Group (HSG) raster
# ------------------------------------------------------------------------------
# HYSOGs250m provides the soil type used in the Curve Number calculation. 
# Outputs read by 08_curve_numbers.R:
#
#   04_hsg.tif          integer HSG, 1=A, 2=B, 3=C, 4=D
#   04_hsg_source.tif
#
# HYSOGs250m codes 1-4 as HSG A-D and 11-14 as the dual groups A/D, B/D, C/D,
# D/D. Dual groups are collapsed to D following Duarte et al. (2024).
# ==============================================================================

source(here::here("R", "00_setup.R"))

# ---- Inputs ------------------------------------------------------------------

lulc <- terra::rast(file.path(paths()$processed, "03_lulc_values.tif"))

hy <- terra::rast(data_path("hysogs250m"))

# align_to_grid() projects onto 03's lulc raster itself, so geometry matches
hy <- hy |>
  align_to_grid(method = "near")

# dual groups (11-14) -> D, keep only A-D, drop water/no-data
hy <- terra::ifel(hy %in% c(11, 12, 13, 14), 4L, hy)
hy <- terra::ifel(hy %in% 1:4, hy, NA)

# D-fill (D = the least infiltrative group / highest runoff) where land cover 
# exists but HYSOG has no soil to avoid silent drops during CN lookup
hsg <- terra::ifel(is.na(hy), 4L, hy)
hsg <- terra::mask(hsg, lulc)

# source: 3 = real HYSOG soil, 4 = D-fill
hsg_src <- terra::ifel(
  !is.na(hy),
  3L,
  terra::ifel(!is.na(hsg), 4L, NA)
)

hsg_src <- terra::mask(hsg_src, lulc)

# ---- Outputs -----------------------------------------------------------------
safe_writeRaster(hsg, file.path(paths()$processed, "04_hsg.tif"))
safe_writeRaster(hsg_src, file.path(paths()$processed, "04_hsg_source.tif"))

# ---- Source diagnostics ------------------------------------------------------
src_freq <- terra::freq(hsg_src)

src_lab <- c(
  `3` = "HYSOGs250m",
  `4` = "D-fill (no soil data)"
)

total_px <- sum(src_freq$count, na.rm = TRUE)

message("HSG source coverage:")

for (i in seq_len(nrow(src_freq))) {
  k <- as.character(src_freq$value[i])
  
  message(sprintf(
    "  · %-24s %10d px  (%5.1f%%)",
    src_lab[k],
    src_freq$count[i],
    100 * src_freq$count[i] / total_px
  ))
}

# ---- QA preview --------------------------------------------------------------
# the soil groups across the study area and how much of the map is D-fill
qa_png("04_hsg_hysog.png", ncol = 2, function() {
  op <- graphics::par(mfrow = c(1, 2))
  on.exit(graphics::par(op), add = TRUE)
  terra::plot(hsg, type = "classes",
              levels = c("A (sand)", "B (loam)", "C (clay loam)", "D (clay)"),
              col = c("#fee08b", "#fdae61", "#d73027", "#7f0000"),
              axes = TRUE, main = "")
  terra::plot(hsg_src, type = "classes", levels = c("HYSOGs250m", "D-fill"),
              col = c("#d95f02", "#7570b3"), axes = TRUE, main = "")
})

message("✓ 04_soils_hysog.R: wrote HSG raster using HYSOGs250m.")