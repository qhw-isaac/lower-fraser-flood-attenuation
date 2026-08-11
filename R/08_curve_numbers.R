# ==============================================================================
# 08_curve_numbers.R: Curve Number raster + Huang slope correction
# ------------------------------------------------------------------------------
# Algorithm (per pixel):
#
#   CN_baseline = CN[lulc, hsg]
#   CN_counterfactual = CN[replace_natural_with_barren(lulc), hsg]
#   a = slope_deg / 100 # Huang et al.
#   CN_slopeadj = if a ≥ 0.14: CN × (322.79 + 15.63*a) / (a + 323.52)
#                 CN otherwise
#
# Inputs (data/processed/):
#   03_lulc_values.tif, 03_lulc_natural_mask.tif
#   04_hsg.tif, 05_slope_deg.tif
#   lookup/cn_lookup.csv
#
# Outputs (data/processed/):
#   08_cn_baseline.tif
#   08_cn_counterfactual.tif
#   08_cn_baseline_slopeadj.tif
#   08_cn_counterfactual_slopeadj.tif
# ==============================================================================

source(here::here("R", "00_setup.R"))

lulc <- terra::rast(file.path(paths()$processed, "03_lulc_values.tif"))
hsg <- terra::rast(file.path(paths()$processed, "04_hsg.tif"))
slope_deg <- terra::rast(file.path(paths()$processed, "05_slope_deg.tif"))

# one row per (land-cover class, soil group) pair
cn_tbl <- readr::read_csv(here::here("lookup", "cn_lookup.csv"),
                          show_col_types = FALSE) |>
  tidyr::pivot_longer(c(CN_A, CN_B, CN_C, CN_D),
                      names_to = "hsg_letter", values_to = "cn") |>
  dplyr::mutate(hsg = match(sub("CN_", "", hsg_letter), c("A", "B", "C", "D"))) |>
  dplyr::select(lulc_code, hsg, cn) |>
  dplyr::filter(!is.na(cn))

# stop if duplicated (lulc_code, hsg) pair are identified
stopifnot(!anyDuplicated(cn_tbl[, c("lulc_code", "hsg")]))

# barren counterfactual: replace every natural pixel with code 16 (Barren)
natural_mask <- terra::rast(file.path(paths()$processed, "03_lulc_natural_mask.tif"))

BARREN_CODE <- 16L

# ---- CN rasters --------------------------------------------------------------
# look up CN by (LULC class, HSG)
cn_from_lulc <- function(lulc_r) {
  key <- lulc_r * 10L + hsg
  terra::classify(key,
                  cbind(cn_tbl$lulc_code * 10 + cn_tbl$hsg, cn_tbl$cn),
                  others = NA)
}

cn_baseline <- cn_from_lulc(lulc)

# counterfactual CN, with natural pixels recoded to barren before lookup
cn_cf <- cn_from_lulc(terra::ifel(natural_mask == 1, BARREN_CODE, lulc))

safe_writeRaster(cn_baseline, file.path(paths()$processed, "08_cn_baseline.tif"))
safe_writeRaster(cn_cf, file.path(paths()$processed, "08_cn_counterfactual.tif"))

# ---- Huang slope correction --------------------------------------------------
alpha <- slope_deg / 100 # Huang convention: degrees / 100 ≈ m/m
slope_factor <- (322.79 + 15.63 * alpha) / (alpha + 323.52)

cn_baseline_adj <- terra::ifel(alpha >= 0.14, cn_baseline * slope_factor, cn_baseline)
cn_baseline_adj <- terra::ifel(cn_baseline_adj > 100, 100, cn_baseline_adj)

cn_cf_adj <- terra::ifel(alpha >= 0.14, cn_cf * slope_factor, cn_cf)
cn_cf_adj <- terra::ifel(cn_cf_adj > 100, 100, cn_cf_adj)

safe_writeRaster(cn_baseline_adj, file.path(paths()$processed, "08_cn_baseline_slopeadj.tif"))
safe_writeRaster(cn_cf_adj, file.path(paths()$processed, "08_cn_counterfactual_slopeadj.tif"))

# ---- QA preview --------------------------------------------------------------
# baseline and bare-ground curve numbers, and the uplift between them, which is
# the quantity every retention estimate downstream is derived from
qa_png("08_cn_slopeadj.png", ncol = 3, function() {
  op <- graphics::par(mfrow = c(1, 3))
  on.exit(graphics::par(op), add = TRUE)
  pal <- grDevices::hcl.colors(50, "YlOrBr")
  terra::plot(cn_baseline_adj, col = pal, range = c(20, 100), axes = TRUE, main = "")
  terra::plot(cn_cf_adj, col = pal, range = c(20, 100), axes = TRUE, main = "")
  terra::plot(sqrt(terra::clamp(cn_cf_adj - cn_baseline_adj, lower = 0)),
              col = grDevices::hcl.colors(50, "Blues", rev = TRUE),
              axes = TRUE, main = "")
})

message("✓ 08_curve_numbers.R: wrote 4 CN rasters ",
        "(baseline/counterfactual, ±slope; AMC II only)")