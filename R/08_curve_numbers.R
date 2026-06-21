# ==============================================================================
# 08_curve_numbers.R — Curve Number raster + Huang slope correction
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
# Antecedent moisture: published CN values are AMC II (average wetness). The
# slope-adjusted CN(II) is also converted to AMC I (dry) and AMC III (wet), so a
# storm can be run on dry, average, or saturated ground. AMC III is the realistic
# state for the Nov 2021 atmospheric river, which fell on soils already saturated
# by an earlier storm.
#
# Outputs (data/processed/):
#   08_cn_baseline.tif
#   08_cn_counterfactual.tif
#   08_cn_baseline_slopeadj.tif                  (AMC II, default)
#   08_cn_counterfactual_slopeadj.tif            (AMC II, default)
#   08_cn_baseline_slopeadj_amc{I,III}.tif
#   08_cn_counterfactual_slopeadj_amc{I,III}.tif
# ==============================================================================

source(here::here("R", "00_setup.R"))

# load in previous script outputs
lulc <- terra::rast(file.path(paths()$processed, "03_lulc_values.tif"))
hsg <- terra::rast(file.path(paths()$processed, "04_hsg.tif"))
slope_deg <- terra::rast(file.path(paths()$processed, "05_slope_deg.tif"))

# Curve Number lookup table formatting
cn_tbl <- readr::read_csv(here::here("lookup", "cn_lookup.csv"),
                          show_col_types = FALSE) |>
  tidyr::pivot_longer(c(CN_A, CN_B, CN_C, CN_D),
                      names_to = "hsg_letter", values_to = "cn") |>
  dplyr::mutate(hsg = match(sub("CN_", "", hsg_letter), c("A", "B", "C", "D"))) |>
  dplyr::select(lulc_code, hsg, cn) |>
  dplyr::filter(!is.na(cn))

# barren counterfactual: replace every natural NALCMS class with code 16 (Barren)
nalcms_codes <- readr::read_csv(here::here("lookup", "nalcms_classes.csv"),
                                show_col_types = FALSE)

natural_codes <- nalcms_codes$nalcms_code[nalcms_codes$is_natural]

BARREN_CODE <- 16L

# ---- CN rasters --------------------------------------------------------------
# baseline CN
df_base <- data.frame(lulc_code = as.vector(lulc), hsg = as.vector(hsg)) |>
  dplyr::left_join(cn_tbl, by = c("lulc_code", "hsg"))

cn_baseline <- lulc
terra::values(cn_baseline) <- df_base$cn

# counterfactual CN — natural pixels recoded to barren before lookup
lulc_cf <- terra::ifel(lulc %in% natural_codes, BARREN_CODE, lulc)
df_cf <- data.frame(lulc_code = as.vector(lulc_cf), hsg = as.vector(hsg)) |>
  dplyr::left_join(cn_tbl, by = c("lulc_code", "hsg"))

cn_cf <- lulc_cf
terra::values(cn_cf) <- df_cf$cn

safe_writeRaster(cn_baseline, file.path(paths()$processed, "08_cn_baseline.tif"))
safe_writeRaster(cn_cf, file.path(paths()$processed, "08_cn_counterfactual.tif"))

# ---- Huang slope correction --------------------------------------------------
alpha <- slope_deg / 100 # Duarte/Huang convention: degrees / 100 ≈ m/m
slope_factor <- (322.79 + 15.63 * alpha) / (alpha + 323.52)

cn_baseline_adj <- terra::ifel(alpha >= 0.14, cn_baseline * slope_factor, cn_baseline)
cn_baseline_adj <- terra::ifel(cn_baseline_adj > 100, 100, cn_baseline_adj)

cn_cf_adj <- terra::ifel(alpha >= 0.14, cn_cf * slope_factor, cn_cf)
cn_cf_adj <- terra::ifel(cn_cf_adj > 100, 100, cn_cf_adj)

safe_writeRaster(cn_baseline_adj, file.path(paths()$processed, "08_cn_baseline_slopeadj.tif"))
safe_writeRaster(cn_cf_adj, file.path(paths()$processed, "08_cn_counterfactual_slopeadj.tif"))

# ---- Antecedent moisture condition (AMC) variants ----------------------------
# Convert the slope-adjusted CN(II) to dry (AMC I) and wet (AMC III) conditions.
# Standard USDA NRCS conversions (Chow, Maidment & Mays 1988; Mishra & Singh
# 2003): AMC III raises CN (less storage, more runoff), AMC I lowers it.
cn_amcI   <- function(cn) terra::clamp(4.2 * cn / (10 - 0.058 * cn), 1, 100)
cn_amcIII <- function(cn) terra::clamp(23  * cn / (10 + 0.13  * cn), 1, 100)

amc_fns <- list(amcI = cn_amcI, amcIII = cn_amcIII)
for (nm in names(amc_fns)) {
  f <- amc_fns[[nm]]
  safe_writeRaster(f(cn_baseline_adj),
    file.path(paths()$processed, paste0("08_cn_baseline_slopeadj_", nm, ".tif")))
  safe_writeRaster(f(cn_cf_adj),
    file.path(paths()$processed, paste0("08_cn_counterfactual_slopeadj_", nm, ".tif")))
}

# ---- QA preview --------------------------------------------------------------
# 3 panel baseline vs barren counterfactual vs difference
qa_png("08_cn_slopeadj.png", ncol = 3, function() {
  op <- graphics::par(mfrow = c(1, 3), mar = c(2, 2, 3, 6))
  on.exit(graphics::par(op), add = TRUE)
  pal <- grDevices::hcl.colors(50, palette = "YlOrBr", rev = FALSE)
  pal_diff <- grDevices::hcl.colors(50, palette = "Blues", rev = TRUE)
  
  terra::plot(cn_baseline_adj, main = "CN baseline (slope-adj)",
              col = pal, range = c(20, 100))
  terra::plot(cn_cf_adj, main = "CN counterfactual (slope-adj)",
              col = pal, range = c(20, 100))
  terra::plot(sqrt(cn_cf_adj - cn_baseline_adj), main = "CN difference (CF - baseline, sqrt scale)",
              col = pal_diff)
})

# Baseline CN under the three antecedent-moisture conditions (dry → wet).
qa_png("08_cn_amc.png", ncol = 3, function() {
  op <- graphics::par(mfrow = c(1, 3), mar = c(2, 2, 3, 6))
  on.exit(graphics::par(op), add = TRUE)
  pal <- grDevices::hcl.colors(50, palette = "YlOrBr", rev = FALSE)
  terra::plot(cn_amcI(cn_baseline_adj),   main = "CN baseline — AMC I (dry)",
              col = pal, range = c(20, 100))
  terra::plot(cn_baseline_adj,            main = "CN baseline — AMC II (average)",
              col = pal, range = c(20, 100))
  terra::plot(cn_amcIII(cn_baseline_adj), main = "CN baseline — AMC III (wet)",
              col = pal, range = c(20, 100))
})

message("✓ 08_curve_numbers.R — wrote 8 CN rasters (baseline/counterfactual, ",
        "±slope, AMC I/II/III)")