# ==============================================================================
# 03_lulc.R: land use / land cover raster + natural-vegetation mask
# ------------------------------------------------------------------------------
# Steps:
#   1. Load NALCMS 2020, align to the working grid
#   2. Stamp AAFC ACI crop classes on top
#   3. Compute the natural-vegetation mask from NALCMS `is_natural`
#
# Inputs (resolved via data_path):
#   nalcms_2020, aafc_aci
#   lookup/nalcms_classes.csv, lookup/aafc_classes.csv
#
# Outputs (data/processed/):
#   03_lulc_values.tif
#   03_lulc_source_coverage.tif
#   03_lulc_natural_mask.tif
#   03_lulc_class_legend.csv
# ==============================================================================

source(here::here("R", "00_setup.R"))

# ---- 1. NALCMS base ----------------------------------------------------------
# this script defines the shared grid, so it builds the template from the AOI
# explicitly; every later script aligns to the raster written here
tmpl <- grid_template()

nalcms <- terra::rast(data_path("nalcms_2020")) |>
  align_to_grid(method = "near", template = tmpl)

# ---- 2. Stamp AAFC crop classes on top ---------------------------------------
aafc_codes <- readr::read_csv(here::here("lookup", "aafc_classes.csv"),
                              show_col_types = FALSE)

aafc <- terra::rast(data_path("aafc_aci")) |>
  align_to_grid(method = "near", template = tmpl)

crop_codes <- aafc_codes$aafc_code[aafc_codes$is_crop]

# override NALCMS with AAFC only on crop codes, then mask to the upstream AOI
lulc <- terra::ifel(!is.na(aafc) & aafc %in% crop_codes, aafc, nalcms)
aoi_mask <- terra::rasterize(terra::vect(read_aoi("upstream")), lulc, field = 1)
lulc <- terra::mask(lulc, aoi_mask)

# track where source data exists, before water/no-data is masked out below
source_coverage <- terra::ifel(!is.na(lulc), 1L, NA_integer_)
safe_writeRaster(source_coverage,
                 file.path(paths()$processed, "03_lulc_source_coverage.tif"))

# warn if the NALCMS/AAFC clips is smaller than the upstream AOI clip
src_px <- terra::global(!is.na(source_coverage), "sum", na.rm = TRUE)[1, 1]
aoi_px <- terra::global(!is.na(aoi_mask), "sum", na.rm = TRUE)[1, 1]
gap <- 1 - src_px / aoi_px
if (is.finite(gap) && gap > 0.001) {
  warning(sprintf(paste0(
    "LULC sources cover only %.1f%% of the upstream AOI (%.1f%% gap); ",
    "sub-basins outside the NALCMS/AAFC footprint will model as NA. ",
    "check the raw clips against the expanded AOI"),
    100 * (1 - gap), 100 * gap), immediate. = TRUE)
}

# drop water (18) and no-data (0); # 1 = NALCMS, 2 = AAFC
lulc <- terra::mask(lulc, lulc %in% c(0L, 18L), maskvalue = TRUE)
lulc_src <- terra::mask(terra::ifel(aafc %in% crop_codes, 2L, 1L), lulc)

safe_writeRaster(lulc, file.path(paths()$processed, "03_lulc_values.tif"))

# ---- 3. Natural-vegetation mask ----------------------------------------------
# 1 = natural NALCMS class, 0 = otherwise
nalcms_codes <- readr::read_csv(here::here("lookup", "nalcms_classes.csv"),
                                show_col_types = FALSE)

nat_codes <- nalcms_codes$nalcms_code[nalcms_codes$is_natural]
nat_mask <- terra::ifel(lulc %in% nat_codes, 1, 0)

safe_writeRaster(nat_mask, file.path(paths()$processed, "03_lulc_natural_mask.tif"))

# ---- 4. Legend CSV -----------------------------------------------------------
legend <- nalcms_codes |>
  dplyr::transmute(class_id = nalcms_code, class_name = nalcms_name,
                   source = "nalcms", is_natural)

aafc_legend <- aafc_codes |>
  dplyr::filter(is_crop) |>
  dplyr::transmute(class_id = aafc_code, class_name = aafc_name,
                   source = "aafc", is_natural = FALSE)

legend <- dplyr::bind_rows(legend, aafc_legend)
readr::write_csv(legend, file.path(paths()$processed, "03_lulc_class_legend.csv"))

# ---- 5. Visualizations -------------------------------------------------------
# NALCMS water (class 18) for river/lake overlay on the classes map
water_qa <- terra::ifel(!is.na(aoi_mask) & nalcms == 18L, 1L, NA_integer_)

# does the AAFC crop override land where expected, and does the natural mask
# pick out the vegetated classes
qa_png("03_lulc_sources.png", ncol = 2, function() {
  op <- graphics::par(mfrow = c(1, 2))
  on.exit(graphics::par(op), add = TRUE)
  terra::plot(lulc_src, type = "classes", levels = c("NALCMS", "AAFC ACI"),
              col = c("forestgreen", "goldenrod"), axes = TRUE, main = "")
  terra::plot(nat_mask, type = "classes", levels = c("not natural", "natural"),
              col = c("grey85", "forestgreen"), axes = TRUE, main = "")
})

# every class present in the model, with mapped water drawn over it
qa_png("03_lulc_classes.png", panel_w = 1400, function() {
  lulc_plot <- lulc
  lulc_plot[lulc_plot %in% aafc_codes$aafc_code[aafc_codes$is_crop]] <- 200L

  ids <- c(0L, nalcms_codes$nalcms_code, 200L)
  cols <- c("lightblue", "darkgreen", "khaki4", "green4", "olivedrab",
            "goldenrod4", "wheat", "rosybrown", "lightgreen", "palegreen3",
            "brown", "yellow3", "grey60", "red", "steelblue",
            "white", "darkorange")
  terra::coltab(lulc_plot) <- data.frame(value = ids, col = cols)
  lulc_plot <- terra::categories(lulc_plot, value = data.frame(
    id = ids,
    label = c("No data", nalcms_codes$nalcms_name, "Cropland (AAFC)")))

  terra::plot(lulc_plot, type = "classes", axes = TRUE, main = "")
  terra::plot(water_qa, add = TRUE, legend = FALSE, col = "#74add1")
})

message("✓ 03_lulc.R: wrote LULC values, natural mask, and legend")