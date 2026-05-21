# ==============================================================================
# 03_lulc.R — land use / land cover raster + natural-vegetation mask
# ------------------------------------------------------------------------------
# Steps:
#   1. Load NALCMS 2020, align to the working grid
#   2. Stamp AAFC ACI crop classes on top
#   3. Compute the natural-vegetation mask from NALCMS `is_natural`; refine
#      inside MVRD with Metro Vancouver RNIN patches and corridors
#
# Inputs (resolved via data_path):
#   nalcms_2020, aafc_aci, mvrd_rnin_patch / mvrd_rnin_corridor
#   lookup/nalcms_classes.csv, lookup/aafc_classes.csv
#
# Outputs (data/processed/):
#   03_lulc_values.tif               — single integer raster on the working grid
#   03_lulc_source_coverage.tif      — 1 where NALCMS/AAFC source data was
#                                      projected (used by 99 to distinguish
#                                      "no source data" from "water" in the map)
#   03_lulc_natural_mask.tif         — 1 = natural ecosystem, 0 = not
#   03_lulc_natural_mask_source.tif  — 1 = NALCMS, 2 = RNIN patch, 3 = RNIN corridor
#   03_lulc_class_legend.csv         — definitive lookup for the combined raster
# ==============================================================================

source(here::here("R", "00_setup.R"))

# ---- 1. NALCMS base ----------------------------------------------------------
# snap nalcms dataset to project grid using nearest neighbour resampling
nalcms <- terra::rast(data_path("nalcms_2020")) |> 
  align_to_grid(method = "near")

# ---- 2. Stamp AAFC crop classes on top ---------------------------------------
# look up AAFC class codes
aafc_codes <- readr::read_csv(here::here("lookup", "aafc_classes.csv"),
                              show_col_types = FALSE)

# snap aafc dataset to project grid using nearest neighbour resampling
aafc <- terra::rast(data_path("aafc_aci")) |> 
  align_to_grid(method = "near")

# find aafc codes that are crops only
crop_codes <- aafc_codes$aafc_code[aafc_codes$is_crop]

# build lulc, overriding nalcms with aafc only on crop codes; mask to upstream
lulc <- terra::ifel(!is.na(aafc) & aafc %in% crop_codes, aafc, nalcms)
aoi_mask <- terra::rasterize(terra::vect(read_aoi("upstream")), lulc, field = 1)
lulc <- terra::mask(lulc, aoi_mask)

# track where source data exists, before water/no-data is masked out below
source_coverage <- terra::ifel(!is.na(lulc), 1L, NA_integer_)
safe_writeRaster(source_coverage,
                 file.path(paths()$processed, "03_lulc_source_coverage.tif"))

# warn if NALCMS/AAFC clips don't cover the full upstream AOI; downstream
# sub-basins outside the source footprint will model as NA
src_px <- terra::global(!is.na(source_coverage), "sum", na.rm = TRUE)[1, 1]
aoi_px <- terra::global(!is.na(aoi_mask),         "sum", na.rm = TRUE)[1, 1]
gap    <- 1 - src_px / aoi_px
if (gap > 0.01) {
  message(sprintf(
    "  ! NALCMS/AAFC source covers %.1f%% of the upstream AOI (%d / %d px); ",
    100 * src_px / aoi_px, as.integer(src_px), as.integer(aoi_px)),
    "re-clip raw landcover to 01_aoi_upstream.gpkg's bbox to fill the gap")
}

# drop water (18) and no-data (0); AAFC built-up (34) / barren (30) stay
lulc <- terra::mask(lulc, lulc %in% c(0L, 18L), maskvalue = TRUE)
lulc_src <- terra::ifel(aafc %in% crop_codes, 2L, 1L)  # 1 = NALCMS, 2 = AAFC

safe_writeRaster(lulc, file.path(paths()$processed, "03_lulc_values.tif"))

# ---- 3. Natural-vegetation mask ----------------------------------------------
nalcms_codes <- readr::read_csv(here::here("lookup", "nalcms_classes.csv"),
                                show_col_types = FALSE)

# identify natural NALCMS class codes
nat_codes <- nalcms_codes$nalcms_code[nalcms_codes$is_natural]
# write as binary; 1 = natural, 0 everything else
nat_mask <- terra::ifel(lulc %in% nat_codes, 1, 0)
# keep track of where the source file is for every natural pixel
nat_src  <- terra::ifel(nat_mask == 1, 1, 0)  # 1 = NALCMS

# RNIN corridors (src = 2)
v <- sf::st_read(data_path("mvrd_rnin_corridor"), quiet = TRUE) |> sf::st_transform(PROJECT_CRS)
r <- terra::rasterize(terra::vect(v), nat_mask, field = 1, background = 0)
nat_mask <- terra::ifel(r == 1, 1, nat_mask)
nat_src  <- terra::ifel(r == 1, 2L, nat_src)

# RNIN patches (src = 3)
v <- sf::st_read(data_path("mvrd_rnin_patch"), quiet = TRUE) |> sf::st_transform(PROJECT_CRS)
r <- terra::rasterize(terra::vect(v), nat_mask, field = 1, background = 0)
nat_mask <- terra::ifel(r == 1, 1, nat_mask)
nat_src  <- terra::ifel(r == 1, 3L, nat_src)

safe_writeRaster(nat_mask, file.path(paths()$processed, "03_lulc_natural_mask.tif"))
safe_writeRaster(nat_src,  file.path(paths()$processed, "03_lulc_natural_mask_source.tif"))

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
# File 1: LULC source and natural mask source
qa_png("03_lulc_sources.png", ncol = 2, function() {
  op <- graphics::par(mfrow = c(1, 2), mar = c(3, 3, 2, 4), oma = c(0, 0, 1, 5))
  on.exit(graphics::par(op), add = TRUE)
  
  terra::plot(lulc_src,
              main = "LULC source\n(NALCMS vs AAFC ACI override)",
              type = "classes",
              levels = c("NALCMS", "AAFC ACI"),
              col = c("forestgreen", "goldenrod"),
              ext = terra::ext(lulc))
  
  terra::plot(nat_src,
              main = "Natural mask source\n(NALCMS + MV RNIN refinement)",
              type = "classes",
              levels = c("none", "NALCMS", "RNIN corridor", "RNIN patch"),
              col = c("grey90", "forestgreen", "steelblue", "darkorange"),
              ext = terra::ext(lulc))
})

# File 2: Full LULC classes
qa_png("03_lulc_classes.png", ncol = 1, function() {
  op <- graphics::par(mfrow = c(1, 1), mar = c(2, 3, 2, 2), oma = c(0, 0, 0, 13))
  on.exit(graphics::par(op), add = TRUE)
  
  lulc_plot <- lulc
  lulc_plot[lulc_plot %in% aafc_codes$aafc_code[aafc_codes$is_crop]] <- 200L
  
  ids <- c(0L, nalcms_codes$nalcms_code, 200L)
  cols <- c("lightblue", "darkgreen", "khaki4", "green4", "olivedrab",
            "goldenrod4", "wheat", "rosybrown", "lightgreen", "palegreen3",
            "brown", "yellow3", "grey60", "red", "steelblue",
            "white", "darkorange")
  
  coltab <- data.frame(value = ids, col = cols)
  terra::coltab(lulc_plot) <- coltab
  
  lulc_labeled <- terra::categories(lulc_plot,
                                    value = data.frame(
                                      id    = ids,
                                      label = c("No data", nalcms_codes$nalcms_name, "Cropland (AAFC)")
                                    ))
  
  terra::plot(lulc_labeled,
              main = "LULC classes\n(NALCMS + AAFC ACI)",
              type = "classes",
              ext = terra::ext(lulc))
})

message("✓ 03_lulc.R — wrote LULC values, natural mask, and legend")