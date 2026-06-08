# ==============================================================================
# 10_demand.R — Service-benefiting "demand" per sub-basin
# ------------------------------------------------------------------------------
# Mirrors Duarte's `scripts_OSF/03_sba/01_beneficiaries_raster.R` and
# `02_beneficiaries_basins.R`. The "demand" layer at each pixel is
# (cropland + built-up) area, masked to the floodplain. Per sub-basin we sum
# that area (km²) and store three flavours:
#
#   ben_both        —  crops + built-up
#   ben_built       —          built-up
#   ben_crops       —  crops
#   ben_crops_w     —  crops weighted by lookup/crop_flood_vulnerability.csv
#   total_demand_w  =  ben_built + ben_crops_w
#
# Note: livestock demand (Census of Agriculture CSD-level head counts) is a
# planned extension once the agri_census_2021 dataset is available.
#
# Demand outside the downstream AOI (MVRD ∪ FVRD — the Hope-to-coast Lower
# Mainland strip) is set to zero so 11_routing_decay.R only credits upstream
# sub-basins for protecting value at risk inside the Lower Mainland.
#
# Inputs (data/processed/):
#   03_lulc_values.tif, 03_lulc_class_legend.csv
#   07_floodplain.tif
#   02_subbasins.gpkg
#   01_aoi_downstream.gpkg
#
# Inputs (lookup/):
#   crop_flood_vulnerability.csv
#
# Outputs (data/processed/):
#   10_demand_pixel.tif        beneficiary indicator × pixel area km²
#   10_demand_subbasin.gpkg    one row per HYBAS_ID with the columns above
# ==============================================================================

source(here::here("R", "00_setup.R"))

lulc <- terra::rast(file.path(paths()$processed, "03_lulc_values.tif"))
fp   <- terra::rast(file.path(paths()$processed, "07_floodplain.tif"))
sb   <- sf::st_read(file.path(paths()$processed, "02_subbasins.gpkg"), quiet = TRUE)
downstream_aoi <- read_aoi("downstream")

# ---- 1. Pixel masks ----------------------------------------------------------
# NALCMS class 17 = Urban and built-up; class 15 = Cropland (replaced by AAFC
# crop codes for any pixel inside the AAFC ACI footprint, in 03_lulc.R).
aafc_codes <- readr::read_csv(here::here("lookup", "aafc_classes.csv"),
                              show_col_types = FALSE)
crop_codes <- c(15L, aafc_codes$aafc_code[aafc_codes$is_crop])

built_mask <- terra::ifel(lulc == 17, 1, 0)
crop_mask  <- terra::ifel(lulc %in% crop_codes, 1, 0)

# Mask to floodplain ∩ downstream AOI.
downstream_r  <- terra::rasterize(terra::vect(downstream_aoi), built_mask, field = 1, background = 0)
in_downstream <- (fp == 1) & (downstream_r == 1)
built_mask <- terra::ifel(in_downstream, built_mask, 0)
crop_mask  <- terra::ifel(in_downstream, crop_mask,  0)

# Pixel area km² (constant at 30 m: 9 × 10⁻⁴).
px_km2 <- (terra::xres(lulc) * terra::yres(lulc)) * 1e-6

ben_b  <- (built_mask + crop_mask) * px_km2
ben_bu <- built_mask               * px_km2
ben_cr <- crop_mask                * px_km2

safe_writeRaster(ben_b, file.path(paths()$processed, "10_demand_pixel.tif"))

# ---- 2. Aggregate per sub-basin ----------------------------------------------
sb$ben_both  <- exactextractr::exact_extract(ben_b,  sb, "sum")
sb$ben_built <- exactextractr::exact_extract(ben_bu, sb, "sum")
sb$ben_crops <- exactextractr::exact_extract(ben_cr, sb, "sum")

# ---- 3. Crop-vulnerability-weighted demand -----------------------------------
vul <- readr::read_csv(here::here("lookup", "crop_flood_vulnerability.csv"),
                       show_col_types = FALSE)
rcl <- as.matrix(vul[, c("aafc_code", "vulnerability_score")])
vul_r <- terra::classify(lulc, rcl, others = 0)
ben_cr_w <- terra::ifel(in_downstream, vul_r * px_km2, 0)
sb$ben_crops_w <- exactextractr::exact_extract(ben_cr_w, sb, "sum")

sb$total_demand_w <- sb$ben_built + sb$ben_crops_w

sf::st_write(sb, file.path(paths()$processed, "10_demand_subbasin.gpkg"),
             delete_dsn = TRUE, quiet = TRUE)

# ---- QA preview --------------------------------------------------------------
# (a) the pixel-level beneficiary mask (built + crops, floodplain ∩ downstream)
# (b) the per-sub-basin total weighted demand that 11_routing_decay.R will
#     propagate upstream.
qa_png("10_demand.png", ncol = 2, function() {
  op <- graphics::par(mfrow = c(1, 2), mar = c(2, 2, 3, 7))
  on.exit(graphics::par(op), add = TRUE)
  terra::plot(ben_b, main = "Built + crop area in floodplain (km²/pixel)",
              col = grDevices::hcl.colors(50, palette = "Reds 3", rev = TRUE))
  plot(sb["total_demand_w"],
       main = "Per-sub-basin total weighted demand",
       border = "grey60", lwd = 0.2, key.pos = 4, reset = FALSE)
  plot(sf::st_geometry(downstream_aoi), add = TRUE, border = "black", lwd = 1.2)
})

message("✓ 10_demand.R — wrote pixel + per-sub-basin demand layers")