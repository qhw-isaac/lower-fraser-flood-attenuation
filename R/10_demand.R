# ==============================================================================
# 10_demand.R — Service-benefiting "demand" per sub-basin
# ------------------------------------------------------------------------------
# Per-pixel "demand" = (cropland + built-up) area masked to the floodplain, 
# summed per sub-basin (km²) in the following variations:
#
#   ben_both — crops + built-up
#   ben_built — built-up
#   ben_crops — crops
#   ben_crops_w — crops weighted by lookup/crop_flood_vulnerability.csv
#   total_demand_w = ben_built + ben_crops_w
#
# Demand outside the downstream AOI (MVRD ∪ FVRD) is zeroed, so 11_routing_decay.R
# only credits upstream sub-basins for value at risk inside the Lower Mainland.
#
# Inputs: 03_lulc_values.tif, 07_floodplain.tif, 02_subbasins.gpkg,
#          01_aoi_downstream.gpkg, lookup/crop_flood_vulnerability.csv
#
# Outputs (data/processed/):
#   10_demand_pixel.tif
#   10_demand_subbasin.gpkg
# ==============================================================================

source(here::here("R", "00_setup.R"))

lulc <- terra::rast(file.path(paths()$processed, "03_lulc_values.tif"))
fp <- terra::rast(file.path(paths()$processed, "07_floodplain.tif"))
sb <- sf::st_read(file.path(paths()$processed, "02_subbasins.gpkg"), quiet = TRUE)
downstream_aoi <- read_aoi("downstream")

# ---- 1. Pixel masks ----------------------------------------------------------
# NALCMS:
# class 17 = Urban and built-up
# class 15 = Cropland (replaced by AAFC crop codes where available)
aafc_codes <- readr::read_csv(here::here("lookup", "aafc_classes.csv"),
                              show_col_types = FALSE)
crop_codes <- c(15L, aafc_codes$aafc_code[aafc_codes$is_crop])

# binary per-pixel masks
built_mask <- terra::ifel(lulc == 17, 1, 0)
crop_mask <- terra::ifel(lulc %in% crop_codes, 1, 0)

# mask to floodplain ∩ downstream AOI
downstream_r <- terra::rasterize(terra::vect(downstream_aoi), built_mask, field = 1, background = 0)
in_downstream <- (fp == 1) & (downstream_r == 1)
built_mask <- terra::ifel(in_downstream, built_mask, 0)
crop_mask <- terra::ifel(in_downstream, crop_mask, 0)

# pixel area km^2
px_km2 <- (terra::xres(lulc) * terra::yres(lulc)) * 1e-6

ben_b <- (built_mask + crop_mask) * px_km2 # beneficiary area km²/pixel (tif + QA)
safe_writeRaster(ben_b, file.path(paths()$processed, "10_demand_pixel.tif"))

# ---- 2. Aggregate per sub-basin ----------------------------------------------
# sum the 0/1 masks per basin, then scale by the constant pixel area
sb$ben_built <- px_km2 * exactextractr::exact_extract(built_mask, sb, "sum")
sb$ben_crops <- px_km2 * exactextractr::exact_extract(crop_mask, sb, "sum")
sb$ben_both <- sb$ben_built + sb$ben_crops

# ---- 3. Crop-vulnerability-weighted demand -----------------------------------
# Extension beyond Duarte (who treat every km^2 of beneficiary land equally)
# Each AAFC crop class carries a flood-vulnerability score (lookup/), so a km^2
# of a high-value/high-loss crop counts for more demand than a km^2 of pasture
vul <- readr::read_csv(here::here("lookup", "crop_flood_vulnerability.csv"),
                       show_col_types = FALSE)
# classify(): remap each AAFC code to its vulnerability score; all else assigned 0,
# giving a per-pixel weight in place of the 0/1 crop mask above
rcl <- as.matrix(vul[, c("aafc_code", "vulnerability_score")])
vul_r <- terra::ifel(in_downstream, terra::classify(lulc, rcl, others = 0), 0)
sb$ben_crops_w <- px_km2 * exactextractr::exact_extract(vul_r, sb, "sum")

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