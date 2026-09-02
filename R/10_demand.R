# ==============================================================================
# 10_demand.R: Service-benefiting "demand" per sub-basin
# ------------------------------------------------------------------------------
# Floodplain ∩ downstream AOI, cropland + built-up, summed per sub-basin (km²):
#
#   ben_both         crops + built-up
#   ben_built        built-up
#   ben_crops        crops
#   ben_crops_w      crops weighted by lookup/crop_flood_vulnerability.csv
#   ben_facilities   schools/hospitals (km²-equivalent; not in total_demand_w)
#   n_facilities     facility count per sub-basin
#   total_demand_w   ben_built + ben_crops_w (the routed column)
#
# Inputs (data/processed/):
#   03_lulc_values.tif, 07_floodplain.tif, 02_subbasins.gpkg
#   01_aoi_downstream.gpkg, lookup/crop_flood_vulnerability.csv,
#   lookup/critical_facilities.csv,
#   raw/downstream/{bc_k12_schools_2026_databc.csv, hlbc_hospitals.csv}
#
# Outputs (data/processed/):
#   10_demand_pixel.tif
#   10_demand_subbasin.gpkg
#   10_demand_da.gpkg   same demand by 2021 dissemination area (exposed DAs)
# ==============================================================================

source(here::here("R", "00_setup.R"))

lulc <- terra::rast(file.path(paths()$processed, "03_lulc_values.tif"))
fp <- terra::rast(file.path(paths()$processed, "07_floodplain.tif"))
fp <- terra::ifel(is.na(fp), 0, fp)  # punched water/foreshore cells are not demand
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
vul <- readr::read_csv(here::here("lookup", "crop_flood_vulnerability.csv"),
                       show_col_types = FALSE)
rcl <- as.matrix(vul[, c("aafc_code", "vulnerability_score")])
vul_r <- terra::ifel(in_downstream, terra::classify(lulc, rcl, others = 0), 0)
sb$ben_crops_w <- px_km2 * exactextractr::exact_extract(vul_r, sb, "sum")

# ---- 4. Critical-facility demand (schools, hospitals) ------------------------
# km²-equivalent within FAC_BUFFER_M of the floodplain ∩ downstream AOI.
#   weight_basis = "per_seat"     -> demand = weight_km2 × capacity (schools)
#   weight_basis = "per_facility" -> demand = weight_km2            (hospitals)
fac_w <- readr::read_csv(here::here("lookup", "critical_facilities.csv"),
                         show_col_types = FALSE)

read_facilities <- function(file, facility_type, capacity_col = NA_character_) {
  d <- readr::read_csv(file.path(paths()$raw, "downstream", file),
                       show_col_types = FALSE)
  d <- d[!is.na(d$LONGITUDE) & !is.na(d$LATITUDE), ]
  cap <- if (!is.na(capacity_col)) suppressWarnings(as.numeric(d[[capacity_col]])) else NA_real_
  sf::st_as_sf(
    data.frame(facility_type = facility_type, capacity = cap,
               lon = d$LONGITUDE, lat = d$LATITUDE),
    coords = c("lon", "lat"), crs = 4326
  )
}

fac <- rbind(
  read_facilities("bc_k12_schools_2026_databc.csv", "school", "DESIGN_CAPACITY_TOTAL"),
  read_facilities("hlbc_hospitals.csv", "hospital")
)

fp_dist <- terra::distance(terra::ifel(in_downstream, 1, NA))
fac <- sf::st_transform(fac, terra::crs(fp_dist))
fac$fp_dist_m <- terra::extract(fp_dist, terra::vect(fac))[, 2]
fac <- sf::st_transform(fac, sf::st_crs(downstream_aoi))
in_aoi <- lengths(sf::st_intersects(fac, downstream_aoi)) > 0
fac <- fac[in_aoi & !is.na(fac$fp_dist_m) & fac$fp_dist_m <= FAC_BUFFER_M, ]

fac <- merge(fac, fac_w, by = "facility_type", all.x = TRUE)
fac <- fac |>
  dplyr::group_by(facility_type) |>
  dplyr::mutate(capacity = dplyr::if_else(
    is.na(capacity) | capacity <= 0,
    stats::median(capacity[capacity > 0], na.rm = TRUE), capacity)) |>
  dplyr::ungroup()
fac$demand_km2 <- dplyr::if_else(fac$weight_basis == "per_seat",
                                 fac$weight_km2 * fac$capacity, fac$weight_km2)

# sum facility demand (and count) into each sub-basin
fac <- sf::st_transform(fac, sf::st_crs(sb))
hit <- sf::st_intersects(sb, fac)
sb$ben_facilities <- vapply(hit, function(i) sum(fac$demand_km2[i], na.rm = TRUE), numeric(1))
sb$n_facilities <- lengths(hit)

sb$total_demand_w <- sb$ben_built + sb$ben_crops_w
message("  · ", sum(sb$n_facilities), " facilities, ",
        round(sum(sb$ben_facilities), 1), " km²-equivalent (not in total_demand_w)")

sf::st_write(sb, file.path(paths()$processed, "10_demand_subbasin.gpkg"),
             delete_dsn = TRUE, quiet = TRUE)

# ---- 5. Same demand, split by census unit ------------------------------------
da <- read_das()
da$ben_built <- px_km2 * exactextractr::exact_extract(built_mask, da, "sum",
                                                      progress = FALSE)
da$ben_crops_w <- px_km2 * exactextractr::exact_extract(vul_r, da, "sum",
                                                        progress = FALSE)
hit_da <- sf::st_intersects(da, fac)
da$ben_facilities <- vapply(hit_da,
                            function(i) sum(fac$demand_km2[i], na.rm = TRUE),
                            numeric(1))
da$n_facilities <- lengths(hit_da)
da$da_area_km2 <- as.numeric(sf::st_area(da)) * 1e-6

exposed <- da$ben_built + da$ben_crops_w + da$ben_facilities > 0
message("  · census split: ", sum(exposed), " of ", nrow(da),
        " ", POP_UNIT, " unit(s) carry exposure (",
        round(sum(da$ben_built[exposed] + da$ben_crops_w[exposed]), 1),
        " km² of land, vs ", round(sum(sb$total_demand_w), 1),
        " km² routed per sub-basin)")
sf::st_write(da[exposed, ], file.path(paths()$processed, "10_demand_da.gpkg"),
             delete_dsn = TRUE, quiet = TRUE)

# ---- QA preview --------------------------------------------------------------
qa_png("10_demand.png", ncol = 2, function() {
  op <- graphics::par(mfrow = c(1, 2))
  on.exit(graphics::par(op), add = TRUE)
  aoi_disp <- aoi_display()

  exposed_px <- terra::ifel(ben_b > 0, 1, NA)
  view <- span_bbox(aoi_disp)
  terra::plot(exposed_px, type = "classes",
              levels = "exposed land", col = "#b30000",
              axes = TRUE, xlim = view$xlim, ylim = view$ylim,
              mar = c(3.1, 3.1, 5.2, 8.5),
              main = paste0("A. Where exposure sits: built-up land and\n",
                            "cropland inside the modelled floodplain"))
  plot(sf::st_geometry(aoi_disp), add = TRUE, border = "grey40", lwd = 1)

  terra::plot(terra::vect(sb), "total_demand_w", type = "interval",
              breaks = c(0, 0.001, 0.1, 1, 5, Inf),
              col = c("grey93", grDevices::hcl.colors(4, "Reds 3", rev = TRUE)),
              border = "grey70", lwd = 0.1, axes = TRUE,
              mar = c(3.1, 3.1, 5.2, 8.5),
              plg = list(title = "km² exposed",
                         legend = c("none", "under 0.1", "0.1 to 1",
                                    "1 to 5", "over 5")),
              main = paste0("B. Exposed land per sub-basin (km²),\n",
                            "the quantity routed downstream in step 11"))
  plot(sf::st_geometry(aoi_disp), add = TRUE, border = "black", lwd = 1.4)
})

message("✓ 10_demand.R: wrote pixel + per-sub-basin demand layers")