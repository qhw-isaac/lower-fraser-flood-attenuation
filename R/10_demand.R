# ==============================================================================
# 10_demand.R: Service-benefiting "demand" per sub-basin
# ------------------------------------------------------------------------------
# Per-pixel "demand" = (cropland + built-up) area masked to the floodplain,
# summed per sub-basin (km²) in the following variations:
#
#   ben_both         crops + built-up
#   ben_built        built-up
#   ben_crops        crops
#   ben_crops_w      crops weighted by lookup/crop_flood_vulnerability.csv
#   ben_facilities   critical-facility demand (schools, hospitals), km²-equivalent
#                    computed but not counted in total_demand_w (see below)
#   n_facilities     facility count per sub-basin (companion to ben_facilities)
#   total_demand_w   ben_built + ben_crops_w (the routed column)
#
# Every km² of beneficiary land counts equally: the vulnerability lookup is
# uniform (every score = 1), so ben_crops_w == ben_crops. The weighted columns
# stay as a re-weighting hook, edit the lookup scores to differentiate crops.
#
# CRITICAL FACILITIES ARE NOT IN THE HEADLINE DEMAND (changed 2026-08-04).
# Schools and hospitals are still located and weighted, but ben_facilities is
# excluded from total_demand_w, so no facility weight reaches 11 -> 12 -> the
# rankings. What a school is worth against a km² of built-up land is a policy
# question, and the weights in lookup/critical_facilities.csv are placeholders,
# not valuations. The columns stay for the map's demand slider. To reinstate,
# add ben_facilities back to total_demand_w and rerun 10 -> 12.
#
# Facilities are tested against a buffer around the floodplain (FAC_BUFFER_M)
# rather than the floodplain itself. A facility is disrupted when the area
# around it floods and cuts off access, power or egress, not only when its own
# footprint goes under water.
#
# Demand outside the downstream AOI (MVRD ∪ FVRD) is zeroed, so 11 only credits
# upstream sub-basins for value at risk inside the Lower Mainland.
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
#   10_demand_da.gpkg   the same demand split by 2021 dissemination area
#                       (exposed DAs only), for reading exposure below the
#                       community level. Not routed; only the map reads it.
# ==============================================================================

source(here::here("R", "00_setup.R"))

lulc <- terra::rast(file.path(paths()$processed, "03_lulc_values.tif"))
fp <- terra::rast(file.path(paths()$processed, "07_floodplain.tif"))
# 07 stores punched cells (water bodies, foreshore) as NA, and terra::ifel()
# takes the yes-branch on an NA condition whenever a branch is a raster, so an
# NA reaching in_downstream below passed ~13.5 km² of built/crop pixels on water
# and foreshore through as demand. Punched cells are dry for demand purposes.
fp <- terra::ifel(is.na(fp), 0, fp)
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

# mask to floodplain ∩ downstream AOI. 07 already clips to the AOI; the
# re-intersection is deliberate, so 10 stays correct if 07 changes
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
# re-weighting hook: each crop class carries a flood-vulnerability score, so a
# km² of a high-loss crop could outweigh a km² of pasture. Scores are all 1
# today, so ben_crops_w == ben_crops.
vul <- readr::read_csv(here::here("lookup", "crop_flood_vulnerability.csv"),
                       show_col_types = FALSE)
# remap each crop code to its score and everything else to 0, giving a per-pixel
# weight in place of the 0/1 crop mask above
rcl <- as.matrix(vul[, c("aafc_code", "vulnerability_score")])
vul_r <- terra::ifel(in_downstream, terra::classify(lulc, rcl, others = 0), 0)
sb$ben_crops_w <- px_km2 * exactextractr::exact_extract(vul_r, sb, "sum")

# ---- 4. Critical-facility demand (schools, hospitals): reported, not routed --
# a facility takes km²-equivalent demand where it sits within FAC_BUFFER_M of
# the floodplain ∩ downstream AOI. It gets its own column and stays out of
# total_demand_w (see the header). Weights are placeholders in
# lookup/critical_facilities.csv.
#   weight_basis = "per_seat"     -> demand = weight_km2 x capacity (schools)
#   weight_basis = "per_facility" -> demand = weight_km2           (hospitals)
fac_w <- readr::read_csv(here::here("lookup", "critical_facilities.csv"),
                         show_col_types = FALSE)

# read a point-facility CSV (LONGITUDE/LATITUDE) into WGS84 sf points, pulling an
# optional capacity column used by per-seat weighting
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

# distance from every cell to the nearest floodplain ∩ downstream cell (0 inside),
# then keep facilities inside the downstream AOI and within the floodplain buffer
fp_dist <- terra::distance(terra::ifel(in_downstream, 1, NA))
fac <- sf::st_transform(fac, terra::crs(fp_dist))
fac$fp_dist_m <- terra::extract(fp_dist, terra::vect(fac))[, 2]
fac <- sf::st_transform(fac, sf::st_crs(downstream_aoi))
in_aoi <- lengths(sf::st_intersects(fac, downstream_aoi)) > 0
fac <- fac[in_aoi & !is.na(fac$fp_dist_m) & fac$fp_dist_m <= FAC_BUFFER_M, ]

# attach weights; fill a missing capacity with the type's median so a school
# with no reported size still counts
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

# headline demand covers exposed land only, so ben_facilities stays out
sb$total_demand_w <- sb$ben_built + sb$ben_crops_w
message("  · critical facilities located but NOT routed: ", sum(sb$n_facilities),
        " facilities, ", round(sum(sb$ben_facilities), 1),
        " km²-equivalent held out of total_demand_w")

sf::st_write(sb, file.path(paths()$processed, "10_demand_subbasin.gpkg"),
             delete_dsn = TRUE, quiet = TRUE)

# ---- 5. The same demand, split by dissemination area -------------------------
# the same pixels counted a second way, over census DAs, so exposure reads below
# the community level: the Township of Langley is one polygon, and its total
# hides that the exposure is all in its north. Nothing downstream consumes this
# (11 and 12 route the sub-basin columns), so it changes no result. Same three
# components, so the map's exposure weightings apply to a DA as to a community.
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

# keep only DAs carrying exposure: they are what any "where is the risk" view
# draws, and a small fraction of the whole set
exposed <- da$ben_built + da$ben_crops_w + da$ben_facilities > 0
message("  · DA split: ", sum(exposed), " of ", nrow(da),
        " dissemination areas carry exposure (",
        round(sum(da$ben_built[exposed] + da$ben_crops_w[exposed]), 1),
        " km² of land, vs ", round(sum(sb$total_demand_w), 1),
        " km² routed per sub-basin)")
sf::st_write(da[exposed, ], file.path(paths()$processed, "10_demand_da.gpkg"),
             delete_dsn = TRUE, quiet = TRUE)

# ---- QA preview --------------------------------------------------------------
# where the exposed land sits per pixel, and the per-sub-basin total 11 routes
qa_png("10_demand.png", ncol = 2, function() {
  op <- graphics::par(mfrow = c(1, 2))
  on.exit(graphics::par(op), add = TRUE)
  aoi_disp <- aoi_display()

  # ben_b holds one constant, the area of a 30 m pixel, so it is a mask rather
  # than a surface: draw the exposed pixels alone against the demand AOI
  exposed_px <- terra::ifel(ben_b > 0, 1, NA)
  # frame panel A on the demand area, not the raster's own extent (the whole
  # 30 m grid), which would leave the map a speck. Panel B keeps the full
  # provider extent, because that is its subject.
  view <- span_bbox(aoi_disp)
  # the extra top margin carries the two-line main, which terra otherwise draws
  # off the top edge of the device
  terra::plot(exposed_px, type = "classes",
              levels = "exposed land", col = "#b30000",
              axes = TRUE, xlim = view$xlim, ylim = view$ylim,
              mar = c(3.1, 3.1, 5.2, 8.5),
              main = paste0("A. Where exposure sits: built-up land and\n",
                            "cropland inside the modelled floodplain"))
  plot(sf::st_geometry(aoi_disp), add = TRUE, border = "grey40", lwd = 1)

  # the same land totalled per sub-basin, which is what 11 routes downstream.
  # Classes are named in the units they carry, and the empty class called empty
  # rather than "0 - 0.001". terra ignores `levels` when type = "interval", so
  # the labels go through plg.
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