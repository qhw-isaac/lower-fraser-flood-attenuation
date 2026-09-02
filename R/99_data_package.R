# ==============================================================================
# 99_data_package.R: Spatial results package
# ------------------------------------------------------------------------------
# GeoPackage of result layers (EPSG:3005) plus input rasters
#
# Outputs (output/data_package/):
#   gpkg/lower_fraser_flood_attenuation.gpkg
#     watersheds, exposure_dissemination_areas, floodplain_extent, areas_of_interest
#   raster/
#     floodplain_extent, exposed_land, natural_land_cover, retention_mm_<scenario>
# ==============================================================================

source(here::here("R", "00_setup.R"))

PKG_ID <- "lower_fraser_flood_attenuation"

pkg_dir <- here::here("output", "data_package")
unlink(pkg_dir, recursive = TRUE)
dir_gpkg <- file.path(pkg_dir, "gpkg")
dir_ras <- file.path(pkg_dir, "raster")
dir.create(dir_gpkg, recursive = TRUE)
dir.create(dir_ras, recursive = TRUE)
gpkg_file <- file.path(dir_gpkg, paste0(PKG_ID, ".gpkg"))

read_processed <- function(f) {
  sf::st_read(file.path(paths()$processed, f), quiet = TRUE)
}

# FWA ids are 13-digit; keep them as text so nothing rounds them
idc <- function(x) ifelse(is.na(x), NA_character_, sprintf("%.0f", as.numeric(x)))

write_layer <- function(x, name) {
  x <- sf::st_make_valid(x)
  sf::st_write(x, gpkg_file, layer = name, append = FALSE, quiet = TRUE)
  message("  · ", name, ": ", nrow(x), " features")
  invisible(x)
}

write_raster_cog <- function(src, out_name) {
  r <- terra::rast(src)
  out <- file.path(dir_ras, paste0(out_name, ".tif"))
  pred <- if (grepl("^FLT", terra::datatype(r)[1])) "3" else "2"
  terra::writeRaster(
    r, out, filetype = "COG", overwrite = TRUE,
    gdal = c("COMPRESS=DEFLATE", paste0("PREDICTOR=", pred),
             "BIGTIFF=IF_SAFER", "RESAMPLING=NEAREST"))
  message("  · ", basename(out), ": ", round(file.size(out) / 1e6, 1), " MB")
}

# ---- Read -------------------------------------------------------------------
sb <- read_processed("02_subbasins.gpkg")
dem <- read_processed("10_demand_subbasin.gpkg") |> sf::st_drop_geometry()
tda <- read_processed("11_tda_subbasin.gpkg") |> sf::st_drop_geometry()

scen_files <- list.files(paths()$processed,
                         pattern = "^12_realised_benefit_.*\\.gpkg$",
                         full.names = TRUE)
scen_ids <- sub("^12_realised_benefit_(.*)\\.gpkg$", "\\1", basename(scen_files))
names(scen_files) <- scen_ids
if (!length(scen_ids)) stop("no realised-benefit layers found; run 12 first")
message("  · scenarios: ", paste(scen_ids, collapse = ", "))

# ---- Watersheds -------------------------------------------------------------
d_idx <- match(sb$HYBAS_ID, dem$HYBAS_ID)
t_idx <- match(sb$HYBAS_ID, tda$HYBAS_ID)

ws <- sb |>
  dplyr::transmute(
    watershed_id = idc(HYBAS_ID),
    downstream_watershed_id = dplyr::if_else(
      is.na(NEXT_DOWN) | NEXT_DOWN == 0 | is_sink %in% TRUE,
      NA_character_, idc(NEXT_DOWN)),
    watershed_name = as.character(GNIS_NAME_1),
    fwa_watershed_code = as.character(FWA_WATERSHED_CODE),
    delineation_source = dplyr::if_else(is.na(FWA_WATERSHED_CODE),
                                        "USGS Watershed Boundary Dataset",
                                        "BC Freshwater Atlas"),
    admin_district = as.character(admin_district),
    area_km2 = round(SUB_AREA, 3),
    reach_km = round(reach_km, 3),
    flow_dist_km = round(flow_dist_km, 3),
    stream_order = as.integer(WATERSHED_ORDER),
    network_magnitude = as.integer(WATERSHED_MAGNITUDE),
    in_demand_area = in_downstream_aoi %in% TRUE,
    is_terminal = is_sink %in% TRUE,
    is_lake_barrier = is_lake_barrier %in% TRUE,
    exposed_built_km2 = round(dplyr::coalesce(dem$ben_built[d_idx], 0), 6),
    exposed_crop_km2 = round(dplyr::coalesce(dem$ben_crops[d_idx], 0), 6),
    exposed_total_km2 = round(dplyr::coalesce(dem$ben_both[d_idx], 0), 6),
    facility_exposed_km2 = round(dplyr::coalesce(dem$ben_facilities[d_idx], 0), 6),
    facilities_n = as.integer(dplyr::coalesce(dem$n_facilities[d_idx], 0)),
    routed_demand_km2 = round(dplyr::coalesce(tda$tda_total_w[t_idx], 0), 4))

nat <- NULL
for (s in scen_ids) {
  d <- sf::st_read(scen_files[[s]], quiet = TRUE) |> sf::st_drop_geometry()
  v <- d$natural_km2[match(ws$watershed_id, idc(d$HYBAS_ID))]
  if (is.null(nat)) nat <- v
}
ws$natural_km2 <- round(dplyr::coalesce(nat, 0), 3)
ws$natural_pct <- round(100 * ws$natural_km2 / ws$area_km2, 1)

for (s in scen_ids) {
  d <- sf::st_read(scen_files[[s]], quiet = TRUE) |> sf::st_drop_geometry()
  dd <- d[match(ws$watershed_id, idc(d$HYBAS_ID)), ]

  retention_mm <- round(dplyr::coalesce(dd$prr_per_nat_mm_km2, 0), 2)
  retention_m3 <- round(dplyr::coalesce(dd$prr_total_mm_km2, 0) * 1000)
  runoff_m3 <- round(dplyr::coalesce(dd$q_baseline_mm_km2, 0) * 1000)
  runoff_increase_pct <- round(ifelse(
    dplyr::coalesce(dd$q_baseline_mm_km2, 0) > 0,
    100 * dplyr::coalesce(dd$prr_total_mm_km2, 0) / dd$q_baseline_mm_km2, 0), 1)
  benefit_index <- round(dplyr::coalesce(dd$ri_index, 0), 4)
  benefit_status <- dplyr::case_when(
    dd$ri_interval %in% "lake_buffer" ~ "lake or reservoir buffer",
    !is.na(dd$ri_interval) ~ "ranked",
    TRUE ~ "not ranked")
  benefit_percentile <- dplyr::if_else(benefit_status == "ranked",
                                       gsub("_", "-", dd$ri_interval),
                                       NA_character_)
  benefit_rank <- rep(NA_integer_, nrow(ws))
  elig <- which(benefit_status == "ranked")
  benefit_rank[elig] <- as.integer(rank(-benefit_index[elig], ties.method = "min"))

  ws[[paste0("retention_mm_", s)]] <- retention_mm
  ws[[paste0("retention_m3_", s)]] <- retention_m3
  ws[[paste0("runoff_m3_", s)]] <- runoff_m3
  ws[[paste0("runoff_increase_pct_", s)]] <- runoff_increase_pct
  ws[[paste0("benefit_index_", s)]] <- benefit_index
  ws[[paste0("benefit_percentile_", s)]] <- benefit_percentile
  ws[[paste0("benefit_rank_", s)]] <- benefit_rank
  ws[[paste0("benefit_status_", s)]] <- benefit_status
}

ws <- sf::st_sf(sf::st_drop_geometry(ws), geometry = sf::st_geometry(sb))

# ---- Dissemination-area exposure --------------------------------------------
da <- read_processed("10_demand_da.gpkg")
fp_r <- terra::rast(file.path(paths()$processed, "07_floodplain.tif"))
px_km2 <- (terra::xres(fp_r) * terra::yres(fp_r)) * 1e-6
da$flood_km2 <- px_km2 * exactextractr::exact_extract(fp_r == 1, da, "sum",
                                                      progress = FALSE)
da_out <- da |>
  dplyr::transmute(
    dauid = as.character(POP_UID),
    area_km2 = round(da_area_km2, 4),
    floodplain_km2 = round(flood_km2, 4),
    floodplain_pct = round(100 * pmin(1, flood_km2 / da_area_km2), 1),
    exposed_built_km2 = round(ben_built, 6),
    exposed_crop_km2 = round(ben_crops_w, 6),
    exposed_total_km2 = round(ben_built + ben_crops_w, 6),
    facility_exposed_km2 = round(ben_facilities, 6),
    facilities_n = as.integer(n_facilities),
    population = as.integer(COUNT_TOTAL),
    population_exposed = round(COUNT_TOTAL * pmin(1, flood_km2 / da_area_km2)))

# ---- Floodplain and study bounds --------------------------------------------
fp_poly <- terra::as.polygons(terra::trim(terra::ifel(fp_r == 1, 1, NA)),
                              dissolve = TRUE) |>
  sf::st_as_sf() |>
  sf::st_make_valid() |>
  sf::st_cast("MULTIPOLYGON")
fp_out <- sf::st_sf(
  extent = "River flood, 100-year return period, no flood protection",
  source = "JRC/CEMS global river flood hazard map, RP100",
  area_km2 = round(as.numeric(sum(sf::st_area(fp_poly))) / 1e6, 2),
  geometry = sf::st_geometry(fp_poly))

one_poly <- function(x) sf::st_union(sf::st_geometry(x))
aoi_out <- rbind(
  sf::st_sf(role = "upstream service-providing area",
            geometry = one_poly(read_aoi("upstream"))),
  sf::st_sf(role = "downstream demand area",
            geometry = one_poly(read_aoi("downstream"))),
  sf::st_sf(role = "downstream demand area, coast-clipped",
            geometry = one_poly(aoi_display())))
aoi_out$area_km2 <- round(as.numeric(sf::st_area(aoi_out)) / 1e6, 1)

# ---- Write vectors ----------------------------------------------------------
write_layer(ws, "watersheds")
write_layer(da_out, "exposure_dissemination_areas")
write_layer(fp_out, "floodplain_extent")
write_layer(aoi_out, "areas_of_interest")

# ---- Rasters ----------------------------------------------------------------
rasters <- c(
  floodplain_extent = file.path(paths()$processed, "07_floodplain.tif"),
  exposed_land = file.path(paths()$processed, "10_demand_pixel.tif"),
  natural_land_cover = file.path(paths()$processed, "03_lulc_natural_mask.tif"))
for (s in scen_ids) {
  rasters[[paste0("retention_mm_", s)]] <-
    file.path(paths()$processed, "runoff", s, "09_prr_mm.tif")
}
message("Writing rasters")
for (nm in names(rasters)) {
  if (file.exists(rasters[[nm]])) write_raster_cog(rasters[[nm]], nm)
}

message("✓ data package written to ", pkg_dir, " (",
        round(sum(file.size(list.files(pkg_dir, recursive = TRUE,
                                       full.names = TRUE))) / 1e6, 1), " MB)")
