# ==============================================================================
# 99_data_package.R: Publishable spatial data package
# ------------------------------------------------------------------------------
# Assembles the committed model layers into a self-describing package for
# release on an open-data portal. Nothing is modelled here: every value is read
# from data/processed and renamed, joined and documented for a reader who has
# not seen the pipeline.
#
# What the package holds:
#   1. one GeoPackage with every vector layer, in EPSG:3005
#   2. the same layers as GeoJSON (EPSG:4326) and as zipped shapefiles
#   3. attribute tables and a long scenario table as CSV
#   4. the result and input rasters as cloud-optimized GeoTIFFs
#   5. a data dictionary, a source and licence table, regional summary
#      statistics, a README, machine-readable metadata, and a file manifest
#
# Scenario-dependent quantities are suffixed with the scenario id, so one row
# describes a watershed under every rainfall scenario the pipeline has run.
#
# Inputs (data/processed/):
#   01_aoi_upstream.gpkg, 01_aoi_downstream.gpkg, 01_aoi_downstream_land.gpkg,
#   01_mvrd.gpkg, 01_fvrd.gpkg
#   02_subbasins.gpkg, 02_topology.csv
#   03_lulc_values.tif, 03_lulc_class_legend.csv, 03_lulc_natural_mask.tif
#   04_hsg.tif, 07_floodplain.tif
#   10_demand_subbasin.gpkg, 10_demand_da.gpkg, 10_demand_pixel.tif
#   11_tda_subbasin.gpkg, 12_realised_benefit_<scenario>.gpkg
#   runoff/<scenario>/09_prr_mm.tif, 09_q_baseline_mm.tif
#   data_sources.csv
#
# Outputs (output/data_package/):
#   README.md, metadata.json, manifest.csv
#   gpkg/, geojson/, shp/, csv/, raster/
# ==============================================================================

source(here::here("R", "00_setup.R"))

# ---- package identity --------------------------------------------------------
PKG_ID <- "lower_fraser_flood_attenuation"
PKG_TITLE <- "Flood-attenuating ecosystem services in the Lower Fraser"
PKG_VERSION <- format(Sys.Date(), "%Y.%m.%d")
PKG_PUBLISHER <- "Metro Vancouver"
PKG_ORIGINATOR <- "UBC Sustainability Scholars Program"
PKG_CONTACT <- "" # set before release

# The derived layers cannot be released under terms looser than the inputs they
# are built from, which are listed with their own licences in sources.csv. Set
# this once the release terms are agreed; the script warns while it is unset.
PKG_LICENCE <- "TO BE CONFIRMED"

# ---- what to emit ------------------------------------------------------------
EMIT_GEOJSON <- TRUE
EMIT_SHAPEFILE <- TRUE
EMIT_CSV <- TRUE
EMIT_RASTER <- TRUE
ZIP_PACKAGE <- FALSE # one archive of the whole package, in addition to the files

# The floodplain is the one published layer whose geometry is a derivative of a
# licensed input rather than of open data: the NHC Lower Fraser 2D model, held
# under provider terms (see sources.csv). Set FALSE to publish the exposure
# statistics derived from it without redistributing the extent itself; the
# floodplain layer and the two rasters that trace it are then withheld.
PUBLISH_FLOODPLAIN_GEOMETRY <- TRUE

GEOJSON_PRECISION <- 6 # decimal places, about 0.1 m at this latitude

# Plain-language scenario names, and short forms for the 10-character field
# names a shapefile allows.
SCENARIO_LABELS <- c(
  ar2021 = "Atmospheric river, 13-16 November 2021 (ECCC RDPA analysis)",
  wettest_month = "Climatological wettest month (PCIC PRISM 1991-2020 normals)")
SCENARIO_SHORT <- c(ar2021 = "ar21", wettest_month = "wm")

# Rasters carried in the package. Elevation, slope, curve numbers and
# precipitation are left out: they are large, and they are inputs a user can
# rebuild from the sources in sources.csv. Baseline runoff is an intermediate
# rather than a result and is off by default; it roughly doubles the package.
RASTER_SET <- tibble::tribble(
  ~src,                       ~name,                   ~role,    ~units,       ~include, ~description,
  "07_floodplain.tif",        "floodplain_extent",     "result", "1 = flooded", PUBLISH_FLOODPLAIN_GEOMETRY, "Modelled floodplain, the demand surface, after removing lakes, double-line rivers and marine foreshore.",
  "10_demand_pixel.tif",      "exposed_land",          "result", "km2 / pixel", PUBLISH_FLOODPLAIN_GEOMETRY, "Built-up land and cropland inside the modelled floodplain, as exposed area per 30 m pixel.",
  "03_lulc_natural_mask.tif", "natural_land_cover",    "input",  "1 = natural", TRUE,    "Forest, shrubland, grassland, lichen, moss and wetland: the cover replaced by bare ground in the counterfactual.",
  "03_lulc_values.tif",       "land_cover",            "input",  "class",       TRUE,    "Harmonized NALCMS land cover with AAFC crop classes. Values are defined in csv/land_cover_legend.csv.",
  "04_hsg.tif",               "hydrologic_soil_group", "input",  "class",       TRUE,    "Hydrologic soil group, A to D, driving infiltration in the curve-number method.")

# per scenario, appended to RASTER_SET below
RASTER_SCENARIO <- tibble::tribble(
  ~src,                   ~name,                ~role,    ~units, ~include, ~description,
  "09_prr_mm.tif",        "retention_mm",       "result", "mm",   TRUE,     "Storm runoff retained by natural land cover: runoff off bare ground less runoff under present cover.",
  "09_q_baseline_mm.tif", "runoff_baseline_mm", "result", "mm",   FALSE,    "Storm runoff generated under present land cover.")

# ==============================================================================
# Directories and helpers
# ==============================================================================

pkg_dir <- here::here("output", "data_package")
unlink(pkg_dir, recursive = TRUE)
new_dir <- function(...) {
  d <- file.path(...)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}
dir_gpkg <- new_dir(pkg_dir, "gpkg")
dir_json <- new_dir(pkg_dir, "geojson")
dir_shp <- new_dir(pkg_dir, "shp")
dir_csv <- new_dir(pkg_dir, "csv")
dir_ras <- new_dir(pkg_dir, "raster")
gpkg_file <- file.path(dir_gpkg, paste0(PKG_ID, ".gpkg"))

read_processed <- function(f) {
  sf::st_read(file.path(paths()$processed, f), quiet = TRUE)
}

# Watershed ids are 13-digit FWA feature ids. They travel as text so that no
# reader's spreadsheet or JSON parser rounds them.
idc <- function(x) ifelse(is.na(x), NA_character_, sprintf("%.0f", as.numeric(x)))

# Shapefiles cap field names at 10 characters and hold no logical type. The
# dictionary carries a chosen short name for every field; this only fills gaps
# and guarantees uniqueness within a layer.
shp_names <- function(nms) {
  out <- substr(nms, 1, 10)
  seen <- character(0)
  for (k in seq_along(out)) {
    stem <- out[k]
    n <- 1L
    while (out[k] %in% seen) {
      sfx <- as.character(n)
      out[k] <- paste0(substr(stem, 1, 10 - nchar(sfx)), sfx)
      n <- n + 1L
    }
    seen <- c(seen, out[k])
  }
  out
}

shp_ready <- function(x) {
  x |> dplyr::mutate(dplyr::across(dplyr::where(is.logical), as.integer))
}

# One writer per layer, so every format leaves with the same contents. The
# GeoPackage is the authoritative copy; GeoJSON is reprojected because RFC 7946
# fixes the coordinate reference system at WGS84.
written_layers <- list()
shp_field_map <- list()

write_layer <- function(x, name, description) {
  stopifnot(inherits(x, "sf"))
  x <- sf::st_make_valid(x)
  sf::st_write(x, gpkg_file, layer = name, append = FALSE, quiet = TRUE)

  if (EMIT_GEOJSON) {
    f <- file.path(dir_json, paste0(name, ".geojson"))
    suppressWarnings(sf::st_write(
      sf::st_transform(x, 4326), f, driver = "GeoJSON", delete_dsn = TRUE,
      quiet = TRUE,
      layer_options = c("RFC7946=YES",
                        paste0("COORDINATE_PRECISION=", GEOJSON_PRECISION))))
  }

  if (EMIT_SHAPEFILE) {
    tmp <- new_dir(paths()$tmp, paste0("shp_", name))
    unlink(list.files(tmp, full.names = TRUE))
    xs <- shp_ready(x)
    short <- dict_shp_names(name, setdiff(names(xs), attr(xs, "sf_column")))
    shp_field_map[[name]] <<- tibble::tibble(
      layer = name, field = names(short), shapefile_field = unname(short))
    names(xs)[match(names(short), names(xs))] <- unname(short)
    suppressWarnings(sf::st_write(xs, file.path(tmp, paste0(name, ".shp")),
                                  delete_dsn = TRUE, quiet = TRUE))
    zf <- file.path(dir_shp, paste0(name, ".zip"))
    ok <- tryCatch({
      utils::zip(zf, list.files(tmp, full.names = TRUE), flags = "-j9Xq")
      file.exists(zf)
    }, error = function(e) FALSE)
    if (!ok) message("  ! could not zip the shapefile for '", name, "'")
    unlink(tmp, recursive = TRUE)
  }

  if (EMIT_CSV) {
    readr::write_csv(sf::st_drop_geometry(x),
                     file.path(dir_csv, paste0(name, ".csv")), na = "")
  }

  written_layers[[name]] <<- tibble::tibble(
    layer = name, features = nrow(x),
    geometry = as.character(unique(sf::st_geometry_type(x)))[1],
    description = description)
  message("  · ", name, ": ", nrow(x), " features")
  invisible(x)
}

write_table <- function(x, name, description) {
  readr::write_csv(x, file.path(dir_csv, paste0(name, ".csv")), na = "")
  written_layers[[name]] <<- tibble::tibble(
    layer = name, features = nrow(x), geometry = "table",
    description = description)
  message("  · ", name, ": ", nrow(x), " rows")
  invisible(x)
}

# Rasters are re-encoded rather than copied: DEFLATE with a horizontal predictor
# roughly halves the LZW files the pipeline writes, and the overviews a COG
# carries let a portal preview them without reading the whole grid.
write_raster_cog <- function(src, out_name) {
  r <- terra::rast(src)
  out <- file.path(dir_ras, paste0(out_name, ".tif"))
  pred <- if (grepl("^FLT", terra::datatype(r)[1])) "3" else "2"
  ok <- tryCatch({
    terra::writeRaster(
      r, out, filetype = "COG", overwrite = TRUE,
      gdal = c("COMPRESS=DEFLATE", paste0("PREDICTOR=", pred),
               "BIGTIFF=IF_SAFER", "RESAMPLING=NEAREST"))
    TRUE
  }, error = function(e) FALSE)
  if (!ok) safe_writeRaster(r, out)
  message("  · ", basename(out), ": ",
          round(file.size(out) / 1e6, 1), " MB")
  invisible(out)
}

# ==============================================================================
# 1. Read the model layers
# ==============================================================================
message("Reading model layers")

sb <- read_processed("02_subbasins.gpkg")
dem <- read_processed("10_demand_subbasin.gpkg") |> sf::st_drop_geometry()
tda <- read_processed("11_tda_subbasin.gpkg") |> sf::st_drop_geometry()
topo <- readr::read_csv(file.path(paths()$processed, "02_topology.csv"),
                        show_col_types = FALSE)

# every join below is positional after a match(), so a duplicate id would move
# values onto the wrong watershed without failing
stopifnot(!anyDuplicated(sb$HYBAS_ID), !anyDuplicated(dem$HYBAS_ID),
          !anyDuplicated(tda$HYBAS_ID))

scen_files <- list.files(paths()$processed,
                         pattern = "^12_realised_benefit_.*\\.gpkg$",
                         full.names = TRUE)
scen_ids <- sub("^12_realised_benefit_(.*)\\.gpkg$", "\\1", basename(scen_files))
names(scen_files) <- scen_ids
if (!length(scen_ids)) stop("no realised-benefit layers found; run 12 first")

scen_label <- function(s) {
  if (s %in% names(SCENARIO_LABELS)) unname(SCENARIO_LABELS[[s]])
  else gsub("_", " ", s)
}
scen_abbr <- function(s) {
  if (s %in% names(SCENARIO_SHORT)) unname(SCENARIO_SHORT[[s]])
  else substr(gsub("[^a-z0-9]", "", tolower(s)), 1, 4)
}
message("  · scenarios: ", paste(scen_ids, collapse = ", "))

# ==============================================================================
# 2. Watersheds: the analysis unit, with every model result joined on
# ==============================================================================
message("Building the watershed layer")

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

# natural area is a property of the land cover, so it is identical in every
# scenario file; take it once and check the rest agree
nat <- NULL
for (s in scen_ids) {
  d <- sf::st_read(scen_files[[s]], quiet = TRUE) |> sf::st_drop_geometry()
  v <- d$natural_km2[match(ws$watershed_id, idc(d$HYBAS_ID))]
  if (is.null(nat)) nat <- v else stopifnot(isTRUE(all.equal(nat, v)))
}
ws$natural_km2 <- round(dplyr::coalesce(nat, 0), 3)
ws$natural_pct <- round(100 * ws$natural_km2 / ws$area_km2, 1)

# per-scenario columns, plus the same values in long form for anyone who would
# rather filter than read across
scen_long <- list()
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

  scen_long[[s]] <- tibble::tibble(
    watershed_id = ws$watershed_id, scenario = s, scenario_label = scen_label(s),
    retention_mm, retention_m3, runoff_m3, runoff_increase_pct,
    benefit_index, benefit_percentile, benefit_rank, benefit_status)
}
scen_long <- dplyr::bind_rows(scen_long)

# geometry last, so the attribute order above is the order a reader sees
ws <- sf::st_sf(sf::st_drop_geometry(ws), geometry = sf::st_geometry(sb))

# ==============================================================================
# 3. Downstream exposure by dissemination area
# ==============================================================================
message("Building the dissemination-area exposure layer")

da <- read_processed("10_demand_da.gpkg")
fp_r <- terra::rast(file.path(paths()$processed, "07_floodplain.tif"))
px_km2 <- (terra::xres(fp_r) * terra::yres(fp_r)) * 1e-6

# Population is apportioned by the share of each area lying in the floodplain,
# the same apportionment 13_population.R reports.
da$flood_km2 <- px_km2 * exactextractr::exact_extract(fp_r == 1, da, "sum",
                                                      progress = FALSE)
da_out <- da |>
  dplyr::transmute(
    dauid = as.character(DAUID),
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
pop_exposed_total <- sum(da$COUNT_TOTAL * pmin(1, da$flood_km2 / da$da_area_km2))

# ==============================================================================
# 4. Floodplain, areas of interest, regional districts
# ==============================================================================
message("Building the boundary layers")

fp_src <- readr::read_csv(here::here("data_sources.csv"), show_col_types = FALSE) |>
  dplyr::filter(id == "nhc_floodplain_1894")
fp_title <- if (nrow(fp_src)) fp_src$title[1] else "NHC Lower Fraser 2D flood model"

fp_bin <- terra::trim(terra::ifel(fp_r == 1, 1, NA))
fp_poly <- terra::as.polygons(fp_bin, dissolve = TRUE) |>
  sf::st_as_sf() |>
  sf::st_make_valid() |>
  sf::st_cast("MULTIPOLYGON")
floodplain_km2 <- as.numeric(sum(sf::st_area(fp_poly))) / 1e6
fp_out <- sf::st_sf(
  extent = "1894 Fraser freshet flood of record",
  source = fp_title,
  area_km2 = round(floodplain_km2, 2),
  geometry = sf::st_geometry(fp_poly))

one_poly <- function(x) sf::st_union(sf::st_geometry(x))
aoi_out <- rbind(
  sf::st_sf(role = "upstream service-providing area",
            description = paste("Assessment watersheds within",
                                MAX_FLOW_DIST_KM,
                                "km of flow distance upstream of the demand area"),
            geometry = one_poly(read_aoi("upstream"))),
  sf::st_sf(role = "downstream demand area",
            description = "Metro Vancouver and Fraser Valley regional districts",
            geometry = one_poly(read_aoi("downstream"))),
  sf::st_sf(role = "downstream demand area, coast-clipped",
            description = "The same boundary closed on the shoreline, for mapping",
            geometry = one_poly(aoi_display())))
aoi_out$area_km2 <- round(as.numeric(sf::st_area(aoi_out)) / 1e6, 1)

rd_out <- rbind(
  sf::st_sf(district = "Metro Vancouver Regional District", code = "MVRD",
            geometry = one_poly(read_processed("01_mvrd.gpkg"))),
  sf::st_sf(district = "Fraser Valley Regional District", code = "FVRD",
            geometry = one_poly(read_processed("01_fvrd.gpkg"))))
rd_out$area_km2 <- round(as.numeric(sf::st_area(rd_out)) / 1e6, 1)

# ==============================================================================
# 5. Data dictionary
# ==============================================================================
# Every published field is described here, including a shapefile-safe short
# name. Fields the dictionary misses are reported before anything is written.

dict_fixed <- tibble::tribble(
  ~layer, ~field, ~units, ~shp, ~description,
  "watersheds", "watershed_id", "", "ws_id", "BC Freshwater Atlas watershed feature id, as text. Unique key of this layer.",
  "watersheds", "downstream_watershed_id", "", "ds_id", "Watershed this one drains into. Empty at a terminal unit or a severed lake or reservoir.",
  "watersheds", "watershed_name", "", "ws_name", "Gazetted name of the watershed's main channel, where one exists.",
  "watersheds", "fwa_watershed_code", "", "fwa_code", "Freshwater Atlas watershed code. Empty for the transboundary units taken from US sources.",
  "watersheds", "delineation_source", "", "src", "Survey the watershed boundary comes from.",
  "watersheds", "admin_district", "", "district", "Regional district or state county the watershed mostly lies in.",
  "watersheds", "area_km2", "km2", "area_km2", "Watershed area.",
  "watersheds", "reach_km", "km", "reach_km", "Measured length of the watershed's main channel.",
  "watersheds", "flow_dist_km", "km", "flowdst_km", "Flow distance along the drainage network to the demand area. Zero inside it.",
  "watersheds", "stream_order", "", "strm_order", "Freshwater Atlas stream order.",
  "watersheds", "network_magnitude", "", "net_magntd", "Number of upstream first-order tributaries.",
  "watersheds", "in_demand_area", "0/1", "in_demand", "Watershed intersects the downstream demand area.",
  "watersheds", "is_terminal", "0/1", "terminal", "Watershed has no downstream link: it reaches the sea, the model edge, or a barrier.",
  "watersheds", "is_lake_barrier", "0/1", "lake_barr", "A lake of at least 10 km2, or a named reservoir, covers the unit and ends downstream routing here.",
  "watersheds", "natural_km2", "km2", "nat_km2", "Natural land cover: forest, shrubland, grassland, lichen, moss and wetland.",
  "watersheds", "natural_pct", "%", "nat_pct", "Natural land cover as a share of watershed area.",
  "watersheds", "exposed_built_km2", "km2", "expblt_km2", "Built-up land inside the modelled floodplain, within this watershed.",
  "watersheds", "exposed_crop_km2", "km2", "expcrp_km2", "Cropland inside the modelled floodplain, within this watershed.",
  "watersheds", "exposed_total_km2", "km2", "exptot_km2", "Built-up land plus cropland inside the modelled floodplain. The demand the model routes.",
  "watersheds", "facility_exposed_km2", "km2", "fac_km2", "Footprint of schools and hospitals within 250 m of the floodplain. Reported only; excluded from demand.",
  "watersheds", "facilities_n", "count", "fac_n", "Number of those facilities.",
  "watersheds", "routed_demand_km2", "km2, weighted", "routed_km2", "Downstream exposed land reachable from this watershed, weighted by an exponential distance decay with a 20 km half-life. In km2 but not a mappable area.",

  "exposure_dissemination_areas", "dauid", "", "dauid", "2021 Census dissemination area identifier.",
  "exposure_dissemination_areas", "area_km2", "km2", "area_km2", "Area of the dissemination area within the demand area.",
  "exposure_dissemination_areas", "floodplain_km2", "km2", "fp_km2", "Area of the modelled floodplain inside it.",
  "exposure_dissemination_areas", "floodplain_pct", "%", "fp_pct", "Share of the dissemination area inside the modelled floodplain.",
  "exposure_dissemination_areas", "exposed_built_km2", "km2", "expblt_km2", "Built-up land inside the modelled floodplain.",
  "exposure_dissemination_areas", "exposed_crop_km2", "km2", "expcrp_km2", "Cropland inside the modelled floodplain.",
  "exposure_dissemination_areas", "exposed_total_km2", "km2", "exptot_km2", "Built-up land plus cropland inside the modelled floodplain.",
  "exposure_dissemination_areas", "facility_exposed_km2", "km2", "fac_km2", "Footprint of schools and hospitals within 250 m of the floodplain.",
  "exposure_dissemination_areas", "facilities_n", "count", "fac_n", "Number of those facilities.",
  "exposure_dissemination_areas", "population", "residents", "pop", "2021 Census population of the dissemination area.",
  "exposure_dissemination_areas", "population_exposed", "residents", "pop_exp", "Population apportioned by the share of the area inside the floodplain. Assumes residents are spread evenly.",

  "floodplain_extent", "extent", "", "extent", "Flood scenario the extent represents.",
  "floodplain_extent", "source", "", "source", "Hydraulic model the extent comes from.",
  "floodplain_extent", "area_km2", "km2", "area_km2", "Mapped floodplain area inside the demand area, after removing lakes, double-line rivers and marine foreshore.",

  "areas_of_interest", "role", "", "role", "Which side of the assessment the boundary defines.",
  "areas_of_interest", "description", "", "descrip", "How the boundary was derived.",
  "areas_of_interest", "area_km2", "km2", "area_km2", "Area of the boundary.",

  "regional_districts", "district", "", "district", "Regional district name.",
  "regional_districts", "code", "", "code", "Short code used in the other layers.",
  "regional_districts", "area_km2", "km2", "area_km2", "Administrative area, including marine extent.")

# the per-scenario block, generated so the dictionary cannot fall behind the data
dict_scen <- purrr::map_dfr(scen_ids, function(s) {
  a <- scen_abbr(s)
  lab <- scen_label(s)
  tibble::tribble(
    ~stem, ~units, ~shp_stem, ~description,
    "retention_mm", "mm", "rmm", "Storm runoff retained per km2 of natural land.",
    "retention_m3", "m3", "rm3", "Total storm runoff retained by the watershed's natural land cover.",
    "runoff_m3", "m3", "qm3", "Storm runoff generated under present land cover.",
    "runoff_increase_pct", "%", "gpc", "Increase in storm runoff if the watershed's natural cover were replaced by bare ground.",
    "benefit_index", "index", "bix", "Retention per km2 of natural land multiplied by routed downstream demand. Relative, not a physical quantity.",
    "benefit_percentile", "percentile band", "bpc", "Percentile band of the benefit index among ranked watersheds. Empty where not ranked.",
    "benefit_rank", "rank", "brk", "Rank of the benefit index among ranked watersheds, 1 highest. Empty where not ranked.",
    "benefit_status", "", "bst", "Whether the watershed is ranked, held out as a lake or reservoir buffer, or unranked for want of natural land, runoff or downstream demand.") |>
    dplyr::transmute(
      layer = "watersheds",
      field = paste0(stem, "_", s),
      units = units,
      shp = paste0(shp_stem, "_", a),
      description = paste0(description, " Scenario: ", lab, "."))
})

dict <- dplyr::bind_rows(dict_fixed, dict_scen)

# used by write_layer() to rename fields on the way into a shapefile
dict_shp_names <- function(layer, fields) {
  d <- dict[dict$layer == layer, ]
  out <- d$shp[match(fields, d$field)]
  out[is.na(out)] <- fields[is.na(out)]
  out <- shp_names(out)
  stats::setNames(out, fields)
}

check_dictionary <- function(x, layer) {
  fields <- setdiff(names(sf::st_drop_geometry(x)), "geometry")
  miss <- setdiff(fields, dict$field[dict$layer == layer])
  if (length(miss))
    warning("undocumented field(s) in '", layer, "': ",
            paste(miss, collapse = ", "), call. = FALSE)
  extra <- setdiff(dict$field[dict$layer == layer], fields)
  if (length(extra))
    warning("dictionary describes missing field(s) in '", layer, "': ",
            paste(extra, collapse = ", "), call. = FALSE)
  invisible(TRUE)
}

for (pair in list(list(ws, "watersheds"), list(da_out, "exposure_dissemination_areas"),
                  list(fp_out, "floodplain_extent"), list(aoi_out, "areas_of_interest"),
                  list(rd_out, "regional_districts"))) {
  check_dictionary(pair[[1]], pair[[2]])
}

# ==============================================================================
# 6. Write the vector layers and tables
# ==============================================================================
message("Writing vector layers")

write_layer(ws, "watersheds",
            "Assessment watersheds with retention, downstream demand and realized benefit.")
write_layer(da_out, "exposure_dissemination_areas",
            "Flood exposure and exposed population by 2021 Census dissemination area.")
if (PUBLISH_FLOODPLAIN_GEOMETRY) {
  write_layer(fp_out, "floodplain_extent",
              "Modelled floodplain used as the demand surface.")
} else {
  message("  · floodplain_extent: withheld (PUBLISH_FLOODPLAIN_GEOMETRY is FALSE)")
}
write_layer(aoi_out, "areas_of_interest",
            "Upstream service-providing area and downstream demand area.")
write_layer(rd_out, "regional_districts",
            "The two regional districts making up the demand area.")

message("Writing tables")

write_table(scen_long, "watershed_scenario_metrics",
            "Per-watershed retention and realized benefit, one row per watershed and scenario.")
write_table(
  topo |> dplyr::transmute(watershed_id = idc(focal_id),
                           downstream_watershed_id = dplyr::if_else(
                             ds_id == 0, NA_character_, idc(ds_id)),
                           reach_km = round(reach_km, 3)),
  "watershed_topology",
  "Directed drainage links between watersheds, with measured reach lengths.")

lulc_legend <- file.path(paths()$processed, "03_lulc_class_legend.csv")
if (file.exists(lulc_legend)) {
  write_table(readr::read_csv(lulc_legend, show_col_types = FALSE),
              "land_cover_legend", "Class values used by the land-cover raster.")
}

# The registry also holds inputs this package does not derive from: the
# floodplain scenarios reserved for sensitivity work, and the designation
# overlays the interactive map draws for context. Listing them here would credit
# sources the data does not use.
SOURCE_IDS <- c("nalcms_2020", "aafc_aci", "hysogs250m", "dem_glo30",
                "prism_pcic", "prism_us", "rdpa_ar2021",
                "fwa_assessment_watersheds", "fwa_stream_networks",
                "wbd_huc12", "nhdplus_flowlines", "nhc_floodplain_1894",
                "da_boundary_2021", "da_population_2021")
sources <- readr::read_csv(here::here("data_sources.csv"), show_col_types = FALSE) |>
  dplyr::filter(id %in% SOURCE_IDS) |>
  dplyr::transmute(id, title, role = upstream_downstream, format,
                   source, attribution, license)
missing_src <- setdiff(SOURCE_IDS, sources$id)
if (length(missing_src))
  warning("source id(s) not in the registry: ", paste(missing_src, collapse = ", "),
          call. = FALSE)
write_table(sources, "sources",
            "External datasets this package derives from, with attribution and licence.")
write_table(dict, "data_dictionary", "Field definitions for every published layer.")

if (EMIT_SHAPEFILE && length(shp_field_map)) {
  write_table(dplyr::bind_rows(shp_field_map), "shapefile_field_names",
              "Field names in the shapefile copies, which allow 10 characters.")
}

# ==============================================================================
# 7. Regional summary statistics
# ==============================================================================
message("Summarising")

stat <- function(metric, value, unit, scenario = NA_character_) {
  tibble::tibble(metric = metric, scenario = scenario,
                 value = value, unit = unit)
}

summary_stats <- dplyr::bind_rows(
  stat("watersheds", nrow(ws), "count"),
  stat("watersheds delineated from the BC Freshwater Atlas",
       sum(ws$delineation_source == "BC Freshwater Atlas"), "count"),
  stat("watersheds delineated from US sources",
       sum(ws$delineation_source != "BC Freshwater Atlas"), "count"),
  stat("watersheds inside the demand area", sum(ws$in_demand_area), "count"),
  stat("watersheds severed at a lake or reservoir", sum(ws$is_lake_barrier), "count"),
  stat("upstream service-providing area", round(sum(ws$area_km2)), "km2"),
  stat("natural land cover", round(sum(ws$natural_km2)), "km2"),
  stat("natural land cover share",
       round(100 * sum(ws$natural_km2) / sum(ws$area_km2), 1), "%"),
  stat("downstream demand area",
       aoi_out$area_km2[aoi_out$role == "downstream demand area"], "km2"),
  stat("modelled floodplain", round(floodplain_km2, 1), "km2"),
  stat("exposed land", round(sum(ws$exposed_total_km2), 2), "km2"),
  stat("exposed built-up land", round(sum(ws$exposed_built_km2), 2), "km2"),
  stat("exposed cropland", round(sum(ws$exposed_crop_km2), 2), "km2"),
  stat("watersheds carrying exposure", sum(ws$exposed_total_km2 > 0), "count"),
  stat("exposed land in severed watersheds",
       round(sum(ws$exposed_total_km2[ws$is_lake_barrier]), 2), "km2"),
  stat("watersheds with non-zero routed demand", sum(ws$routed_demand_km2 > 0), "count"),
  stat("dissemination areas carrying exposure", nrow(da_out), "count"),
  stat("population of those dissemination areas", sum(da_out$population), "residents"),
  stat("population inside the modelled floodplain",
       round(pop_exposed_total), "residents"),
  purrr::map_dfr(scen_ids, function(s) {
    d <- scen_long[scen_long$scenario == s, ]
    dplyr::bind_rows(
      stat("runoff retained by natural land cover", round(sum(d$retention_m3)), "m3", s),
      stat("storm runoff under present land cover", round(sum(d$runoff_m3)), "m3", s),
      stat("increase in storm runoff if natural cover were lost",
           round(100 * sum(d$retention_m3) / sum(d$runoff_m3), 1), "%", s),
      stat("watersheds ranked", sum(d$benefit_status == "ranked"), "count", s),
      stat("watersheds held out as lake or reservoir buffers",
           sum(d$benefit_status == "lake or reservoir buffer"), "count", s),
      stat("watersheds not ranked", sum(d$benefit_status == "not ranked"), "count", s))
  }))

write_table(summary_stats, "summary_statistics",
            "Regional totals, by scenario where the quantity depends on rainfall.")

# ==============================================================================
# 8. Rasters
# ==============================================================================
ras_cols <- c("name", "role", "units", "scenario", "description")
ras <- dplyr::bind_rows(
  RASTER_SET |> dplyr::mutate(path = file.path(paths()$processed, src),
                              scenario = NA_character_),
  purrr::map_dfr(scen_ids, function(s) {
    RASTER_SCENARIO |>
      dplyr::mutate(path = file.path(paths()$processed, "runoff", s, src),
                    name = paste0(name, "_", s), scenario = s,
                    description = paste0(description, " Scenario: ",
                                         scen_label(s), "."))
  })) |>
  dplyr::filter(include, file.exists(path))

if (EMIT_RASTER && nrow(ras)) {
  message("Writing rasters")
  for (i in seq_len(nrow(ras))) write_raster_cog(ras$path[i], ras$name[i])
  write_table(ras[, ras_cols], "raster_layers",
              "The raster layers in raster/, with units and what each holds.")
} else {
  ras <- ras[0, ]
}

# ==============================================================================
# 9. README, metadata, manifest
# ==============================================================================
message("Writing documentation")

pick <- function(m, s = NA) {
  v <- summary_stats$value[summary_stats$metric == m &
                             (if (is.na(s)) is.na(summary_stats$scenario)
                              else summary_stats$scenario %in% s)]
  if (!length(v)) NA else v[1]
}
fmt <- function(x, digits = 0) formatC(x, format = "f", big.mark = ",", digits = digits)

layer_table <- dplyr::bind_rows(written_layers)
scen_lines <- paste0(
  vapply(scen_ids, function(s) as.character(glue::glue(
    "- **{s}**, {scen_label(s)}. ",
    "{fmt(pick('runoff retained by natural land cover', s) / 1e6, 1)} million m3 retained; ",
    "{pick('watersheds ranked', s)} watersheds ranked.")), character(1)),
  collapse = "\n")

readme <- glue::glue("
# {PKG_TITLE}

Version {PKG_VERSION}. Produced by the {PKG_ORIGINATOR} for {PKG_PUBLISHER}.

This package maps the flood-attenuation service that natural land cover provides
to the Lower Fraser: how much storm runoff upstream ecosystems hold back, where
downstream flood exposure sits, and which upstream watersheds are connected to
that exposure through the drainage network.

The unit of analysis is the assessment watershed. Runoff is estimated at 30 m
under present land cover and again under a bare-ground counterfactual using the
SCS curve-number method; the difference is the retention credited to natural
land cover. Downstream demand is built-up land and cropland inside the modelled
floodplain, accumulated upstream through the drainage network and discounted by
flow distance on a 20 km half-life. Realized benefit is retention per km2 of
natural land multiplied by that routed demand.

## Coverage

- {fmt(pick('watersheds'))} assessment watersheds over {fmt(pick('upstream service-providing area'))} km2, of which {fmt(pick('natural land cover'))} km2 ({pick('natural land cover share')}%) carries natural land cover
- Demand area: Metro Vancouver and the Fraser Valley Regional District, {fmt(pick('downstream demand area'))} km2
- Modelled floodplain: {pick('modelled floodplain')} km2, carrying {pick('exposed land')} km2 of exposed built-up land and cropland
- {fmt(pick('population inside the modelled floodplain'))} residents inside the modelled floodplain, in {pick('dissemination areas carrying exposure')} dissemination areas holding {fmt(pick('population of those dissemination areas'))} residents

## Rainfall scenarios

Retention and realized benefit are reported for each rainfall scenario, as
columns suffixed with the scenario id.

{scen_lines}

Rankings are far steadier across scenarios than volumes. Compare watersheds
within a scenario, and use volumes to compare between them.

## Formats

- `gpkg/{PKG_ID}.gpkg`: every vector layer, EPSG:3005 (NAD83 / BC Albers). The authoritative copy.
- `geojson/`: the same layers in EPSG:4326, for web use.
- `shp/`: the same layers as zipped shapefiles. Shapefiles allow 10-character field names; `csv/shapefile_field_names.csv` maps them back.
- `csv/`: attribute tables without geometry, the long scenario table, the drainage topology, the data dictionary, the source list and the summary statistics.
- `raster/`: cloud-optimized GeoTIFFs on the shared 30 m grid, EPSG:3005.

## Contents

{paste(sprintf('- `%s` (%s, %s): %s', layer_table$layer, layer_table$geometry,
               fmt(layer_table$features), layer_table$description), collapse = '\n')}

Raster layers:

{if (nrow(ras)) paste(sprintf('- `raster/%s.tif` (%s, %s): %s', ras$name, ras$role, ras$units, ras$description), collapse = '\n') else '- none'}

Every vector field is defined in `csv/data_dictionary.csv`, and every raster in
`csv/raster_layers.csv`.

## How to read the results

- **Retention** is what the land already does. It is a modelled comparison
  between present cover and bare ground, not a measured volume.
- **Realized benefit** combines retention with the downstream exposure a
  watershed drains toward. It is a relative index for comparing watersheds, not
  a physical quantity, and its magnitude depends on the 20 km decay half-life.
- **Routed demand** is expressed in km2 but is distance-weighted, so it is not a
  mappable area.
- Watersheds severed at a lake or reservoir keep their retention and are held
  out of the ranking, because that retention buffers the water body rather than
  reaching downstream communities.

## Limitations

- This is a regional screening product. It estimates no flood depth, extent or
  timing, and does not replace site-specific hydraulic modelling.
- Exposure comes from the Fraser freshet model domain, Hope to the Salish Sea.
  Coastal storm surge, the Serpentine and Nicomekl lowlands and small-tributary
  flooding are not represented, so the watersheds upstream of them receive no
  credit.
- Every km2 of exposed land counts equally. No depth-damage or dollar loss is
  estimated, and crops are not ranked against each other.
- Schools and hospitals are located but excluded from demand. Which assets
  belong in a critical-infrastructure set is a policy decision the model does
  not make.
- Population is reported after the ranking, to describe who lives downstream. It
  does not drive the ranking.
- Land cover is a single snapshot. No land-use change, growth or disturbance is
  modelled apart from the bare-ground comparison.
- The curve-number method estimates runoff generation only. It routes no water
  through channels, so it says nothing about peak timing or how a flood wave
  attenuates as it travels.

## Provenance and licence

Generated by `R/99_data_package.R` from the processed outputs of the analysis
pipeline in `R/`. Input datasets, with their own attribution and licence terms,
are listed in `csv/sources.csv`. The derived layers cannot be released under
terms looser than those inputs allow.

Every input is open data except the floodplain. The modelled floodplain comes
from the Northwest Hydraulic Consultants Lower Fraser River 2D flood model
(2019, for the Fraser Basin Council), held under provider terms. The
`floodplain_extent` layer and the `floodplain_extent` and `exposed_land` rasters
are derivatives of it, reprojected, rasterized to 30 m, clipped to the demand
area and with mapped water bodies removed. Confirm the provider's terms before
redistributing them. The exposure statistics elsewhere in the package are
aggregates and do not reproduce the extent.

Licence: {PKG_LICENCE}
Contact: {if (nzchar(PKG_CONTACT)) PKG_CONTACT else 'to be set before release'}
Generated: {Sys.Date()}
")
writeLines(readme, file.path(pkg_dir, "README.md"))

bb <- sf::st_bbox(sf::st_transform(aoi_out, 4326))
jsonlite::write_json(list(
  id = PKG_ID, title = PKG_TITLE, version = PKG_VERSION,
  publisher = PKG_PUBLISHER, originator = PKG_ORIGINATOR,
  contact = PKG_CONTACT, licence = PKG_LICENCE,
  generated = as.character(Sys.Date()),
  spatial_reference = list(authoritative = "EPSG:3005", geojson = "EPSG:4326"),
  resolution_m = WORKING_RES_M,
  bounding_box_wgs84 = list(xmin = bb[["xmin"]], ymin = bb[["ymin"]],
                            xmax = bb[["xmax"]], ymax = bb[["ymax"]]),
  parameters = list(distance_decay_halflife_km = HALFLIFE_KM,
                    max_flow_distance_km = MAX_FLOW_DIST_KM,
                    lake_barrier_km2 = LAKE_BARRIER_KM2,
                    facility_buffer_m = FAC_BUFFER_M),
  scenarios = lapply(scen_ids, function(s) list(id = s, label = scen_label(s))),
  layers = layer_table, rasters = ras[, ras_cols],
  lineage = "Generated by R/99_data_package.R from data/processed"
), file.path(pkg_dir, "metadata.json"), auto_unbox = TRUE, pretty = TRUE,
   digits = NA, na = "null")

files <- list.files(pkg_dir, recursive = TRUE, full.names = TRUE)
files <- files[basename(files) != "manifest.csv"]
readr::write_csv(
  tibble::tibble(
    file = sub(paste0("^", pkg_dir, "/"), "", files),
    bytes = file.size(files),
    md5 = unname(tools::md5sum(files))) |>
    dplyr::arrange(file),
  file.path(pkg_dir, "manifest.csv"))

if (ZIP_PACKAGE) {
  zf <- file.path(here::here("output"), paste0(PKG_ID, "_", PKG_VERSION, ".zip"))
  ok <- tryCatch({
    utils::zip(zf, list.files(pkg_dir, recursive = TRUE, full.names = TRUE),
               flags = "-r9Xq")
    file.exists(zf)
  }, error = function(e) FALSE)
  message(if (ok) paste0("  · archive: ", zf) else "  ! could not build the archive")
}

if (identical(PKG_LICENCE, "TO BE CONFIRMED"))
  warning("PKG_LICENCE is unset. Agree the release terms before publishing.",
          call. = FALSE)

message("✓ data package written to ", pkg_dir, " (",
        round(sum(file.size(list.files(pkg_dir, recursive = TRUE,
                                       full.names = TRUE))) / 1e6, 1), " MB)")
