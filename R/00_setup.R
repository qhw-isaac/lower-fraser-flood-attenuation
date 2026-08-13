# ==============================================================================
# 00_setup.R: Project Set-Up
# ==============================================================================
# Loaded by every numbered script via `source(here::here("R", "00_setup.R"))`.
# Defines constants, directory structure, and helper functions.
# ==============================================================================

# ---- Libraries ---------------------------------------------------------------
suppressPackageStartupMessages({
  # core
  library(here)           # relative paths
  library(tidyverse)
  library(glue)           # string interpolation

  # spatial
  library(sf)             # vector data
  library(terra)          # raster data
  library(exactextractr)  # fast zonal stats

  # remote data / boundaries
  library(bcmaps)         # BC administrative boundaries
  library(bcdata)         # BC Data Catalogue (FWA hydrography, etc.)
  library(nhdplusTools)   # stream-network navigation

  # network / topology
  library(igraph)         # routing graph over subbasins

  # interactive map
  library(jsonlite)
})

# ---- terra scratch space -----------------------------------------------------
# keep terra's scratch files in one temp folder for easy deletion
local({
  terra_tmp <- here::here("data", "tmp")
  dir.create(terra_tmp, showWarnings = FALSE, recursive = TRUE)
  terra::terraOptions(
    tempdir = terra_tmp,
    memfrac = 0.7,
    todisk = FALSE
  )
})

# ==============================================================================
# CONSTANTS
# ==============================================================================

PROJECT_CRS <- "EPSG:3005" # NAD83 / BC Albers
WORKING_RES_M <- 30 # project resolution (built off NALCMS)
HALFLIFE_KM <- 20 # distance-decay half-life for routing
MAX_FLOW_DIST_KM <- 100 # upstream flow-distance cap (~5 half-lives)
EDGE_BUFFER_M <- 5000 # raster-clip buffer around supply AOI
FAC_BUFFER_M <- 250 # floodplain buffer for coordinate based demand
LAKE_BARRIER_KM2 <- 10 # lake area that breaks routing (storage absorbs inflow)
LAKE_COVER_FRAC <- 0.5 # min share of a unit a barrier lake must cover

# dammed reservoirs that break routing regardless of size
RESERVOIR_NAMES <- c("Capilano Lake", "Seymour Lake", "Coquitlam Lake",
                     "Buntzen Lake", "Alouette Lake", "Stave Lake",
                     "Hayward Lake", "Wahleach Lake")

# default precipitation scenario = atmospheric river event of 2021
DEFAULT_SCENARIO <- Sys.getenv("FLOOD_SCENARIO", "ar2021")

options(scipen = 999)
sf::sf_use_s2(FALSE)

# ==============================================================================
# PATHS
# ==============================================================================

paths <- function() {
  list(
    raw = here::here("data", "raw"),
    processed = here::here("data", "processed"),
    lookup = here::here("lookup"),
    figures = here::here("output", "figures"),
    tables = here::here("output", "tables"),
    shp_out = here::here("output", "shapefiles"),
    tmp = here::here("data", "tmp")
  )
}

invisible(lapply(paths(), dir.create, showWarnings = FALSE, recursive = TRUE))

# ==============================================================================
# HELPERS: data registry
# ==============================================================================

# absolute paths for raw data files, looked up by id in data_sources.csv
data_path <- function(id, ...) {
  reg <- readr::read_csv(here::here("data_sources.csv"), show_col_types = FALSE)
  file.path(paths()$raw, reg$local_path[reg$id == id])
}

# ==============================================================================
# HELPERS: AOI / raster grid
# ==============================================================================

# load an area-of-interest polygon: "upstream" = contributing watersheds (built
# by 02), "downstream" = Lower Mainland demand area (built by 01),
# "downstream_land" = the same demand area clipped to the coast (also 01)
read_aoi <- function(role = c("upstream", "downstream", "downstream_land")) {
  role <- match.arg(role)
  sf::st_read(file.path(paths()$processed, paste0("01_aoi_", role, ".gpkg")), quiet = TRUE)
}

# the demand AOI as it should be drawn. Regional-district boundaries run far out
# into the Strait of Georgia, so on a map framed on land data the border leaves
# the panel mid-line. 01 writes a coast-clipped copy that closes on the
# shoreline. Falls back to the administrative polygon if that copy is missing.
aoi_display <- function() {
  f <- file.path(paths()$processed, "01_aoi_downstream_land.gpkg")
  if (file.exists(f)) sf::st_read(f, quiet = TRUE) else read_aoi("downstream")
}

# bounding box spanning several layers, so a map framed on one does not clip
# another drawn on top. Accepts sf and terra objects, returns xlim/ylim for
# terra::plot().
span_bbox <- function(..., pad = 0.02) {
  as_box <- function(x) {
    if (inherits(x, c("SpatRaster", "SpatVector", "SpatExtent"))) {
      v <- as.vector(terra::ext(x))
      c(xmin = v[["xmin"]], xmax = v[["xmax"]],
        ymin = v[["ymin"]], ymax = v[["ymax"]])
    } else {
      b <- sf::st_bbox(x)
      c(xmin = b[["xmin"]], xmax = b[["xmax"]],
        ymin = b[["ymin"]], ymax = b[["ymax"]])
    }
  }
  bbs <- lapply(Filter(Negate(is.null), list(...)), as_box)
  grab <- function(k, f) f(vapply(bbs, function(b) b[[k]], numeric(1)))
  x <- c(grab("xmin", min), grab("xmax", max))
  y <- c(grab("ymin", min), grab("ymax", max))
  list(xlim = x + c(-1, 1) * pad * diff(x),
       ylim = y + c(-1, 1) * pad * diff(y))
}

# empty 30 m raster over the AOI, snapped to clean WORKING_RES_M multiples
# 03_lulc.R is the only script that builds the grid from the AOI
grid_template <- function(aoi = read_aoi("upstream")) {
  e <- terra::ext(terra::vect(aoi))
  res <- WORKING_RES_M
  terra::rast(
    extent = terra::ext(floor( terra::xmin(e) / res) * res,
                        ceiling(terra::xmax(e) / res) * res,
                        floor( terra::ymin(e) / res) * res,
                        ceiling(terra::ymax(e) / res) * res),
    resolution = res, crs = PROJECT_CRS
  )
}

# reproject/resample any raster onto the shared 30 m grid
align_to_grid <- function(r, method = "near", template = NULL) {
  if (is.null(template)) {
    lulc_f <- file.path(paths()$processed, "03_lulc_values.tif")
    template <- if (file.exists(lulc_f)) terra::rast(lulc_f) else grid_template()
  }
  terra::project(r, template, method = method, align_only = FALSE)
}

# ==============================================================================
# HELPERS: census geography
# ==============================================================================
# 2021 Census dissemination areas (DAs) over the downstream AOI. DAs are the
# finest standard census geography (400-700 people), finer than the census
# subdivisions that span many watersheds.
#
# Built from two StatCan downloads: the national DA boundary file, and the
# geographic attribute file, which counts population by dissemination block and
# is summed to DA here.

read_das <- function() {
  pick <- function(nms, prefix, what) {
    hit <- grep(paste0("^", prefix), nms, value = TRUE, ignore.case = TRUE)[1]
    if (is.na(hit)) {
      stop("no ", what, " column (expected ", prefix,
           "...) in the StatCan download", call. = FALSE)
    }
    hit
  }

  # DA geometry, clipped to the downstream AOI
  bnd <- sf::st_read(data_path("da_boundary_2021"), quiet = TRUE) |>
    sf::st_transform(PROJECT_CRS)
  bnd$DAUID <- as.character(bnd[[pick(names(bnd), "DAUID", "DA id")]])
  keep <- lengths(sf::st_intersects(
    bnd, sf::st_union(read_aoi("downstream")))) > 0
  bnd <- bnd[keep, "DAUID"]

  gaf <- readr::read_csv(data_path("da_population_2021"),
                         show_col_types = FALSE, guess_max = 50000)
  pop_tbl <- dplyr::tibble(
      DAUID = as.character(gaf[[pick(names(gaf), "DAUID", "DA id")]]),
      pop = suppressWarnings(
        as.numeric(gaf[[pick(names(gaf), "DBPOP", "block population")]]))) |>
    dplyr::group_by(DAUID) |>
    dplyr::summarise(COUNT_TOTAL = sum(pop, na.rm = TRUE), .groups = "drop")

  # a DA with no block record holds no residents, so counted as zero
  bnd |>
    dplyr::left_join(pop_tbl, by = "DAUID") |>
    dplyr::mutate(COUNT_TOTAL = dplyr::coalesce(COUNT_TOTAL, 0)) |>
    sf::st_make_valid()
}

# ==============================================================================
# HELPERS: files
# ==============================================================================

# write a raster as a GeoTIFF
safe_writeRaster <- function(r, file, ...) {
  terra::writeRaster(
    r, file,
    filetype = "GTiff",
    overwrite = TRUE,
    gdal = c("COMPRESS=LZW", "TILED=YES", "BIGTIFF=IF_SAFER"),
    ...
  )
}

# ==============================================================================
# HELPERS: QA previews
# ==============================================================================
qa_dir <- function() file.path(paths()$figures, "qa")
invisible(dir.create(qa_dir(), showWarnings = FALSE, recursive = TRUE))

# render plot_fn() to a PNG in the QA folder (ncol widens the canvas for panels)
qa_png <- function(file, plot_fn, ncol = 1, panel_w = 1100, panel_h = 1050, res = 150) {
  out <- file.path(qa_dir(), file)
  grDevices::png(out, width = panel_w * ncol, height = panel_h, res = res)
  on.exit(grDevices::dev.off(), add = TRUE)
  plot_fn()
  message("  \u00b7 qa preview: ", out)
  invisible(out)
}

# ==============================================================================
message(glue::glue(
  "\u2713 setup loaded: CRS = {PROJECT_CRS}  res = {WORKING_RES_M} m  ",
  "decay half-life = {HALFLIFE_KM} km"
))