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
FLOOD_DEPTH_MIN_M <- 0 # min depth counted as flooded (JRC floor is 0.1 m)
BORDER_LAT <- 49 # the admin border between BC and the US
LAKE_COVER_FRAC <- 0.5 # min share of a unit a barrier lake must cover
POP_UNIT <- "DA" # census unit for benefiting population: "DA" or "DB"

# waterbodies that break routing
RESERVOIR_NAMES <- c("Capilano Lake", "Seymour Lake", "Coquitlam Lake", # dams
                     "Buntzen Lake", "Alouette Lake", "Stave Lake",
                     "Hayward Lake", "Wahleach Lake",
                     "Carpenter Lake", "Seton Lake") # BC Hydro-regulated
BARRIER_LAKE_NAMES <- c("Harrison Lake")
BARRIER_WATERBODIES <- c(RESERVOIR_NAMES, BARRIER_LAKE_NAMES)

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

# the demand AOI clipped to the coast
aoi_display <- function() {
  f <- file.path(paths()$processed, "01_aoi_downstream_land.gpkg")
  if (file.exists(f)) sf::st_read(f, quiet = TRUE) else read_aoi("downstream")
}

# bounding box to prevent cross layer clipping
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
# 2021 Census population units over the downstream AOI, at whichever geography
# POP_UNIT names. "DA" is dissemination areas, the default and what every
# published figure was built on. "DB" is dissemination blocks: 18,293 over this
# AOI against 4,139 DAs, so the floodplain split in 13 apportions over a much
# smaller unit. Block counts are noisier, and the AOI clip bites tighter, since
# a DA on the edge brings its whole population where only the blocks that
# actually touch are kept.
#
# Built from the national boundary file and the geographic attribute file,
# which counts population by block: summed to DA in DA mode, joined straight on
# in DB. Returns POP_UID (the unit's own id) and DAUID for labelling.
read_das <- function(unit = POP_UNIT) {
  unit <- match.arg(toupper(unit), c("DA", "DB"))

  pick <- function(nms, prefix, what) {
    hit <- grep(paste0("^", prefix), nms, value = TRUE, ignore.case = TRUE)[1]
    if (is.na(hit)) {
      stop("no ", what, " column (expected ", prefix,
           "...) in the StatCan download", call. = FALSE)
    }
    hit
  }

  src <- if (unit == "DA") "da_boundary_2021" else "db_boundary_2021"
  id_prefix <- if (unit == "DA") "DAUID" else "DBUID"

  # Unit geometry, clipped to the downstream AOI. The block file is 708 MB over
  # 498k features, so let GDAL do the extent clip via wkt_filter.
  f <- data_path(src)
  aoi <- sf::st_union(read_aoi("downstream"))
  lyr_crs <- sf::st_layers(f)$crs[[1]]
  filt <- sf::st_as_text(sf::st_as_sfc(sf::st_bbox(
    if (is.na(lyr_crs)) aoi else sf::st_transform(aoi, lyr_crs))))
  bnd <- sf::st_read(f, wkt_filter = filt, quiet = TRUE) |>
    sf::st_transform(PROJECT_CRS)
  bnd$POP_UID <- as.character(bnd[[pick(names(bnd), id_prefix, paste(unit, "id"))]])
  # the filter works on bounding boxes, so trim to the AOI itself
  keep <- lengths(sf::st_intersects(bnd, aoi)) > 0
  bnd <- bnd[keep, "POP_UID"]

  gaf <- readr::read_csv(data_path("da_population_2021"),
                         show_col_types = FALSE, guess_max = 50000)
  pop_tbl <- dplyr::tibble(
      POP_UID = as.character(gaf[[pick(names(gaf), id_prefix, paste(unit, "id"))]]),
      pop = suppressWarnings(
        as.numeric(gaf[[pick(names(gaf), "DBPOP", "block population")]]))) |>
    dplyr::group_by(POP_UID) |>
    dplyr::summarise(COUNT_TOTAL = sum(pop, na.rm = TRUE), .groups = "drop")

  # a unit with no block record holds no residents, so counted as zero
  out <- bnd |>
    dplyr::left_join(pop_tbl, by = "POP_UID") |>
    dplyr::mutate(COUNT_TOTAL = dplyr::coalesce(COUNT_TOTAL, 0)) |>
    sf::st_make_valid()
  # a DBUID is its DAUID plus a block suffix, so the DA falls out of the id
  out$DAUID <- if (unit == "DA") out$POP_UID else substr(out$POP_UID, 1, 8)
  message("  \u00b7 population unit = ", unit, ": ", nrow(out), " polygon(s), ",
          format(sum(out$COUNT_TOTAL), big.mark = ","), " residents in the AOI")
  out
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