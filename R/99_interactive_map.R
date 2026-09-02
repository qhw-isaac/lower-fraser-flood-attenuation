# ==============================================================================
# 99_interactive_map.R: Web data for the upstream-ecosystems interactive map
# ------------------------------------------------------------------------------
# Builds the data behind output/interactive_map/index.html
#
# Inputs (data/processed/):
#   02_subbasins.gpkg, 10_demand_subbasin.gpkg, 10_demand_da.gpkg,
#   11_tda_subbasin.gpkg, 12_realised_benefit_<scenario>.gpkg
#   01_aoi_downstream.gpkg, 01_mvrd.gpkg, 01_fvrd.gpkg
#   bcmaps::census_subdivision()
#   optional overlays: alr_polygons / rnin_patch / rnin_corridor
#
# Outputs (output/interactive_map/data.js):
#   META, SUBBASINS, MUNIS, DAS, AOI, FIRSTNATIONS, LAKES, RIVERS, STREAMS, ALR,
#   RNIN_PATCH, RNIN_CORRIDOR, GRAPH, MUNI_NAMES, MUNI_BASINS, MUNI_DEMAND,
#   DEMAND, DA_DEMAND, RETAIN, RETAIN_MM, RETAIN_PCT, NATTYPES
# ==============================================================================

source(here::here("R", "00_setup.R"))

SCENARIO <- DEFAULT_SCENARIO
SIMPLIFY_M <- 60 # geometry simplification tolerance (m) for web payload
MIN_EXPOSURE_KM2 <- 0.001 # a member is served by a basin from this much exposure

out_dir <- here::here("output", "interactive_map")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- 1. load model layers ----------------------------------------------------
read_processed <- function(f) {
  sf::st_read(file.path(paths()$processed, f), quiet = TRUE)
}

sb <- read_processed("02_subbasins.gpkg")
tda <- read_processed("11_tda_subbasin.gpkg") |> sf::st_drop_geometry()
ri <- read_processed(glue::glue("12_realised_benefit_{SCENARIO}.gpkg")) |>
  sf::st_drop_geometry()

# components for the map's demand-mix sliders
dem <- read_processed("10_demand_subbasin.gpkg") |> sf::st_drop_geometry()

downstream_aoi <- read_aoi("downstream")

idc <- function(x) sprintf("%.0f", as.numeric(x))

sb$id <- idc(sb$HYBAS_ID)
sb$next_id <- idc(sb$NEXT_DOWN)
sb$next_id[sb$NEXT_DOWN == 0 | sb$is_sink] <- NA_character_

# merge in the model metrics the map surfaces
ri$id <- idc(ri$HYBAS_ID)
tda$id <- idc(tda$HYBAS_ID)
sb <- sb |>
  dplyr::left_join(
    ri |> dplyr::select(id, natural_km2, prr_total_mm_km2, prr_per_nat_mm_km2,
                        q_baseline_mm_km2),
    by = "id") |>
  dplyr::left_join(tda |> dplyr::select(id, own_demand, tda_total_w), by = "id")

is_us <- if ("FWA_WATERSHED_CODE" %in% names(sb))
  is.na(sb$FWA_WATERSHED_CODE) else rep(FALSE, nrow(sb))

# ---- 1a. Clip US/FWA overlap -------------------------------------------------
if (any(is_us) && any(!is_us)) {
  s2_prev <- sf::sf_use_s2(FALSE)
  bc_union <- sf::st_union(sf::st_geometry(sb[!is_us, ]))
  us_geom <- suppressWarnings(
    sf::st_difference(sf::st_geometry(sb[is_us, ]), bc_union))
  before <- as.numeric(sf::st_area(sf::st_geometry(sb[is_us, ]))) / 1e6
  after <- as.numeric(sf::st_area(us_geom)) / 1e6
  SLIVER_KM2 <- 5
  trimmed <- (before - after) > 0.01
  if (any(trimmed)) {
    g <- sf::st_geometry(sb)
    g[which(is_us)] <- sf::st_make_valid(us_geom)
    sf::st_geometry(sb) <- g
    sb$SUB_AREA[is_us] <- after
    if (sum(before - after) > SLIVER_KM2)
      warning(sprintf("%.1f km2 of US/FWA overlap", sum(before - after)),
              immediate. = TRUE)
  }
  sf::sf_use_s2(s2_prev)
}

# ---- 1c. Natural cover-type composition per sub-basin ------------------------
lulc_r <- terra::rast(file.path(paths()$processed, "03_lulc_values.tif"))
nat_grp <- c(`1` = "forest", `2` = "forest", `5` = "forest", `6` = "forest",
             `8` = "shrub", `11` = "shrub", `10` = "grass", `12` = "grass",
             `13` = "alpine", `14` = "wetland")
lulc_cell_km2 <- prod(terra::res(lulc_r)) / 1e6
nat_comp <- lapply(
  exactextractr::exact_extract(lulc_r, sf::st_transform(sb, terra::crs(lulc_r)),
                               progress = FALSE),
  function(df) {
    df <- df[!is.na(df$value) & as.character(df$value) %in% names(nat_grp), ]
    if (!nrow(df)) return(NULL)
    a <- tapply(df$coverage_fraction, nat_grp[as.character(df$value)], sum)
    as.list(round(a * lulc_cell_km2, 1))
  })
names(nat_comp) <- sb$id
nat_comp <- nat_comp[!vapply(nat_comp, is.null, logical(1))]
message("  · natural cover composition for ", length(nat_comp), " sub-basins")

# ---- 2. Member jurisdictions across the downstream AOI (MVRD ∪ FVRD) ---------
csd_all <- bcmaps::census_subdivision() |>
  sf::st_transform(PROJECT_CRS) |>
  sf::st_make_valid()

# Treaty First Nations (e.g. Tsawwassen) have a blank TYPE_DESC; fall back to code
tdesc <- csd_all$CENSUS_SUBDIVISION_TYPE_DESC
tcode <- csd_all$CENSUS_SUBDIVISION_TYPE_CODE
csd_all$member_type <- dplyr::case_when(
  tdesc == "City" ~ "City",
  tdesc == "District Municipality" ~ "District",
  tdesc == "Town" ~ "Town",
  tdesc == "Village" ~ "Village",
  tdesc == "Island Municipality" ~ "Island Muni.",
  tdesc == "Regional Municipality" ~ "Regional Muni.",
  tdesc == "Regional District Electoral Area" ~ "Electoral Area",
  tdesc == "Indian Reserve" ~ "Reserve",
  tcode %in% c("TWL", "NL", "TAL", "IGD", "S-É") ~ "Treaty Land",
  TRUE ~ "Other")
csd_all$member_class <- dplyr::case_when(
  csd_all$member_type %in% c("City", "District", "Town", "Village",
                             "Island Muni.", "Regional Muni.") ~ "Municipality",
  csd_all$member_type == "Electoral Area" ~ "Electoral Area",
  csd_all$member_type %in% c("Reserve", "Treaty Land") ~ "First Nation",
  TRUE ~ "Other")

# tag each CSD by regional district (point-on-surface; drop those outside MVRD/FVRD)
mvrd_poly <- read_processed("01_mvrd.gpkg") |> sf::st_geometry() |> sf::st_union()
fvrd_poly <- read_processed("01_fvrd.gpkg") |> sf::st_geometry() |> sf::st_union()
ctr <- suppressWarnings(sf::st_point_on_surface(csd_all))
in_m <- as.logical(sf::st_within(ctr, mvrd_poly, sparse = FALSE)[, 1])
in_f <- as.logical(sf::st_within(ctr, fvrd_poly, sparse = FALSE)[, 1])
csd_all$rd <- dplyr::case_when(in_m ~ "MVRD", in_f ~ "FVRD", TRUE ~ NA_character_)

muni <- csd_all[!is.na(csd_all$rd), ]
muni$name <- muni$CENSUS_SUBDIVISION_NAME
muni <- muni |>
  dplyr::group_by(name, rd, member_type, member_class) |>
  dplyr::summarise(.groups = "drop") |>
  sf::st_make_valid()

fn_lands <- muni[muni$member_class == "First Nation", c("name", "member_type", "rd")]
message("  · First Nation lands layer: ", nrow(fn_lands), " CSD(s) (",
        sum(fn_lands$rd == "MVRD"), " MVRD / ", sum(fn_lands$rd == "FVRD"), " FVRD)")

# dissolve reserves into the member containing them (avoids hole rings in Leaflet)
is_res <- muni$member_type == "Reserve"
if (any(is_res) && any(!is_res)) {
  host_rows <- which(!is_res); res_rows <- which(is_res)
  nn <- sf::st_nearest_feature(muni[res_rows, ], muni[host_rows, ])
  g <- sf::st_geometry(muni)
  for (i in seq_along(res_rows)) {
    gi <- host_rows[nn[i]]
    g[gi] <- sf::st_union(g[c(gi, res_rows[i])])
  }
  sf::st_geometry(muni) <- sf::st_make_valid(g)
  muni <- muni[host_rows, ]           # drop the now-absorbed reserves
  message("  · dissolved ", length(res_rows), " Indian Reserve(s) into host members")
}

# display label: bare name, with a type suffix only where names collide
dup <- muni$name %in% muni$name[duplicated(muni$name)]
muni$muni <- ifelse(dup, paste0(muni$name, " (", muni$member_type, ")"), muni$name)
muni$muni_id <- seq_len(nrow(muni))

# ---- 2b. Snap the community layer to the regional-district boundary ----------
# clip to the RD, then fold uncovered slivers into the neighbour sharing the
# longest boundary
best_host <- function(sliver, g, cand) {
  if (length(cand) < 2) return(cand)
  shared <- vapply(cand, function(k) sum(as.numeric(sf::st_length(
    suppressWarnings(sf::st_intersection(sf::st_boundary(sliver), g[k]))))),
    numeric(1))
  cand[which.max(shared)]
}

rd_union <- sf::st_union(c(mvrd_poly, fvrd_poly))
muni <- suppressWarnings(sf::st_intersection(muni, rd_union)) |> sf::st_make_valid()
gaps <- suppressWarnings(sf::st_difference(rd_union, sf::st_union(sf::st_geometry(muni))))
gaps <- suppressWarnings(sf::st_cast(sf::st_make_valid(gaps), "POLYGON"))
gaps <- gaps[as.numeric(sf::st_area(gaps)) > 1e3]
if (length(gaps) > 0) {
  g <- sf::st_geometry(muni)
  hits <- sf::st_intersects(gaps, muni)
  n_fold <- 0L
  for (i in seq_along(gaps)) {
    host <- best_host(gaps[i], g, hits[[i]])
    if (!length(host)) next
    g[host] <- sf::st_union(g[host], gaps[i])
    n_fold <- n_fold + 1L
  }
  sf::st_geometry(muni) <- sf::st_make_valid(g)
  message("  · snapped communities to the RD boundary (clipped overhangs; folded ",
          n_fold, " boundary sliver(s) into their neighbours)")
}
message("  · ", nrow(muni), " member jurisdictions in downstream AOI (",
        sum(muni$rd == "MVRD"), " MVRD / ", sum(muni$rd == "FVRD"), " FVRD)")

# ---- 3. Link demand basins to jurisdictions ----------------------------------
# a basin serves a member where they overlap and the basin has floodplain demand.
# drop lake barriers and units whose downstream chain never reaches the demand area.
drains_to_demand <- local({
  nxt <- stats::setNames(sb$next_id, sb$id)
  in_aoi <- stats::setNames(sb$in_downstream_aoi %in% TRUE, sb$id)
  vapply(sb$id, function(start) {
    cur <- start
    seen <- character(0)
    while (!is.na(cur) && cur %in% names(in_aoi) && !cur %in% seen) {
      if (in_aoi[[cur]]) return(TRUE)
      seen <- c(seen, cur)
      cur <- nxt[[cur]]
    }
    FALSE
  }, logical(1))
})
if (any(!drains_to_demand)) {
  message("  · ", sum(!drains_to_demand), " unit(s) drain out of the region ",
          "(no path to the demand area); not credited as providers")
}

demand_basins <- sb$id[!is.na(sb$own_demand) & sb$own_demand > 0 &
                         !sb$is_lake_barrier &
                         drains_to_demand[sb$id]]
dem$id <- idc(dem$HYBAS_ID)

# split demand basins along community boundaries; link only where the piece holds exposure
inter <- suppressWarnings(
  sf::st_intersection(
    sb[sb$id %in% demand_basins, c("id")],
    muni[, c("muni_id", "muni")]
  )
) |> sf::st_make_valid()
if (any(sf::st_geometry_type(inter) == "GEOMETRYCOLLECTION"))
  inter <- sf::st_collection_extract(inter, "POLYGON", warn = FALSE)
inter <- inter[sf::st_geometry_type(inter) %in% c("POLYGON", "MULTIPOLYGON"), ]
dem_pix <- terra::rast(file.path(paths()$processed, "10_demand_pixel.tif"))
inter$exp_km2 <- exactextractr::exact_extract(dem_pix, inter, "sum",
                                              progress = FALSE)
basin_muni <- inter |>
  sf::st_drop_geometry() |>
  dplyr::filter(exp_km2 >= MIN_EXPOSURE_KM2)
message("  · basin-member links: ", nrow(basin_muni), " of ", nrow(inter),
        " overlaps carry exposure (", nrow(inter) - nrow(basin_muni),
        " dropped as overlap without exposure)")

# basin id -> integer muni ids it directly contains demand for
direct_munis <- split(basin_muni$muni_id, basin_muni$id)

# apportion each basin's demand components by the member's share of that basin's exposure
muni_dem <- basin_muni |>
  dplyr::group_by(id) |>
  dplyr::mutate(share = exp_km2 / sum(exp_km2)) |>
  dplyr::ungroup() |>
  dplyr::left_join(dplyr::select(dem, id, ben_built, ben_crops_w), by = "id") |>
  dplyr::group_by(muni_id) |>
  dplyr::summarise(built = sum(share * ben_built),
                   crops = sum(share * ben_crops_w),
                   .groups = "drop")

dem_out <- dem[dem$id %in% demand_basins,
               c("id", "ben_built", "ben_crops_w", "ben_facilities")]
demand_json <- jsonlite::toJSON(
  stats::setNames(
    lapply(seq_len(nrow(dem_out)), function(i) list(
      built = round(dem_out$ben_built[i], 3),
      crops = round(dem_out$ben_crops_w[i], 3),
      fac   = round(dem_out$ben_facilities[i], 3))),
    dem_out$id),
  auto_unbox = TRUE)
message("  · demand components for ", nrow(dem_out), " demand basin(s)")

# ---- 4. Flow graph: ancestors (upstream) / descendants (downstream) ----------
edges <- sb |>
  sf::st_drop_geometry() |>
  dplyr::filter(!is.na(next_id), next_id %in% sb$id) |>
  dplyr::select(from = id, to = next_id)
g <- igraph::graph_from_data_frame(edges, directed = TRUE,
                                   vertices = data.frame(name = sb$id))

up_list <- vector("list", nrow(sb)); names(up_list) <- sb$id
contributes <- vector("list", nrow(sb)); names(contributes) <- sb$id

barrier_ids <- sb$id[sb$is_lake_barrier]

for (v in sb$id) {
  desc <- setdiff(names(igraph::subcomponent(g, v, mode = "out")), v)
  up_list[[v]] <- setdiff(names(igraph::subcomponent(g, v, mode = "in")), v)
  reach <- setdiff(c(v, desc), barrier_ids)
  pm <- unique(unlist(direct_munis[intersect(reach, names(direct_munis))]))
  contributes[[v]] <- sort(as.integer(pm))
}

sb$n_contrib <- vapply(contributes, length, integer(1))
sb$n_up <- vapply(up_list, length, integer(1))

# split contributor counts by regional district
rd_of_id <- stats::setNames(muni$rd, as.character(muni$muni_id))
sb$n_contrib_mvrd <- vapply(contributes,
  function(p) sum(rd_of_id[as.character(p)] == "MVRD"), integer(1))
sb$n_contrib_fvrd <- vapply(contributes,
  function(p) sum(rd_of_id[as.character(p)] == "FVRD"), integer(1))

# ---- 4b. Drop cross-border units that serve nobody and hold no demand-area ground
local({
  bx <- lapply(sf::st_geometry(sf::st_transform(sb, 4326)), sf::st_bbox)
  crosses <- vapply(bx, function(b) b[["ymin"]] < BORDER_LAT - 0.005,
                    logical(1))
  serves_nobody <- sb$n_contrib == 0 & !(sb$id %in% names(direct_munis))

  # test geometry, not in_downstream_aoi (centroid can sit south of the border)
  dom <- sf::st_union(downstream_aoi)
  near <- lengths(sf::st_intersects(sb, dom)) > 0
  in_demand <- rep(FALSE, nrow(sb))
  in_demand[near] <- vapply(which(near), function(i) {
    g <- suppressWarnings(sf::st_intersection(sf::st_geometry(sb)[i], dom))
    length(g) > 0 && sum(as.numeric(sf::st_area(g))) / 1e6 > 0.01
  }, logical(1))

  inert <- serves_nobody & !in_demand
  drop <- crosses & inert
  if (any(drop)) {
    message("  \u00b7 dropped ", sum(drop), " cross-border unit(s) delivering ",
            "nothing and holding no exposure (",
            sprintf("%.0f", sum(sb$SUB_AREA[drop], na.rm = TRUE)), " km2)")
    sb <<- sb[!drop, ]
  }
})
contributes <- contributes[names(contributes) %in% sb$id]
up_list <- up_list[names(up_list) %in% sb$id]
nat_comp <- nat_comp[names(nat_comp) %in% sb$id]
is_us <- if ("FWA_WATERSHED_CODE" %in% names(sb))
  is.na(sb$FWA_WATERSHED_CODE) else rep(FALSE, nrow(sb))

# ---- 5. Build web geometries (WGS84, simplified) -----------------------------

# reproject + simplify. coverage = TRUE uses rmapshaper so shared borders stay aligned
to_web <- function(x, coverage = FALSE) {
  x <- sf::st_transform(x, PROJECT_CRS)
  if (any(sf::st_geometry_type(x) == "GEOMETRYCOLLECTION"))
    x <- sf::st_collection_extract(x, "POLYGON", warn = FALSE)
  if (coverage && requireNamespace("rmapshaper", quietly = TRUE)) {
    x <- rmapshaper::ms_simplify(x, keep = 0.12, keep_shapes = TRUE,
                                 explode = FALSE)
  } else {
    x <- sf::st_simplify(x, dTolerance = SIMPLIFY_M, preserveTopology = TRUE)
  }
  x |>
    sf::st_transform("EPSG:4326") |>
    sf::st_make_valid()
}

add_centroids <- function(x) {
  pt <- suppressWarnings(sf::st_point_on_surface(sf::st_geometry(x)))
  co <- sf::st_coordinates(pt)
  x$cx <- round(co[, 1], 5)
  x$cy <- round(co[, 2], 5)
  x
}

sb$ws_name <- if ("GNIS_NAME_1" %in% names(sb)) as.character(sb$GNIS_NAME_1) else NA_character_

# ---- 4c. Close pinholes in the watershed coverage (display only) -------------
HOLE_FILL_KM2 <- 1

local({
  s2_prev <- sf::sf_use_s2(FALSE)
  on.exit(sf::sf_use_s2(s2_prev), add = TRUE)

  # interior rings of the dissolved coverage = gaps between watersheds
  u <- sf::st_make_valid(sf::st_union(sf::st_geometry(sb)))
  rings <- list()
  for (poly in sf::st_cast(u, "POLYGON", warn = FALSE)) {
    if (length(poly) > 1)
      for (k in 2:length(poly)) rings[[length(rings) + 1]] <- sf::st_polygon(poly[k])
  }
  if (!length(rings)) return(invisible(NULL))

  gaps <- sf::st_sfc(rings, crs = sf::st_crs(sb))
  gaps <- gaps[as.numeric(sf::st_area(gaps)) / 1e6 < HOLE_FILL_KM2]
  if (!length(gaps)) return(invisible(NULL))

  g <- sf::st_geometry(sb)
  hits <- sf::st_intersects(gaps, sb)
  n <- 0L
  for (i in seq_along(gaps)) {
    host <- best_host(gaps[i], g, hits[[i]])
    if (!length(host)) next
    g[host] <- sf::st_union(g[host], gaps[i])
    n <- n + 1L
  }
  sf::st_geometry(sb) <<- sf::st_make_valid(g)
  message("  · closed ", n, " pinhole(s) in the watershed coverage (display only)")
})

sb_web <- to_web(sb, coverage = TRUE) |> add_centroids()

sb_out <- sb_web |>
  dplyr::transmute(
    id, next_id,
    name = ws_name,
    natural_km2 = round(natural_km2, 1),
    retain_m3 = round(dplyr::coalesce(prr_total_mm_km2, 0) * 1000),
    retain_mm = round(dplyr::coalesce(prr_per_nat_mm_km2, 0), 1),
    runoff_gain_pct = round(
      ifelse(dplyr::coalesce(q_baseline_mm_km2, 0) > 0,
             100 * dplyr::coalesce(prr_total_mm_km2, 0) / q_baseline_mm_km2, 0)),
    reach_km = round(dplyr::coalesce(reach_km, 0), 2),
    n_contrib, n_contrib_mvrd, n_contrib_fvrd, n_up,
    barrier = is_lake_barrier,
    barrier_type = dplyr::case_when(
      !is_lake_barrier ~ NA_character_,
      sub(" .*", "", ws_name) %in% sub(" .*", "", RESERVOIR_NAMES) ~ "Reservoir",
      TRUE ~ "Lake"),
    drains_south = {
      code <- if ("FWA_WATERSHED_CODE" %in% names(sb))
        dplyr::coalesce(as.character(sb$FWA_WATERSHED_CODE), "") else ""
      startsWith(code, "970") |
        grepl("Nooksack", dplyr::coalesce(ws_name, ""), ignore.case = TRUE)
    },
    cx, cy
  )

muni_out <- to_web(muni, coverage = TRUE) |>
  add_centroids() |>
  dplyr::transmute(muni_id, muni, rd, member_type, member_class, cx, cy)

aoi_web <- downstream_aoi |> sf::st_union() |> to_web()

fn_web <- to_web(fn_lands)

# ---- 5b. Orientation hydrography (lakes, rivers, streams) --------------------
# cache WFS pulls; reuse the cache if a later pull fails
fwa_cache <- file.path(paths()$processed, "99_hydrography_cache")
dir.create(fwa_cache, showWarnings = FALSE, recursive = TRUE)

fwa_pull <- function(layer, label, ...) {
  cache_f <- file.path(fwa_cache, paste0(layer, ".gpkg"))
  got <- tryCatch({
    if (!requireNamespace("bcdata", quietly = TRUE)) stop("bcdata not installed")
    suppressPackageStartupMessages(library(bcdata))
    dom <- sf::st_as_sfc(sf::st_bbox(sb)) # PROJECT_CRS domain bbox
    bcdc_query_geodata(paste0("WHSE_BASEMAPPING.", layer)) |>
      filter(INTERSECTS(dom), ...) |> collect()
  }, error = function(e) {
    message("  ! ", label, " pull failed: ", conditionMessage(e)); NULL
  })

  if (!is.null(got) && nrow(got) > 0) {
    try(sf::st_write(got, cache_f, delete_dsn = TRUE, quiet = TRUE), silent = TRUE)
    return(got)
  }
  if (file.exists(cache_f)) {
    message("  · ", label, ": reusing cached pull (", basename(cache_f), ")")
    return(sf::st_read(cache_f, quiet = TRUE))
  }
  message("  ! ", label, " unavailable and never cached, so the layer is empty")
  NULL
}

empty_layer <- function(...) sf::st_sf(..., geometry = sf::st_sfc(crs = "EPSG:4326"))

COARSE_M <- 150

lakes_web <- empty_layer()
lk <- fwa_pull("FWA_LAKES_POLY", "lakes")
if (!is.null(lk) && nrow(lk) > 0)
  lakes_web <- lk[lk$AREA_HA > 150, ] |>  # major lakes only (>1.5 km²)
    sf::st_transform(PROJECT_CRS) |> dplyr::summarise() |> to_web()

# punch lakes out of sub-basins, one unit at a time (st_difference drops emptied rows)
if (nrow(lakes_web) > 0) {
  lake_u <- sf::st_union(lakes_web)
  g <- sf::st_geometry(sb_out)
  for (i in which(lengths(sf::st_intersects(sb_out, lake_u)) > 0)) {
    cut <- sf::st_difference(g[i], lake_u)
    if (length(cut)) g[i] <- cut
  }
  sf::st_geometry(sb_out) <- g
}

rivers_web <- empty_layer()
rv <- fwa_pull("FWA_RIVERS_POLY", "rivers")
if (!is.null(rv) && nrow(rv) > 0)
  rivers_web <- rv[rv$AREA_HA > 30, ] |>  # prominent rivers only
    sf::st_transform(PROJECT_CRS) |> dplyr::summarise() |>
    sf::st_simplify(dTolerance = COARSE_M, preserveTopology = TRUE) |> to_web()

streams_web <- empty_layer(order = integer(0))
sn <- fwa_pull("FWA_STREAM_NETWORKS_SP", "streams", STREAM_ORDER >= 4)

us_sn <- NULL
if (any(is_us) && requireNamespace("nhdplusTools", quietly = TRUE)) {
  us_sn <- tryCatch({
    bb <- sf::st_bbox(sf::st_transform(sb[is_us, ], 4326))
    bb[1:2] <- floor(bb[1:2] * 1e4) / 1e4
    bb[3:4] <- ceiling(bb[3:4] * 1e4) / 1e4
    fl <- nhdplusTools::get_nhdplus(AOI = sf::st_as_sfc(bb), realization = "flowline")
    if (!inherits(fl, "sf") || !nrow(fl)) stop("NHDPlus returned no flowlines")
    oc <- grep("streamorde", names(fl), ignore.case = TRUE, value = TRUE)[1]
    fl$order <- if (!is.na(oc)) as.integer(fl[[oc]]) else 4L
    fl <- sf::st_zm(fl[fl$order >= 2, c("order")]) |> sf::st_transform(PROJECT_CRS)
    fl[lengths(sf::st_intersects(fl, sf::st_union(sb[is_us, ]))) > 0, ]
  }, error = function(e) {
    message("  ! US flowline pull skipped: ", conditionMessage(e)); NULL
  })
  if (!is.null(us_sn))
    message("  · US flowlines for the cascade: ", nrow(us_sn), " segment(s)")
}

if (!is.null(sn) && nrow(sn) > 0) {
  clean_lines <- function(x) {
    x <- sf::st_zm(x, drop = TRUE, what = "ZM")
    x <- sf::st_make_valid(x)
    x <- x[!sf::st_is_empty(sf::st_geometry(x)), ]
    x <- sf::st_cast(x, "MULTILINESTRING", warn = FALSE)
    x <- x[, intersect(c("order", "src"), names(x))]
    names(x)[names(x) == attr(x, "sf_column")] <- "geometry"
    sf::st_geometry(x) <- "geometry"
    x
  }
  sn_all <- sn |>
    sf::st_transform(PROJECT_CRS) |>
    dplyr::mutate(order = pmin(as.integer(STREAM_ORDER), 7L)) |>
    dplyr::mutate(src = "fwa") |>
    dplyr::select(order, src) |>
    clean_lines()
  if (!is.null(us_sn) && nrow(us_sn) > 0)
    sn_all <- rbind(sn_all, clean_lines(
      dplyr::mutate(us_sn, order = pmin(order, 7L), src = "nhd")))

  streams_web <- sn_all |>
    sf::st_simplify(dTolerance = COARSE_M, preserveTopology = TRUE) |>
    to_web()
  parts <- sf::st_cast(streams_web, "LINESTRING", warn = FALSE)
  keep <- as.numeric(sf::st_length(sf::st_transform(parts, PROJECT_CRS))) > 5
  parts <- parts[keep, ]

  # orient parts downstream: DEM where the drop is decisive, else source digitisation
  dem_f <- file.path(paths()$processed, "05_dem.tif")
  if (file.exists(dem_f) && nrow(parts) > 0) {
    dem <- terra::rast(dem_f)
    ends <- function(which_end) {
      cs <- lapply(sf::st_geometry(parts), function(g)
        g[if (which_end == 1) 1 else nrow(g), 1:2, drop = FALSE])
      m <- do.call(rbind, cs)
      terra::extract(dem, terra::project(
        terra::vect(m, crs = "EPSG:4326", type = "points"), terra::crs(dem)),
        ID = FALSE)[, 1]
    }
    z1 <- ends(1); z2 <- ends(2)
    dz <- z1 - z2
    flip <- is.finite(dz) & dz < 0

    DZ_TRUST_M <- 5
    decisive <- is.finite(dz) & abs(dz) >= DZ_TRUST_M
    by_source <- parts$src == "fwa"
    flip <- ifelse(decisive, dz < 0, by_source)

    idx <- which(flip)
    if (length(idx) > 0) {
      g <- sf::st_geometry(parts)
      g[idx] <- lapply(g[idx], function(l) sf::st_linestring(l[nrow(l):1, ]))
      sf::st_geometry(parts) <- g
    }
    message("  \u00b7 oriented ", nrow(parts), " stream part(s) downstream: ",
            length(idx), " reversed (", sum(decisive),
            " settled by the DEM, ", sum(!decisive),
            " by the source's own digitisation)")
  }

  STITCH_M <- 1000
  if (any(parts$src == "nhd") && any(parts$src == "fwa")) {
    us_i <- which(parts$src == "nhd"); bc_i <- which(parts$src == "fwa")
    gm <- sf::st_transform(sf::st_geometry(parts), PROJECT_CRS)
    pt <- function(idx, which_end) sf::st_sfc(lapply(idx, function(i) {
      l <- gm[[i]]; sf::st_point(l[if (which_end == 1) 1 else nrow(l), 1:2])
    }), crs = PROJECT_CRS)
    dn <- pt(us_i, 2)   # US downstream ends
    up <- pt(bc_i, 1)   # BC upstream ends
    nn <- sf::st_nearest_feature(dn, up)
    gap <- as.numeric(sf::st_distance(dn, up[nn], by_element = TRUE))
    ok <- which(is.finite(gap) & gap <= STITCH_M)
    if (length(ok) > 0) {
      g <- sf::st_geometry(parts)
      for (k in ok) {
        i <- us_i[k]; j <- bc_i[nn[k]]
        tgt <- unclass(g[[j]])[1, , drop = FALSE]
        li <- unclass(g[[i]])
        if (!identical(as.numeric(li[nrow(li), ]), as.numeric(tgt)))
          g[[i]] <- sf::st_linestring(rbind(li, tgt))
      }
      sf::st_geometry(parts) <- g

      off <- round(stats::median(
        parts$order[bc_i[nn[ok]]] - parts$order[us_i[ok]], na.rm = TRUE))
      if (is.finite(off) && off > 0) {
        parts$order[us_i] <- pmin(parts$order[us_i] + off, 7L)
        message("  \u00b7 raised US stream order by ", off,
                " to match the FWA scale across the border")
      }
    }
    message("  \u00b7 stitched ", length(ok), " of ", length(us_i),
            " US channel(s) onto the BC network at the border (median gap ",
            if (length(ok)) round(stats::median(gap[ok])) else NA, " m)")
  }

  streams_web <- parts |>
    dplyr::group_by(order) |>
    dplyr::summarise(do_union = FALSE, .groups = "drop")
}

# ---- 5c. Per-scenario retained runoff ----------------------------------------
scen_files <- list.files(paths()$processed,
                         pattern = "^12_realised_benefit_.*\\.gpkg$",
                         full.names = TRUE)
scen_ids <- sub("^12_realised_benefit_(.*)\\.gpkg$", "\\1", basename(scen_files))
scen_label_map <- c(wettest_month = "Wettest Month (800\u00a0m grid)",
                    ar2021        = "Atmospheric River Event (Nov 2021)")
pretty_scen <- function(s) if (s %in% names(scen_label_map))
  unname(scen_label_map[[s]]) else tools::toTitleCase(gsub("_", " ", s))

scen_data <- lapply(stats::setNames(scen_files, scen_ids), function(f) {
  d <- sf::st_read(f, quiet = TRUE) |> sf::st_drop_geometry()
  idx <- match(sb_out$id, idc(d$HYBAS_ID))

  retain <- round(dplyr::coalesce(d$prr_total_mm_km2, 0) * 1000)[idx]
  retain[is.na(retain)] <- 0
  depth <- dplyr::coalesce(d$prr_per_nat_mm_km2, 0)[idx]
  depth[is.na(depth)] <- 0
  gain <- round(ifelse(dplyr::coalesce(d$q_baseline_mm_km2, 0) > 0,
                       100 * dplyr::coalesce(d$prr_total_mm_km2, 0) /
                         d$q_baseline_mm_km2, 0))[idx]
  gain[is.na(gain)] <- 0

  blank <- sum(retain <= 0)
  if (blank > 0)
    message("  ! scenario '", tools::file_path_sans_ext(basename(f)),
            "': ", blank, " watershed(s) with no retention")

  list(m3 = stats::setNames(as.list(retain), sb_out$id),
       mm = stats::setNames(as.list(round(depth, 1)), sb_out$id),
       pct = stats::setNames(as.list(gain), sb_out$id))
})

retain_json <- jsonlite::toJSON(lapply(scen_data, `[[`, "m3"), auto_unbox = TRUE)
retain_mm_json <- jsonlite::toJSON(lapply(scen_data, `[[`, "mm"), auto_unbox = TRUE)
retain_pct_json <- jsonlite::toJSON(lapply(scen_data, `[[`, "pct"), auto_unbox = TRUE)

scen_order <- c(SCENARIO, setdiff(scen_ids, SCENARIO))
scen_order <- scen_order[scen_order %in% scen_ids]

# ---- 5d. Downstream designation overlays (display only) ----------------------
read_overlay <- function(id, class_field = NULL, clip = FALSE) {
  empty <- empty_layer(cls = character(0))
  tryCatch({
    p <- data_path(id)
    if (!file.exists(p)) stop("file not found: ", p)
    x <- sf::st_read(p, quiet = TRUE) |>
      sf::st_transform(PROJECT_CRS) |>
      sf::st_make_valid()
    x$cls <- if (!is.null(class_field)) {
      k <- names(x)[tolower(names(x)) == tolower(class_field)]
      if (length(k)) as.character(x[[k[1]]]) else NA_character_
    } else NA_character_
    if (clip) {
      x <- suppressWarnings(sf::st_intersection(x, sf::st_union(downstream_aoi)))
      x <- sf::st_make_valid(x)
    }
    x <- x[!sf::st_is_empty(sf::st_geometry(x)), ]
    if (!nrow(x)) stop("no features inside the downstream AOI")
    x <- x[, "cls"] |> to_web()
    message("  · overlay '", id, "': ", nrow(x), " feature(s)")
    x
  }, error = function(e) {
    message("  ! overlay '", id, "' skipped: ", conditionMessage(e)); empty
  })
}
alr_web  <- read_overlay("alr_polygons",  clip = TRUE)
rnin_p_web <- read_overlay("rnin_patch",    class_field = "STATUS")
rnin_c_web <- read_overlay("rnin_corridor", class_field = "Status")

# ---- 5e. Exposure by dissemination area --------------------------------------
das_web <- empty_layer(da_id = character(0), muni_id = integer(0))
da_demand_json <- "{}"
da_f <- file.path(paths()$processed, "10_demand_da.gpkg")
if (!file.exists(da_f)) {
  message("  ! DA exposure split skipped: ", basename(da_f),
          " not found (rerun 10_demand.R)")
} else {
  da <- sf::st_read(da_f, quiet = TRUE) |>
    sf::st_transform(PROJECT_CRS) |>
    sf::st_make_valid()
  da$da_id <- as.character(da$POP_UID)

  fp_r <- terra::rast(file.path(paths()$processed, "07_floodplain.tif"))
  fp_r <- terra::ifel(is.na(fp_r), 0, fp_r)
  fp_px <- (terra::xres(fp_r) * terra::yres(fp_r)) * 1e-6
  da$flood_km2 <- fp_px * exactextractr::exact_extract(fp_r, da, "sum",
                                                       progress = FALSE)
  da$pop_exposed <- da$COUNT_TOTAL *
    pmin(1, da$flood_km2 / (as.numeric(sf::st_area(da)) / 1e6))

  da <- suppressWarnings(sf::st_intersection(da, rd_union)) |> sf::st_make_valid()
  da <- da[!sf::st_is_empty(sf::st_geometry(da)), ]
  da$da_km2 <- as.numeric(sf::st_area(da)) / 1e6

  ov <- suppressWarnings(sf::st_intersection(da[, "da_id"], muni[, "muni_id"]))
  ov$a <- as.numeric(sf::st_area(ov))
  host <- ov |>
    sf::st_drop_geometry() |>
    dplyr::group_by(da_id) |>
    dplyr::slice_max(a, n = 1, with_ties = FALSE) |>
    dplyr::ungroup()
  da$muni_id <- host$muni_id[match(da$da_id, host$da_id)]
  da <- da[!is.na(da$muni_id), ]

  das_web <- to_web(da[, c("da_id", "muni_id")], coverage = TRUE)
  da_demand_json <- jsonlite::toJSON(
    stats::setNames(
      lapply(seq_len(nrow(da)), function(i) list(
        built = round(da$ben_built[i], 4),
        crops = round(da$ben_crops_w[i], 4),
        fac   = round(da$ben_facilities[i], 4),
        pop   = round(da$COUNT_TOTAL[i]),
        pop_x = round(da$pop_exposed[i]),
        km2   = round(da$da_km2[i], 3),
        muni  = as.integer(da$muni_id[i]))),
      da$da_id),
    auto_unbox = TRUE)
  message("  · DA exposure split: ", nrow(da),
          " exposed dissemination areas across ",
          dplyr::n_distinct(da$muni_id), " communities")

  da_by_muni <- da |>
    sf::st_drop_geometry() |>
    dplyr::group_by(muni_id) |>
    dplyr::summarise(fac = sum(ben_facilities), pop = sum(pop_exposed),
                     .groups = "drop")
  muni_dem <- muni_dem |>
    dplyr::full_join(da_by_muni, by = "muni_id")
}

muni_dem <- muni_dem |>
  dplyr::right_join(dplyr::tibble(muni_id = muni$muni_id), by = "muni_id") |>
  dplyr::mutate(dplyr::across(-muni_id, ~ dplyr::coalesce(.x, 0)))
muni_demand_json <- jsonlite::toJSON(
  stats::setNames(
    lapply(seq_len(nrow(muni_dem)), function(i) list(
      built = round(muni_dem$built[i], 3),
      crops = round(muni_dem$crops[i], 3),
      fac   = round(if (is.null(muni_dem$fac)) 0 else muni_dem$fac[i], 3),
      pop   = round(if (is.null(muni_dem$pop)) 0 else muni_dem$pop[i]))),
    muni_dem$muni_id),
  auto_unbox = TRUE)
message("  · member exposure: ", sum(muni_dem$built + muni_dem$crops > 0),
        " of ", nrow(muni_dem), " members carry exposure, ",
        round(sum(muni_dem$built + muni_dem$crops), 1), " km² in total")

# ---- 6. Emit data.js ---------------------------------------------------------

unwrap_collections <- function(gj) {
  plain <- vapply(gj$features,
                  function(f) if (is.null(f$geometry)) "" else f$geometry$type,
                  character(1))
  want <- if (any(plain %in% c("Polygon", "MultiPolygon")))
    c("Polygon", "MultiPolygon") else c("LineString", "MultiLineString")

  gj$features <- Filter(Negate(is.null), lapply(gj$features, function(f) {
    if (is.null(f$geometry) || !identical(f$geometry$type, "GeometryCollection"))
      return(f)
    keep <- Filter(function(gm) gm$type %in% want, f$geometry$geometries)
    if (!length(keep)) return(NULL)
    f$geometry <- keep[[1]]
    f
  }))
  gj
}

write_geojson <- function(x) {
  tf <- tempfile(fileext = ".geojson")
  on.exit(unlink(tf), add = TRUE)
  sf::st_write(x, tf, driver = "GeoJSON", quiet = TRUE,
               delete_dsn = TRUE,
               layer_options = c("RFC7946=YES", "COORDINATE_PRECISION=5"))
  txt <- paste(readLines(tf, warn = FALSE), collapse = "\n")

  if (grepl("GeometryCollection", txt, fixed = TRUE)) {
    gj <- unwrap_collections(jsonlite::fromJSON(txt, simplifyVector = FALSE))
    txt <- as.character(jsonlite::toJSON(gj, auto_unbox = TRUE, digits = NA,
                                         null = "null"))
  }
  txt
}

graph_json <- jsonlite::toJSON(
  stats::setNames(
    lapply(sb$id, function(v) list(contributes = as.integer(contributes[[v]]))),
    sb$id),
  auto_unbox = FALSE)

muni_lookup <- jsonlite::toJSON(
  stats::setNames(as.list(muni$muni), as.character(muni$muni_id)),
  auto_unbox = TRUE)

muni_basins <- jsonlite::toJSON(
  stats::setNames(
    lapply(muni$muni_id,
           function(m) as.character(basin_muni$id[basin_muni$muni_id == m])),
    muni$muni_id),
  auto_unbox = FALSE)

unit_label <- if ("FWA_WATERSHED_CODE" %in% names(sb))
  "BC Freshwater Atlas assessment watersheds" else "HydroBASINS L12 sub-basins"

meta <- jsonlite::toJSON(list(
  scenario = SCENARIO,
  scenarios = lapply(scen_order, function(s) list(id = s, label = pretty_scen(s))),
  scenario_labels = stats::setNames(as.list(vapply(scen_ids, pretty_scen, "")), scen_ids),
  unit = unit_label,
  n_subbasin = nrow(sb_out),
  n_muni = nrow(muni_out),
  n_muni_mvrd = sum(muni$rd == "MVRD"),
  n_muni_fvrd = sum(muni$rd == "FVRD"),
  n_da = nrow(das_web),
  generated = as.character(Sys.Date())
), auto_unbox = TRUE)

data_js <- paste0(
  "// Auto-generated by R/99_interactive_map.R. Do not edit by hand.\n",
  "const META = ", meta, ";\n\n",
  "const SUBBASINS = ", write_geojson(sb_out), ";\n\n",
  "const MUNIS = ", write_geojson(muni_out), ";\n\n",
  "const DAS = ", write_geojson(das_web), ";\n\n",
  "const AOI = ", write_geojson(sf::st_sf(geometry = aoi_web)), ";\n\n",
  "const FIRSTNATIONS = ", write_geojson(fn_web), ";\n\n",
  "const LAKES = ", write_geojson(lakes_web), ";\n\n",
  "const RIVERS = ", write_geojson(rivers_web), ";\n\n",
  "const STREAMS = ", write_geojson(streams_web), ";\n\n",
  "const ALR = ", write_geojson(alr_web), ";\n\n",
  "const RNIN_PATCH = ", write_geojson(rnin_p_web), ";\n\n",
  "const RNIN_CORRIDOR = ", write_geojson(rnin_c_web), ";\n\n",
  "const GRAPH = ", graph_json, ";\n\n",
  "const MUNI_NAMES = ", muni_lookup, ";\n\n",
  "const MUNI_BASINS = ", muni_basins, ";\n\n",
  "const MUNI_DEMAND = ", muni_demand_json, ";\n\n",
  "const DEMAND = ", demand_json, ";\n\n",
  "const DA_DEMAND = ", da_demand_json, ";\n\n",
  "const RETAIN = ", retain_json, ";\n\n",
  "const RETAIN_MM = ", retain_mm_json, ";\n\n",
  "const RETAIN_PCT = ", retain_pct_json, ";\n\n",
  "const NATTYPES = ", jsonlite::toJSON(nat_comp, auto_unbox = TRUE), ";\n"
)

writeLines(data_js, file.path(out_dir, "data.js"))

message("✓ 99_interactive_map.R: wrote ", file.path(out_dir, "data.js"),
        " (", nrow(sb_out), " sub-basins, ", nrow(muni_out), " municipalities)")
