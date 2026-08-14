# ==============================================================================
# 99_interactive_map.R: Web data for the upstream-ecosystems interactive map
# ------------------------------------------------------------------------------
# Builds the data behind output/interactive_map/index.html, a page for exploring
# which upstream watersheds attenuate flooding for which downstream communities,
# and which of them contribute to several communities at once.
#
# The model already supplies, per sub-basin, realised benefit, provision
# intensity (PRR), demand, and the NEXT_DOWN flow graph. This script adds the
# municipal dimension:
#   1. link each downstream demand basin to the jurisdictions it overlaps;
#   2. walk the flow graph so every upstream basin knows the full set of
#      jurisdictions its outflow eventually reaches;
#   3. emit GeoJSON (WGS84, simplified) inlined into data.js
#
# Inputs (data/processed/):
#   02_subbasins.gpkg, 10_demand_subbasin.gpkg, 10_demand_da.gpkg,
#   11_tda_subbasin.gpkg
#   12_realised_benefit_<scenario>.gpkg
#   01_aoi_downstream.gpkg, 01_mvrd.gpkg, 01_fvrd.gpkg
#   bcmaps::census_subdivision() - members of each regional district
#   data_sources.csv - ids
#   optional display overlays - alr_polygons / rnin_patch / rnin_corridor
#                               (a missing file drops the layer)
#
# Outputs (output/interactive_map/data.js):
#   META, SUBBASINS, MUNIS, DAS, AOI, FIRSTNATIONS, LAKES, RIVERS, STREAMS, ALR,
#   RNIN_PATCH, RNIN_CORRIDOR, GRAPH, MUNI_NAMES, MUNI_BASINS, MUNI_DEMAND,
#   DEMAND, DA_DEMAND, RETAIN, RETAIN_MM, RETAIN_PCT, NATTYPES
#
# MUNI_BASINS and MUNI_DEMAND answer two different questions and must not be
# used for each other's: MUNI_BASINS is where a member's flood risk sits, so the
# map can route drainage lines to it; MUNI_DEMAND is how much exposure the
# member holds. Summing basin demand over MUNI_BASINS credits a member with
# every basin it touches and overstates the region by 3.9x.
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

# ---- 1a. Territory claimed by both the US and BC sub-basin sets --------------
# 02 selects US (WBD HUC12) units on their unclipped shapes, so a transboundary
# unit can still hold the Canadian half of its watershed, ground the FWA maps as
# an assessment watershed of its own. Clip the overlap off, or the border strip
# is drawn twice and its runoff credited twice.
if (any(is_us) && any(!is_us)) {
  s2_prev <- sf::sf_use_s2(FALSE)
  bc_union <- sf::st_union(sf::st_geometry(sb[!is_us, ]))
  us_geom <- suppressWarnings(
    sf::st_difference(sf::st_geometry(sb[is_us, ]), bc_union))
  before <- as.numeric(sf::st_area(sf::st_geometry(sb[is_us, ]))) / 1e6
  after <- as.numeric(sf::st_area(us_geom)) / 1e6
  # sliver tolerance: unioning polygons digitised by two national surveys leaves
  # hairline overlaps along the join, 0.4 km² across the domain. Clip those
  # silently; a total this large instead means 02's border fix is missing.
  SLIVER_KM2 <- 5
  trimmed <- (before - after) > 0.01
  if (any(trimmed)) {
    g <- sf::st_geometry(sb)
    g[which(is_us)] <- sf::st_make_valid(us_geom)
    sf::st_geometry(sb) <- g
    sb$SUB_AREA[is_us] <- after
    if (sum(before - after) > SLIVER_KM2) warning(sprintf(paste0(
      "%s US basin(s) overlapped the FWA by %.1f km2: 02_subbasins.gpkg predates ",
      "the border fix in 02_subbasins_fwa.R. Their land cover and retention were ",
      "extracted over BC ground that belongs to the neighbouring FWA watershed ",
      "and are double-counted. Rerun 02 -> 12 before trusting this map."),
      sum(trimmed), sum(before - after)), immediate. = TRUE)
  }
  sf::sf_use_s2(s2_prev)
}

# ---- 1c. Natural cover-type composition per sub-basin ------------------------
# which KIND of natural land does the attenuating: forest / shrub / grass /
# wetland / alpine km² per sub-basin, so the map can tell a community what its
# upstream contributing land is made of
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
# census subdivisions are the legal partition of an RD into municipalities,
# electoral areas and reserve/treaty lands, and they tile it with no gaps.
# bcmaps::municipalities() keeps only incorporated municipalities, which left
# several members off the map.
csd_all <- bcmaps::census_subdivision() |>
  sf::st_transform(PROJECT_CRS) |>
  sf::st_make_valid()

# Short type label + class from the CSD type. Treaty first nations (e.g.
# Tsawwassen) carry a blank TYPE_DESC but a treaty-land TYPE_CODE, so those fall
# back to the code.
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

# Tag each member by regional district. A CSD nests within one regional district,
# so a point-on-surface test assigns each cleanly and drops CSDs outside the two
# study RDs. The RD polygons come from 01_aoi.R, so the split matches
# downstream_aoi (= MVRD ∪ FVRD).
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

# First Nation lands are dissolved into their host members below, so the
# community polygons stay solid. The map still shows them as their own
# toggleable layer.
fn_lands <- muni[muni$member_class == "First Nation", c("name", "member_type", "rd")]
message("  · First Nation lands layer: ", nrow(fn_lands), " CSD(s) (",
        sum(fn_lands$rd == "MVRD"), " MVRD / ", sum(fn_lands$rd == "FVRD"), " FVRD)")

# dissolve reserves into the member containing them: CSDs carve reserves OUT of
# the surrounding member, and Leaflet strokes the resulting holes as inner rings
# even with the reserve layer hidden. Each reserve's land and floodplain demand
# goes to its host, so reserves are not their own communities in this view.
is_res <- muni$member_type == "Reserve"
if (any(is_res) && any(!is_res)) {
  host_rows <- which(!is_res); res_rows <- which(is_res)
  # a carved-out reserve shares its whole hole boundary with its host, so the
  # nearest host feature (distance 0) is the one it is embedded in
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
# the CSD and regional-district products are generalized differently, so their
# edges disagree by thin slivers (~80 km² of RD with no CSD cover, ~110 km² of
# CSD overhang plus a marine wedge). The RD boundary is the model's demand AOI,
# so it wins: clip the communities to it, then fold each uncovered sliver into
# the community it shares the most boundary with.

# of the members a sliver touches, the one sharing the longest boundary with it
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
# a sub-basin serves a member where it overlaps that member AND has demand of
# its own (built-up or cropland in the floodplain).
#
# Barriers are excluded even when they qualify: an FWA unit like "Capilano
# River" holds both the reservoir and urban demand at its mouth, so seeding the
# graph there would let the upstream walk start AT the barrier and credit the
# lake's whole catchment. Section 4 excludes barriers the same way.
demand_basins <- sb$id[!is.na(sb$own_demand) & sb$own_demand > 0 &
                         !sb$is_lake_barrier]
dem$id <- idc(dem$HYBAS_ID)

# Split each demand basin along the community boundaries and measure the
# exposure inside every piece, from the per-pixel raster 10 wrote. A member is
# served by a basin when it holds exposure THERE, not when the two polygons
# merely overlap: an area test credited a member clipping the corner of a
# valley-floor watershed with the whole of that watershed's exposure, which
# summed across members to 3.9x the region's total and put electoral areas at
# the top of the ranking.
inter <- suppressWarnings(
  sf::st_intersection(
    sb[sb$id %in% demand_basins, c("id")],
    muni[, c("muni_id", "muni")]
  )
) |> sf::st_make_valid()
# members that only touch along an edge intersect as lines or collections, and
# the zonal step takes polygons only
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

# Exposure per member, apportioned from each basin's own components by the
# share of that basin's exposure lying inside the member. Apportioning the
# components rather than re-measuring them keeps the ledger exact: each
# component sums across members to the same total it sums to across basins, so
# the region total is the sum of its parts and nothing is counted twice.
# The apportioning weight comes from the unweighted pixel raster, so a future
# non-uniform crop weighting would shift the split slightly but never the total.
muni_dem <- basin_muni |>
  dplyr::group_by(id) |>
  dplyr::mutate(share = exp_km2 / sum(exp_km2)) |>
  dplyr::ungroup() |>
  dplyr::left_join(dplyr::select(dem, id, ben_built, ben_crops_w), by = "id") |>
  dplyr::group_by(muni_id) |>
  dplyr::summarise(built = sum(share * ben_built),
                   crops = sum(share * ben_crops_w),
                   .groups = "drop")

# demand components per demand basin, km² (crops use the weighted column the
# model routes). Only demand basins ship; the rest are zero on all three. No
# basin carries facility demand without land demand, so the sliders change how
# much each basin counts, never which ones seed the network.
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

# Lake and reservoir units are terminal: 02 severs their routing, so they pass
# no benefit downstream and are not credited as served demand either. Crediting
# the held-back runoff to the lake itself is a deferred refinement.
barrier_ids <- sb$id[sb$is_lake_barrier]

for (v in sb$id) {
  desc <- setdiff(names(igraph::subcomponent(g, v, mode = "out")), v)
  up_list[[v]] <- setdiff(names(igraph::subcomponent(g, v, mode = "in")), v)
  # members this watershed contributes to = members of every demand basin
  # downstream of it (itself included), barriers removed
  reach <- setdiff(c(v, desc), barrier_ids)
  pm <- unique(unlist(direct_munis[intersect(reach, names(direct_munis))]))
  contributes[[v]] <- sort(as.integer(pm))
}

sb$n_contrib <- vapply(contributes, length, integer(1))
sb$n_up <- vapply(up_list, length, integer(1))

# Split that count by regional district, so the map can compare where an upstream
# watershed's benefit lands and flag the cross-jurisdiction watersheds that reach
# communities in both.
rd_of_id <- stats::setNames(muni$rd, as.character(muni$muni_id))
sb$n_contrib_mvrd <- vapply(contributes,
  function(p) sum(rd_of_id[as.character(p)] == "MVRD"), integer(1))
sb$n_contrib_fvrd <- vapply(contributes,
  function(p) sum(rd_of_id[as.character(p)] == "FVRD"), integer(1))

# ---- 5. Build web geometries (WGS84, simplified) -----------------------------

# reproject + simplify for the browser. coverage = TRUE for polygon layers that
# TILE the domain (sub-basins, members): simplified independently, each side of
# a shared border moves differently and the basemap shows through as white
# seams. rmapshaper::ms_simplify simplifies each shared arc once instead;
# without it, fall back to st_simplify and accept the seams.
to_web <- function(x, coverage = FALSE) {
  x <- sf::st_transform(x, PROJECT_CRS)
  if (any(sf::st_geometry_type(x) == "GEOMETRYCOLLECTION"))
    x <- sf::st_collection_extract(x, "POLYGON", warn = FALSE)
  if (coverage && requireNamespace("rmapshaper", quietly = TRUE)) {
    # keep = fraction of vertices retained; keep_shapes drops no small unit
    x <- rmapshaper::ms_simplify(x, keep = 0.12, keep_shapes = TRUE,
                                 explode = FALSE)
  } else {
    x <- sf::st_simplify(x, dTolerance = SIMPLIFY_M, preserveTopology = TRUE)
  }
  x |>
    sf::st_transform("EPSG:4326") |>
    sf::st_make_valid()
}

# Label anchor for each polygon, as plain lon/lat columns the web app can read
# without touching the geometry.
add_centroids <- function(x) {
  pt <- suppressWarnings(sf::st_point_on_surface(sf::st_geometry(x)))
  co <- sf::st_coordinates(pt)
  x$cx <- round(co[, 1], 5)
  x$cy <- round(co[, 2], 5)
  x
}

# FWA names most assessment watersheds after a stream or lake (Alouette River,
# Big Silver Creek, …), which labels the otherwise anonymous upstream areas on
# hover. HydroBASINS has no equivalent.
sb$ws_name <- if ("GNIS_NAME_1" %in% names(sb)) as.character(sb$GNIS_NAME_1) else NA_character_

# ---- 4b. Close the pinholes in the watershed coverage ------------------------
# FWA assessment watersheds cover land only, so the atlas carves out water
# surfaces and leaves about 160 gaps, most under a tenth of a km². At the scale
# this map is read at they render as white hairlines through the shading.
#
# Fold every gap below HOLE_FILL_KM2 into the watershed it shares the most
# boundary with, as 2b does for communities. Display only: this runs after every
# model quantity is extracted, so no figure moves. Larger holes are genuine
# water bodies and stay, and 5b punches the mapped lakes out so they read as
# water.
HOLE_FILL_KM2 <- 1

local({
  s2_prev <- sf::sf_use_s2(FALSE)
  on.exit(sf::sf_use_s2(s2_prev), add = TRUE)

  # interior rings of the dissolved coverage = the gaps between watersheds
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

# only the fields the web app reads. Percentiles, TDA and own-demand stay in the
# GPKGs. The map covers retention volume and drainage connectivity.
sb_out <- sb_web |>
  dplyr::transmute(
    id, next_id,
    name = ws_name,
    natural_km2 = round(natural_km2, 1),
    retain_m3 = round(dplyr::coalesce(prr_total_mm_km2, 0) * 1000),
    # retention DEPTH: mm held back per km² of natural land, the same quantity
    # 99_fig2 maps per pixel. Does not scale with watershed size, so it reads as
    # a property of the ecosystem rather than of the polygon.
    retain_mm = round(dplyr::coalesce(prr_per_nat_mm_km2, 0), 1),
    # retained runoff as a share of the runoff this watershed produces anyway:
    # "losing this land's natural cover would increase storm runoff from here by
    # N%". Both terms come from 09 via 12, where PRR is already
    # Q_bare - Q_natural, so this restates the model's own comparison.
    runoff_gain_pct = round(
      ifelse(dplyr::coalesce(q_baseline_mm_km2, 0) > 0,
             100 * dplyr::coalesce(prr_total_mm_km2, 0) / q_baseline_mm_km2, 0)),
    reach_km = round(dplyr::coalesce(reach_km, 0), 2),
    n_contrib, n_contrib_mvrd, n_contrib_fvrd, n_up,
    # ship the barrier flag + type so the map can label these units rather than
    # leave them as unexplained faint patches
    barrier = is_lake_barrier,
    # a barrier unit is named for its stream ("Capilano River") while the
    # reservoir list is lake-named ("Capilano Lake"), so match on the base name
    barrier_type = dplyr::case_when(
      !is_lake_barrier ~ NA_character_,
      sub(" .*", "", ws_name) %in% sub(" .*", "", RESERVOIR_NAMES) ~ "Reservoir",
      TRUE ~ "Lake"),
    cx, cy
  )

muni_out <- to_web(muni, coverage = TRUE) |>
  add_centroids() |>
  dplyr::transmute(muni_id, muni, rd, member_type, member_class, cx, cy)

# downstream AOI outline for context
aoi_web <- downstream_aoi |> sf::st_union() |> to_web()

# First Nation lands, a display overlay outside the watershed coverage
fn_web <- to_web(fn_lands)

# ---- 5b. Orientation hydrography (lakes, rivers, streams) --------------------

# pull an FWA layer over the domain bbox. The three hydrography layers are
# decorative, so a failed pull is not fatal, but the DataBC WFS goes down
# intermittently and an empty layer silently drops the stream network from
# data.js. Cache each success and fall back to it, so only a first run is empty.
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

# these three are for orientation only, so they are simplified harder than the
# model layers, trading vertex count for payload size and render smoothness.
COARSE_M <- 150

# Major lakes, so their labels (Harrison, Pitt, Stave …) sit on a visible water
# body instead of floating among drainage lines.
lakes_web <- empty_layer()
lk <- fwa_pull("FWA_LAKES_POLY", "lakes")
if (!is.null(lk) && nrow(lk) > 0)
  lakes_web <- lk[lk$AREA_HA > 150, ] |>  # major lakes only (>1.5 km²)
    sf::st_transform(PROJECT_CRS) |> dplyr::summarise() |> to_web()

# Punch the lakes out of the sub-basins so a lake reads as water rather than
# shaded land. cx/cy are already computed, so labels and markers are unaffected.
#
# Cut basin by basin: st_difference drops a geometry it erases completely, and
# two units are entirely lake, so cutting the layer in one call would take those
# off the map and leave the rest misaligned.
if (nrow(lakes_web) > 0) {
  lake_u <- sf::st_union(lakes_web)
  g <- sf::st_geometry(sb_out)
  for (i in which(lengths(sf::st_intersects(sb_out, lake_u)) > 0)) {
    cut <- sf::st_difference(g[i], lake_u)
    if (length(cut)) g[i] <- cut
  }
  sf::st_geometry(sb_out) <- g
}

# The Fraser main stem and its major tributaries as real river polygons, so the
# trunk reads as an actual river beneath the drainage network and flow animation.
rivers_web <- empty_layer()
rv <- fwa_pull("FWA_RIVERS_POLY", "rivers")
if (!is.null(rv) && nrow(rv) > 0)
  rivers_web <- rv[rv$AREA_HA > 30, ] |>  # prominent rivers only
    sf::st_transform(PROJECT_CRS) |> dplyr::summarise() |>
    sf::st_simplify(dTolerance = COARSE_M, preserveTopology = TRUE) |> to_web()

# Real stream centrelines for "where the water flows", replacing a synthetic
# centroid-to-centroid skeleton that drew every confluence as a radial star.
# Order >= 4 keeps the major channels and drops the order-3 mesh, which roughly
# doubled payload for little orientation value. Dissolving per order lets the
# app width-code them.
#
# The FWA stops at 49°N, so on its own this layer ends dead on the border while
# the model routes US watersheds across it, and the drainage lines and cascade
# animation both ride on this layer. NHDPlus carries the same channels on the US
# side (02 already uses it for reach lengths), so pull the matching orders and
# merge them in before the per-order dissolve.
streams_web <- empty_layer(order = integer(0))
sn <- fwa_pull("FWA_STREAM_NETWORKS_SP", "streams", STREAM_ORDER >= 4)

us_sn <- NULL
if (any(is_us) && requireNamespace("nhdplusTools", quietly = TRUE)) {
  us_sn <- tryCatch({
    aoi_us <- sf::st_as_sfc(sf::st_bbox(sf::st_transform(sb[is_us, ], 4326)))
    fl <- nhdplusTools::get_nhdplus(AOI = aoi_us, realization = "flowline")
    oc <- grep("streamorde", names(fl), ignore.case = TRUE, value = TRUE)[1]
    fl$order <- if (!is.na(oc)) as.integer(fl[[oc]]) else 4L
    # small headwater systems, so the >= 4 cut that keeps the BC payload sane
    # would discard almost all of them. The US footprint is a sliver of the
    # domain, so take order >= 2 and let the channels reach the border.
    fl <- sf::st_zm(fl[fl$order >= 2, c("order")]) |> sf::st_transform(PROJECT_CRS)
    # keep only what actually lies in a US watershed: the NHD bbox overlaps BC,
    # and the FWA already owns everything north of the line
    fl[lengths(sf::st_intersects(fl, sf::st_union(sb[is_us, ]))) > 0, ]
  }, error = function(e) {
    message("  ! US flowline pull skipped: ", conditionMessage(e)); NULL
  })
  if (!is.null(us_sn))
    message("  · US flowlines for the cascade: ", nrow(us_sn), " segment(s)")
}

if (!is.null(sn) && nrow(sn) > 0) {
  # the two sources disagree on the details GEOS cares about (FWA carries Z, NHD
  # carries ZM, and the two mix LINESTRING with MULTILINESTRING), and the union
  # inside summarise() fails on the mixture. Normalise both to plain 2D
  # multilines first.
  clean_lines <- function(x) {
    x <- sf::st_zm(x, drop = TRUE, what = "ZM")
    x <- sf::st_make_valid(x)
    x <- x[!sf::st_is_empty(sf::st_geometry(x)), ]
    sf::st_cast(x, "MULTILINESTRING", warn = FALSE)
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

  # NO union here. summarise() without do_union = FALSE dissolves and re-nodes
  # the geometry, which throws away each feature's vertex order, and that order
  # is the FWA telling us which way the water goes. Simplify per feature instead
  # and keep them separate until they have been oriented.
  streams_web <- sn_all |>
    sf::st_simplify(dTolerance = COARSE_M, preserveTopology = TRUE) |>
    to_web()
  # Drop near-zero-length parts: 5-decimal coordinate rounding collapses them to
  # stray POINTs that L.geoJSON would draw as default teardrop markers. Explode,
  # drop the sub-5 m fragments, recombine per order (no union / noding).
  parts <- sf::st_cast(streams_web, "LINESTRING", warn = FALSE)
  keep <- as.numeric(sf::st_length(sf::st_transform(parts, PROJECT_CRS))) > 5
  parts <- parts[keep, ]

  # Store every part in downstream order, decided by the DEM. Sample elevation
  # at both ends and reverse the ones that run uphill. The app animates its
  # cascade from the start of a line to its end, so settling direction once here
  # replaces its old guess from sub-basin flow distances, which drew 15% of the
  # chains uphill by a median of 247 m.
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
    dz <- z1 - z2                        # positive = already downstream
    flip <- is.finite(dz) & dz < 0

    # On flat ground the gradient says nothing, which covers the whole lower
    # Fraser and is why the main stem used to animate upstream. There the source
    # settles it: FWA segments are digitised downstream-first (98.6% of the
    # 4,805 with an unambiguous >20 m drop), NHDPlus flowlines the other way.
    # The DEM overrides wherever it IS decisive (>5 m), catching the 1.4% of FWA
    # segments digitised against the convention.
    DZ_TRUST_M <- 5
    decisive <- is.finite(dz) & abs(dz) >= DZ_TRUST_M
    # want first vertex = upstream, so FWA (downstream-first) flips, NHD does not
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

  # ---- stitch the two stream networks across the 49th parallel ---------------
  # The model routes US watersheds straight into BC ones, but FWA and NHDPlus
  # are separate survey products whose lines share no vertices, and the web app
  # chains segments by exact shared endpoints, so a cascade the model says
  # continues stops dead on the border.
  #
  # Give them the vertex they are missing. Every US channel ends at its
  # downstream end and every BC channel starts at its upstream end, so where
  # such a pair sits within STITCH_M, extend the US line to end exactly on the
  # BC line's start. The coordinate is copied, so it survives the 5-decimal
  # rounding and the app chains straight through.
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

      # each product computes stream order on its own network density: the Sumas
      # arrives as a wide FWA order-5 channel and continues as an NHD order-2
      # thread, so an uncorrected join makes a real river peter out at the
      # boundary. Order only width-codes the line, so shift the US orders by the
      # offset the stitched pairs imply and the channel keeps its weight.
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

# ---- 5c. Per-scenario retained runoff (for the in-map scenario toggle) -------
# Ship every available scenario's retained-runoff volume, keyed by sub-basin id,
# so the browser can switch rainfall scenario and recolour without a rebuild.
# retain_m3 = PRR (mm·km²) × 1000.
scen_files <- list.files(paths()$processed,
                         pattern = "^12_realised_benefit_.*\\.gpkg$",
                         full.names = TRUE)
scen_ids <- sub("^12_realised_benefit_(.*)\\.gpkg$", "\\1", basename(scen_files))
# labels name what the reader is choosing between, a typical heavy-rain month
# against one real storm
scen_label_map <- c(wettest_month = "Wettest Month (800\u00a0m grid)",
                    ar2021        = "Atmospheric River Event (Nov 2021)")
pretty_scen <- function(s) if (s %in% names(scen_label_map))
  unname(scen_label_map[[s]]) else tools::toTitleCase(gsub("_", " ", s))

# read each basin's retention straight from that scenario's file. Every scenario
# covers the whole AOI, so a zero means no natural land left to retain runoff,
# not a gap. sb_out preserves sb's row order, so sb's columns index it directly.
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

  # a scenario that leaves whole watersheds empty is a coverage failure
  # upstream, so report it and leave the gap in place
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

# active scenario first, so it leads the toggle
scen_order <- c(SCENARIO, setdiff(scen_ids, SCENARIO))
scen_order <- scen_order[scen_order %in% scen_ids]

# ---- 5d. Downstream designation overlays (display only) ----------------------
# Context layers over the downstream corridor: the Agricultural Land Reserve and
# Metro Vancouver's Regional Natural Infrastructure Network (patches +
# corridors). These are not model inputs. No upstream script reads them, and
# they ship only so the map can show where they sit against the modelled
# contribution.
# Linking them to the model is future work.
#
# The ALR is province-wide (~46,000 km²) so it is clipped to the downstream AOI;
# RNIN is already inside it and is only reprojected. Same fallback as the
# hydrography pulls: an unreadable file leaves the layer empty and the web app
# omits its toggle.
read_overlay <- function(id, class_field = NULL, clip = FALSE) {
  empty <- empty_layer(cls = character(0))
  tryCatch({
    p <- data_path(id)
    if (!file.exists(p)) stop("file not found: ", p)
    x <- sf::st_read(p, quiet = TRUE) |>
      sf::st_transform(PROJECT_CRS) |>
      sf::st_make_valid()
    # the class label drives styling + legend; sources disagree on case (patch
    # uses STATUS, corridor uses Status), so match case-insensitively
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

# ---- 5e. Exposure below the community level (dissemination areas) ------------
# a finer view of the same demand rather than a finer model. Ranking, selection
# and routing stay per community. Shading a large community by its total
# exposure says nothing about where that exposure sits, and the Township of
# Langley's
# floodplain is all north.
#
# 10 splits the same demand pixels over census DAs (400-700 people) with the
# same three components, so the map's exposure weightings apply unchanged. Only
# exposed DAs ship, a few hundred of roughly four thousand, and the community
# fill still shows elsewhere. A DA straddling two communities is drawn under
# whichever holds most of it, and DA totals do not add to community totals,
# which count a shared floodplain once per community. A missing file drops the
# layer.
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
  da$da_id <- as.character(da$DAUID)

  # Exposed population, measured before the clip below so the count and the area
  # it is apportioned over describe the same polygon. 13_population.R splits a
  # DA's residents by the share of it lying in the FLOODPLAIN, so this does too;
  # apportioning by the built-up and cropland share instead would report a
  # different number (25,553 against 34,581) from the same census input.
  fp_r <- terra::rast(file.path(paths()$processed, "07_floodplain.tif"))
  fp_r <- terra::ifel(is.na(fp_r), 0, fp_r) # punched cells are dry, as in 10
  fp_px <- (terra::xres(fp_r) * terra::yres(fp_r)) * 1e-6
  da$flood_km2 <- fp_px * exactextractr::exact_extract(fp_r, da, "sum",
                                                       progress = FALSE)
  da$pop_exposed <- da$COUNT_TOTAL *
    pmin(1, da$flood_km2 / (as.numeric(sf::st_area(da)) / 1e6))

  # DAs tile the census land base, not the regional districts, so clip them to
  # the same boundary the community layer was snapped to in 2b. Demand is
  # already zeroed outside the downstream AOI by 10, so the clipped-off part
  # carries none and the exposure values stand.
  da <- suppressWarnings(sf::st_intersection(da, rd_union)) |> sf::st_make_valid()
  da <- da[!sf::st_is_empty(sf::st_geometry(da)), ]
  da$da_km2 <- as.numeric(sf::st_area(da)) / 1e6

  # attribute each DA to the community holding the most of it
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

  # Facilities and residents per member come from the DA table rather than the
  # basin apportionment above: a facility is a point sitting in one DA, and
  # residents are counted per DA, so both land in a member directly instead of
  # being spread over a watershed. A DA is counted whole, under the member
  # holding most of it.
  da_by_muni <- da |>
    sf::st_drop_geometry() |>
    dplyr::group_by(muni_id) |>
    dplyr::summarise(fac = sum(ben_facilities), pop = sum(pop_exposed),
                     .groups = "drop")
  muni_dem <- muni_dem |>
    dplyr::full_join(da_by_muni, by = "muni_id")
}

# every member ships a row, so the map never has to guess at a missing one
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

# COORDINATE_PRECISION rounding collapses thin slivers into lines or points and
# wraps the feature in a GeometryCollection, which L.geoJSON draws as stray
# teardrop markers. Keep only the member matching the layer's own dimension, and
# drop features where nothing survived (two RNIN patches did exactly that).
unwrap_collections <- function(gj) {
  # the layer's dimension, judged from the features that are not collections
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
  # via temp file rather than an sf-version-dependent string writer
  tf <- tempfile(fileext = ".geojson")
  on.exit(unlink(tf), add = TRUE)
  sf::st_write(x, tf, driver = "GeoJSON", quiet = TRUE,
               delete_dsn = TRUE,
               layer_options = c("RFC7946=YES", "COORDINATE_PRECISION=5"))
  txt <- paste(readLines(tf, warn = FALSE), collapse = "\n")

  # only pay the parse cost when a collection is actually present
  if (grepl("GeometryCollection", txt, fixed = TRUE)) {
    gj <- unwrap_collections(jsonlite::fromJSON(txt, simplifyVector = FALSE))
    txt <- as.character(jsonlite::toJSON(gj, auto_unbox = TRUE, digits = NA,
                                         null = "null"))
  }
  txt
}

# GRAPH[id].contributes = the members each sub-basin contributes to. The web app
# reads nothing else from the graph; the n_up / n_contrib counts ride on the
# SUBBASINS features.
graph_json <- jsonlite::toJSON(
  stats::setNames(
    lapply(sb$id, function(v) list(contributes = as.integer(contributes[[v]]))),
    sb$id),
  auto_unbox = FALSE)

muni_lookup <- jsonlite::toJSON(
  stats::setNames(as.list(muni$muni), as.character(muni$muni_id)),
  auto_unbox = TRUE)

# Terminal (demand) sub-basins per member: where that community's flood risk
# sits. The web app converges its drainage lines on these so the flow actually
# reaches the community.
muni_basins <- jsonlite::toJSON(
  stats::setNames(
    lapply(muni$muni_id,
           function(m) as.character(basin_muni$id[basin_muni$muni_id == m])),
    muni$muni_id),
  auto_unbox = FALSE)

# Name the analysis unit for the footer. Only the FWA layer carries
# FWA_WATERSHED_CODE, so this reads the layer instead of taking a flag.
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
