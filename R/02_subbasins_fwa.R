# ==============================================================================
# 02_subbasins_fwa.R — sub-basin polygons + downstream topology + upstream AOI
# ------------------------------------------------------------------------------
# Replaces the national-study HydroBASINS L12 with the BC Freshwater Atlas
# Assessment Watersheds data product. FWA has no NEXT_DOWN field, so we derive 
# the downstream topology from the FWA watershed-code system.
#
# Inputs:
#   WHSE_BASEMAPPING.FWA_ASSESSMENT_WATERSHEDS_POLY (DataBC WFS, via bcdata)
#   fwa_stream_networks (FWA_STREAM_NETWORKS_SP.gdb)
#   hydrobasins_l12 — (hydrobasins_l12.shp)
#   data/processed/01_aoi_downstream.gpkg
#
# Outputs (data/processed/):
#   02_subbasins.gpkg — sub-basin polygons + tags
#   02_topology.csv — focal_id, ds_id, reach_km
#   01_aoi_upstream.gpkg
# ==============================================================================

source(here::here("R", "00_setup.R"))

# Lower Mainland demand polygon (MVRD ∪ FVRD)
downstream_aoi <- read_aoi("downstream")

# ---- 1. FWA assessment watersheds --------------------------------------------
suppressPackageStartupMessages(library(bcdata))

# expand the AOI by the flow-distance cap + edge buffer
fwa_buf <- sf::st_buffer(sf::st_union(downstream_aoi),
                         MAX_FLOW_DIST_KM * 1000 + EDGE_BUFFER_M) |> sf::st_geometry()

# pull every watershed intersecting that buffer from the DataBC WFS
message("pulling FWA assessment watersheds from DataBC WFS…")

fwa <- bcdc_query_geodata("WHSE_BASEMAPPING.FWA_ASSESSMENT_WATERSHEDS_POLY") |>
  filter(INTERSECTS(fwa_buf)) |>
  collect()

# drop the CAD-annotation blob column, re-project to BC Albers
fwa <- fwa[, setdiff(names(fwa), "SE_ANNO_CAD_DATA")] |>
  sf::st_transform(PROJECT_CRS)

# ---- 2. Parse watershed code -------------------------------------------------
# FWA_WATERSHED_CODE "100-190442-244975-...-000000" = a base drainage (100) plus
# six-digit junction segments ordered DOWNSTREAM -> UPSTREAM. LOCAL_WATERSHED_CODE
# adds one segment giving position along the watershed's own stream (smaller =
# further downstream)

# split an FWA code into its non-zero tributary segments (drop 000000 segments)
none_zero <- function(code) {
  s <- strsplit(code, "-", fixed = TRUE)[[1]]
  s[s != "000000"]
}

# keep only the attributes needed for the upstream/downstream calculations
base <- fwa |>
  sf::st_drop_geometry() |>
  dplyr::transmute(id = WATERSHED_FEATURE_ID,
                   magnitude = WATERSHED_MAGNITUDE,
                   area_km2 = AREA_HA / 100,
                   fwa = FWA_WATERSHED_CODE, 
                   loc = LOCAL_WATERSHED_CODE)

# parse each watershed's FWA code into its tributary hierarchy
segments <- lapply(base$fwa, none_zero)

# depth = number of levels in the hierarchy (higher = farther upstream)
depth <- lengths(segments)

# local measure = watershed's position along its own tributary
measure <- as.numeric(
  mapply(\(l, d)
         strsplit(l, "-", fixed = TRUE)[[1]][d + 1],
         base$loc, depth))

# no local subdivision -> 0
measure[is.na(measure)] <- 0

# lookup: watershed id -> row number, for faster indexing
row_of_id <- setNames(seq_len(nrow(base)), base$id)

# ---- 3. Pairwise downstream test ---------------------------------------------
# is_downstream(a, b): TRUE if b lies downstream of a. The shared FWA address
# tests whether they sit on the same drainage path; the local measure tests
# whether b is lower than a on that path.

is_downstream <- function(a, b) {

  # convert watershed IDs to row positions
  row_a <- row_of_id[[as.character(a)]]
  row_b <- row_of_id[[as.character(b)]]

  # FWA tributary addresses for a and b
  addr_a <- segments[[row_a]]
  addr_b <- segments[[row_b]]

  # downstream watershed cannot be deeper in the tributary hierarchy than a
  if (length(addr_b) > length(addr_a)) return(FALSE)

  # b must share the beginning of a's FWA address
  if (!all(addr_b == addr_a[seq_along(addr_b)])) return(FALSE)

  # two independent coastal catchments may both share the ocean as their outlet,
  # but they never flow into one another (e.g. Sunshine Coast vs North Shore)
  if (base$magnitude[row_a] == 0 && base$magnitude[row_b] == 0) return(FALSE)

  # flow should only accumulate downstream, so b should be > or = upstream network magnitude as a
  if (base$magnitude[row_b] < base$magnitude[row_a]) return(FALSE)

  if (length(addr_b) == length(addr_a)) {
    # case 1: same tributary, where b is downstream if its local position is lower than a
    measure[row_b] < measure[row_a]

  } else {
    # case 2: b is on a downstream trunk/main stem.
    # the next segment in a's FWA address marks where a's tributary joins b's stream
    confluence <- as.numeric(addr_a[length(addr_b) + 1])

    # b is downstream only if it is at or below that confluence.
    measure[row_b] <= confluence
  }
}

# ---- 4. NEXT_DOWN = adjacent downstream neighbour carrying the most flow ------
# find shared boundaries between watersheds
touch <- sf::st_touches(fwa)

# expand to one row per (a, neighbour b) pair
pairs <- data.frame(a = rep(base$id, lengths(touch)),
                    b = base$id[unlist(touch)])

# check if a flow into b (or if b is downstream)
pairs$flows <- mapply(is_downstream, pairs$a, pairs$b)

# identify b's magnitude to later pick the trunk among neighbours
pairs$b_magnitude <- base$magnitude[row_of_id[as.character(pairs$b)]]

# choose one downstream neighbour per watershed: the highest-magnitude = the trunk
next_down <- pairs |>
  dplyr::filter(flows) |> # keep only pairs where a -> b
  dplyr::group_by(a) |>
  dplyr::slice_max(b_magnitude, n = 1, with_ties = FALSE) |> # choose the downstream neighbour most likely to be the main trunk
  dplyr::ungroup() |>
  dplyr::select(id = a, next_down = b)

# join next down to base
base <- base |>
  dplyr::left_join(next_down, by = "id") |>
  dplyr::mutate(next_down = dplyr::coalesce(next_down, 0)) # replaces missing values with 0

# ---- 4b. Reach length = measured stream length from FWA stream networks -------
# preferred: clipped length of the watershed's own-code mainstem
# fallback:  dominant BLUE_LINE_KEY inside the watershed
reach_from_streams <- function(fwa_sf, code_by_id, aoi_buf) { # polygon, lookup table, study area
  gdb <- data_path("fwa_stream_networks")
  message("reading AOI streams from FWA_STREAM_NETWORKS gdb…")

  # read only stream segments intersecting the buffered AOI
  bbox_wkt <- sf::st_as_text(sf::st_as_sfc(sf::st_bbox(aoi_buf)))
  layers <- sf::st_layers(gdb)
  spatial_layers <- layers$name[!is.na(layers$geomtype)]

  # keep only fields needed for stream-length assignment
  keep_cols <- c("FWA_WATERSHED_CODE", "BLUE_LINE_KEY", "LENGTH_METRE")
  stream_layers <- list()

  for (layer in spatial_layers) {
    layer_streams <- sf::st_read(gdb, layer = layer, quiet = TRUE, wkt_filter = bbox_wkt)

    if (nrow(layer_streams) > 0) {
      stream_layers[[layer]] <- sf::st_zm(layer_streams[, intersect(keep_cols, names(layer_streams))])
    }
  }

  streams <- do.call(rbind, stream_layers)

  # clip stream segments to watershed boundaries
  message("clipping ", nrow(streams), " stream segments to watershed boundaries…")

  unit_ids <- fwa_sf["WATERSHED_FEATURE_ID"]
  names(unit_ids)[1] <- "id"
  clip <- suppressWarnings(sf::st_intersection(streams, unit_ids))

  # keep only line geometries produced by the intersection
  clip <- clip[
    sf::st_geometry_type(clip) %in% c("LINESTRING", "MULTILINESTRING"),
  ]

  # convert clipped stream pieces to a plain table and measure their length
  clip_tbl <- sf::st_drop_geometry(clip)
  clip_tbl$clip_km <- as.numeric(sf::st_length(clip)) / 1000

  # primary reach length:
  # sum clipped stream pieces whose FWA code matches the watershed's own code
  own_reach <- clip_tbl[
    clip_tbl$FWA_WATERSHED_CODE == code_by_id[as.character(clip_tbl$id)],
  ] |>
    dplyr::group_by(id) |>
    dplyr::summarise(km = sum(clip_km), .groups = "drop")

  # fallback reach length:
  # use the single BLUE_LINE_KEY with the most clipped length inside the watershed
  dom_blue <- clip_tbl |>
    dplyr::group_by(id, BLUE_LINE_KEY) |>
    dplyr::summarise(l = sum(clip_km), .groups = "drop") |>
    dplyr::group_by(id) |>
    dplyr::summarise(km = max(l), .groups = "drop")

  # combine primary and fallback reach lengths for every watershed.
  reach_tbl <- data.frame(id = fwa_sf$WATERSHED_FEATURE_ID) |>
    dplyr::left_join(own_reach, by = "id") |>
    dplyr::left_join(dom_blue,  by = "id", suffix = c(".own", ".dom"))

  message(
    "reach: ", sum(!is.na(reach_tbl$km.own)), " own-mainstem, ",
    sum(is.na(reach_tbl$km.own) & !is.na(reach_tbl$km.dom)), " dominant-blue-line fallback, ",
    sum(is.na(reach_tbl$km.own) & is.na(reach_tbl$km.dom)), " streamless (reach 0)"
  )

  reach <- dplyr::coalesce(reach_tbl$km.own, reach_tbl$km.dom, 0)
  stats::setNames(reach, as.character(reach_tbl$id))
}

# measure reach length and join to base
measured_reach <- reach_from_streams(
  fwa,
  stats::setNames(base$fwa, base$id),
  fwa_buf
)

base$reach_km <- as.numeric(measured_reach[as.character(base$id)])

message(
  "measured reach length for ",
  sum(base$reach_km > 0), "/", nrow(base),
  " watersheds"
)

# ---- 4c. Splice in US (HydroBASINS) providers draining into the BC AOI --------
# FWA coverage stops at the 49th parallel, but several Whatcom County (US)
# watersheds drain north into the BC demand area and so provide real upstream
# attenuation. FWA is limited to BC, so we take those units from HydroBASINS L12,
# keep only the ones that (a) sit south of the border and (b) actually route into
# the demand AOI, and append them to base/fwa as extra providers.

# columns carried on the spatial layer (must match the FWA layer's schema below)
us_cols <- c("WATERSHED_FEATURE_ID", "WATERSHED_ORDER", "WATERSHED_MAGNITUDE",
             "FWA_WATERSHED_CODE", "GNIS_NAME_1", "AREA_HA")

# read HydroBASINS, keep only units intersecting the buffered AOI
hydrobasins <- sf::st_read(data_path("hydrobasins_l12"), quiet = TRUE) |>
  sf::st_transform(PROJECT_CRS)
hydrobasins <- hydrobasins[as.logical(
  sf::st_intersects(hydrobasins, sf::st_union(fwa_buf), sparse = FALSE)[, 1]), ]

# flag US units by centroid latitude (HydroBASINS spans both sides of the border)
centroid_lat <- sf::st_coordinates(
  sf::st_transform(suppressWarnings(sf::st_centroid(sf::st_geometry(hydrobasins))), 4326))[, 2]
is_us <- centroid_lat < 49.0

# build the HydroBASINS flow graph (HYBAS_ID -> NEXT_DOWN), then filter for units
# that can reach a unit that touches the demand AOI
hydrobasins_graph <- igraph::graph_from_data_frame(
  hydrobasins |>
    sf::st_drop_geometry() |>
    dplyr::filter(NEXT_DOWN != 0, NEXT_DOWN %in% HYBAS_ID) |>
    dplyr::transmute(from = as.character(HYBAS_ID), to = as.character(NEXT_DOWN)),
  vertices = data.frame(name = as.character(hydrobasins$HYBAS_ID)))

# units touching the demand AOI = the targets we test reachability against
aoi_unit_ids <- hydrobasins$HYBAS_ID[as.logical(
  sf::st_intersects(hydrobasins, downstream_aoi, sparse = FALSE)[, 1])]

# reaches_aoi = finite downstream distance to at least one AOI-touching unit
reachable <- is.finite(apply(igraph::distances(
  hydrobasins_graph, to = igraph::V(hydrobasins_graph)[name %in% as.character(aoi_unit_ids)],
  mode = "out"), 1, min))
hydrobasins$reaches_aoi <- reachable[match(as.character(hydrobasins$HYBAS_ID),
                                           igraph::V(hydrobasins_graph)$name)]

# the providers we keep: south of the border AND draining into the demand AOI
us_ids   <- hydrobasins$HYBAS_ID[is_us & hydrobasins$reaches_aoi %in% TRUE]
us_units <- hydrobasins[hydrobasins$HYBAS_ID %in% us_ids, ]

if (nrow(us_units) > 0) {
  # Cross-border NEXT_DOWN: keep a US unit if its downstream neighbour is a kept
  # US unit. Otherwise it flows into BC, so re-point it at the
  # FWA watershed that overlaps its (BC-side) HydroBASINS neighbour the most
  us_next_down <- us_units$NEXT_DOWN
  for (i in which(!(us_next_down %in% us_ids))) {
    target <- hydrobasins[hydrobasins$HYBAS_ID == us_units$NEXT_DOWN[i], ]
    overlap <- if (nrow(target)) suppressWarnings(
      sf::st_intersection(fwa[, "WATERSHED_FEATURE_ID"], sf::st_geometry(target))) else target[0, ]
    us_next_down[i] <- if (nrow(overlap))
      overlap$WATERSHED_FEATURE_ID[which.max(as.numeric(sf::st_area(overlap)))] else 0
  }
  # append providers to the attribute table; no FWA streams in the US, so reach
  # length is estimated from area via Hack's law (reach_km ≈ 1.4 * area_km2^0.6)
  base <- dplyr::bind_rows(base, data.frame(
    id = us_units$HYBAS_ID, magnitude = NA_real_, area_km2 = us_units$SUB_AREA,
    fwa = NA_character_, loc = NA_character_, next_down = us_next_down,
    reach_km = 1.4 * us_units$SUB_AREA^0.6))
  # and to the spatial layer, mapped into the FWA schema (US-only fields left NA)
  us_polys <- us_units |>
    dplyr::transmute(
    WATERSHED_FEATURE_ID = HYBAS_ID, WATERSHED_ORDER = NA_real_,
    WATERSHED_MAGNITUDE = NA_real_, FWA_WATERSHED_CODE = NA_character_,
    GNIS_NAME_1 = NA_character_, AREA_HA = SUB_AREA * 100)
  fwa <- rbind(fwa[, us_cols], us_polys)
  message("spliced ", nrow(us_units), " US (HydroBASINS) providers; flow-distance trim decides which survive")
}

# ---- 5. Trim to ≤ MAX_FLOW_DIST_KM upstream of the downstream AOI -------------
# copy NEXT_DOWN and reach length onto the spatial layer
fwa$NEXT_DOWN <- base$next_down[match(fwa$WATERSHED_FEATURE_ID, base$id)]
fwa$reach_km <- base$reach_km[match(fwa$WATERSHED_FEATURE_ID, base$id)]

# all watershed ids = graph vertices
all_ids <- base$id

# build directed flow edges: each watershed points to its immediate downstream neighbour, with edge weight equal to the reach length (km)
edges <- base |>
  dplyr::filter(next_down != 0, next_down %in% all_ids) |>
  dplyr::transmute(from = as.character(id),
                   to = as.character(next_down),
                   weight = reach_km)

# flow network
flow_graph <- igraph::graph_from_data_frame(edges, vertices = data.frame(name = as.character(all_ids)))

# watersheds touching the downstream AOI
outlet_idx <- as.logical(sf::st_intersects(fwa, downstream_aoi, sparse = FALSE)[, 1])
outlet_ids <- fwa$WATERSHED_FEATURE_ID[outlet_idx]

# network distance from every watershed down to each outlet
dist_matrix <- igraph::distances(
  flow_graph,
  to = igraph::V(flow_graph)[name %in% as.character(outlet_ids)],
  mode = "out")

# shortest downstream distance to any outlet
flow_dist <- apply(dist_matrix, 1, min, na.rm = TRUE)

# align back to base rows
base$flow_dist_km <- flow_dist[match(base$id, as.numeric(igraph::V(flow_graph)$name))]

# keep only watersheds within the flow-distance cap
keep <- is.finite(base$flow_dist_km) &
  base$flow_dist_km <= MAX_FLOW_DIST_KM

# trim the spatial layer
fwa <- fwa[match(base$id[keep], fwa$WATERSHED_FEATURE_ID), ]

# trim the attribute table to match
base <- base[keep, ]

# NEXT_DOWN now pointing at a dropped watershed
dangling <- !base$next_down %in% c(0, base$id)

if (any(dangling)) {
  message("rewrote ", sum(dangling), " out-of-domain NEXT_DOWN references to 0")
  # treat them as domain outlets
  base$next_down[dangling] <- 0
}

# refresh NEXT_DOWN after the rewrite
fwa$NEXT_DOWN <- base$next_down[match(fwa$WATERSHED_FEATURE_ID, base$id)]

# ---- 5b. Drop drainage-isolated islands --------------------------------------
# watershed sharing no boundary with any other = island/sliver
island <- lengths(sf::st_touches(fwa)) == 0

if (any(island)) {
  message("dropped ", sum(island), " drainage-isolated island/sliver watershed(s)")
  # drop islands as they do not attenuate any downstream areas
  fwa <- fwa[!island, ]
  # trim the attribute table to match (nothing routes into an island, so no dangling refs)
  base <- base[base$id %in% fwa$WATERSHED_FEATURE_ID, ]
}

# ---- 5c. Lake barriers: sever routing through large lakes --------------------
# Restores the lake-as-sink rule lost in the HydroBASINS -> FWA migration (FWA has
# no LAKE/ENDO field). Big lakes/reservoirs absorb the inflow flood pulse in their
# own storage, so upstream attenuation gives ~no marginal benefit downstream: flag
# each qualifying lake's unit(s) and cut their downstream edge.
base$is_lake_barrier <- FALSE
lakes <- tryCatch({
  lakes_all <- bcdc_query_geodata("WHSE_BASEMAPPING.FWA_LAKES_POLY") |>
    filter(INTERSECTS(fwa_buf)) |> collect()
  # qualifying barriers: large lakes (surface area, ha -> km^2) OR named reservoirs
  lakes_all[(lakes_all$AREA_HA / 100 >= LAKE_BARRIER_KM2) | (lakes_all$GNIS_NAME_1 %in% RESERVOIR_NAMES), ]
}, error = function(e) { message("  ! FWA lakes pull skipped: ", conditionMessage(e)); NULL })

if (!is.null(lakes) && nrow(lakes) > 0) {
  # dissolve by lake name so multi-polygon lakes count as one waterbody
  lakes <- lakes |> 
    sf::st_transform(PROJECT_CRS) |>
    dplyr::group_by(GNIS_NAME_1) |> 
    dplyr::summarise(.groups = "drop") |>
    sf::st_make_valid()
  
  unit_polys <- fwa["WATERSHED_FEATURE_ID"]; 
  names(unit_polys)[1] <- "id"
  lake_overlap <- suppressWarnings(sf::st_intersection(unit_polys, lakes))
  lake_overlap$area <- as.numeric(sf::st_area(lake_overlap))
  overlap_tbl <- sf::st_drop_geometry(lake_overlap)
  unit_area <- stats::setNames(as.numeric(sf::st_area(unit_polys)), unit_polys$id)
  
  # A unit is a barrier if a qualifying lake (a) covers >= LAKE_COVER_FRAC of it
  # (lake-dominated unit, handles multi-unit lakes) OR (b) is the unit a NAMED
  # lake overlaps most (its outlet/core unit — catches small reservoirs sitting
  # inside a larger mixed unit, e.g. Coquitlam/Alouette/Capilano).
  
  covered_ids <- unique(overlap_tbl$id[(overlap_tbl$area / unit_area[as.character(overlap_tbl$id)]) >= LAKE_COVER_FRAC])
  core_ids <- overlap_tbl |>
    dplyr::filter(!is.na(GNIS_NAME_1)) |>
    dplyr::group_by(GNIS_NAME_1) |>
    dplyr::slice_max(area, n = 1, with_ties = FALSE) |> 
    dplyr::pull(id)
  
  base$is_lake_barrier <- base$id %in% unique(c(covered_ids, core_ids))
}

if (any(base$is_lake_barrier)) {
  message("lake barriers: severed ", sum(base$is_lake_barrier),
          " large-lake unit(s) from downstream routing (",
          paste(stats::na.omit(fwa$GNIS_NAME_1[match(base$id[base$is_lake_barrier],
                fwa$WATERSHED_FEATURE_ID)]), collapse = ", "), ")")
  base$next_down[base$is_lake_barrier] <- 0
}
# carry the flag + refreshed NEXT_DOWN onto the spatial layer
fwa$is_lake_barrier <- base$is_lake_barrier[match(fwa$WATERSHED_FEATURE_ID, base$id)]
fwa$NEXT_DOWN <- base$next_down[match(fwa$WATERSHED_FEATURE_ID, base$id)]

# ---- 6. AOI flags + sink flag + admin district -------------------------------
# use the FWA watershed ID as the analysis-unit ID used in subsequent scripts
fwa$HYBAS_ID <- fwa$WATERSHED_FEATURE_ID

# watershed area in km^2
fwa$SUB_AREA <- fwa$AREA_HA / 100

# attach downstream flow distance calculated from the routing graph
fwa$flow_dist_km <- base$flow_dist_km[match(fwa$WATERSHED_FEATURE_ID, base$id)]

# flag watersheds whose centroid falls inside the downstream/demand AOI
fwa$in_downstream_aoi <- as.logical(
  sf::st_intersects(sf::st_centroid(fwa), downstream_aoi, sparse = FALSE)[, 1]
)

# flag domain outlets (watersheds with no downstream neighbour inside the model)
fwa$is_sink <- fwa$NEXT_DOWN == 0

# tag each watershed by regional district
if (requireNamespace("bcmaps", quietly = TRUE)) {

  # load regional district polygons and project them to the working CRS
  districts <- bcmaps::regional_districts() |>
    sf::st_transform(PROJECT_CRS)

  # default label for watersheds outside the three focal districts
  fwa$admin_district <- "Other"

  # label watersheds whose centroid falls inside one of the focal districts
  for (district_name in c(
    "Metro Vancouver Regional District",
    "Fraser Valley Regional District",
    "Squamish-Lillooet Regional District"
  )) {
    poly <- dplyr::filter(districts, ADMIN_AREA_NAME == district_name)

    if (nrow(poly) == 1) {
      in_district <- as.logical(
        sf::st_intersects(sf::st_centroid(fwa), poly, sparse = FALSE)[, 1]
      )

      # store a shorter district name by removing " Regional District"
      fwa$admin_district[in_district] <- gsub(" Regional District$", "", district_name)
    }
  }
}

# spliced US providers (no FWA watershed code) are tagged distinctly
fwa$admin_district[is.na(fwa$FWA_WATERSHED_CODE)] <- "Whatcom (US)"

# ---- 7. Persist sub-basin layer ----------------------------------------------
# columns needed by downstream scripts, plus selected original FWA identifiers
keep_cols <- c("HYBAS_ID", "NEXT_DOWN", "SUB_AREA", "flow_dist_km", "reach_km",
               "in_downstream_aoi", "is_sink", "is_lake_barrier", "admin_district",
               "WATERSHED_FEATURE_ID", "WATERSHED_ORDER", "WATERSHED_MAGNITUDE",
               "FWA_WATERSHED_CODE", "GNIS_NAME_1")

# write the final routed sub-basin polygon layer
sf::st_write(fwa[, intersect(keep_cols, names(fwa))],
             file.path(paths()$processed, "02_subbasins.gpkg"),
             delete_dsn = TRUE, quiet = TRUE)

# ---- 8. Topology CSV ---------------------------------------------------------
# flow topology as a table
topo <- base |>
  dplyr::transmute(focal_id = id, ds_id = next_down, reach_km = reach_km)

# output processed data
readr::write_csv(topo, file.path(paths()$processed, "02_topology.csv"))

# ---- 9. Refined upstream AOI -------------------------------------------------
# dissolve all kept watersheds, then buffer the edge
refined <- fwa |>
  sf::st_union() |>
  # mitre join preserves sharp watershed corners, avoiding the rounded edge cut-offs produced by the default buffer
  sf::st_buffer(
    EDGE_BUFFER_M,
    joinStyle = "MITRE",
    mitreLimit = 2) |>

  # convert back to an sf object
  sf::st_sf(role = "upstream", crs = PROJECT_CRS)

# overwrite the placeholder upstream AOI from 01_aoi.R
sf::st_write(refined, file.path(paths()$processed, "01_aoi_upstream.gpkg"),
             delete_dsn = TRUE, quiet = TRUE)

# ---- QA preview --------------------------------------------------------------
# Three single-purpose panels:
#   1. Topology   — the DERIVED NEXT_DOWN edges over a plain basemap; edges should
#                   form a dendritic tree converging on the terminal units, which
#                   are classified by colour (demand outlet / domain edge / lake)
#   2. Flow dist  — network distance down to the admin AOI (drives the step-5 trim)
#   3. Reach km   — measured main-channel length per unit (step 4b)
# ---- small helpers shared by the panels --------------------------------------
# per-feature colours + value range for a continuous fill
qa_contcol <- function(v, pal) {
  rng <- range(v, finite = TRUE)
  idx <- cut(v, seq(rng[1], rng[2], length.out = length(pal) + 1),
             include.lowest = TRUE, labels = FALSE)
  list(col = pal[idx], rng = rng)
}
# vertical continuous colour bar in its own layout cell
qa_draw_bar <- function(rng, pal, title) {
  graphics::par(mar = c(3, 0.4, 4, 3.2))
  ys <- seq(rng[1], rng[2], length.out = length(pal) + 1)
  graphics::image(x = c(0, 1), y = ys, z = matrix(seq_along(pal), nrow = 1),
                  col = pal, axes = FALSE, xlab = "", ylab = "")
  graphics::axis(4, las = 1, cex.axis = 0.95, lwd = 0, lwd.ticks = 1)
  graphics::box(); graphics::mtext(title, side = 3, line = 0.6, cex = 0.9, font = 2)
}
# draw a filled map panel with the AOI outline on top
qa_fill_panel <- function(geom, fill, title) {
  graphics::par(mar = c(2, 1, 4, 1))
  plot(geom, col = fill, border = "grey85", lwd = 0.1, reset = FALSE, main = title)
  plot(sf::st_geometry(downstream_aoi), add = TRUE, border = "black", lwd = 1.8)
}

qa_png("02_subbasins_fwa.png", ncol = 3, panel_w = 1200, panel_h = 1150,
       plot_fn = function() {
  geom      <- sf::st_geometry(fwa)
  flow_pal  <- hcl.colors(100, "Greens", rev = TRUE)
  reach_pal <- hcl.colors(100, "Blues",  rev = TRUE)

  # demand footprint = floodplain ∩ downstream AOI:
  # shows the demand is the Hope->Salish Sea valley, not the mountainous admin AOI
  fp_path <- file.path(paths()$processed, "07_floodplain.tif")
  fp <- NULL
  if (file.exists(fp_path)) {
    fp <- terra::aggregate(terra::mask(terra::rast(fp_path),
                                       terra::vect(downstream_aoi)),
                           8, fun = "max", na.rm = TRUE)
    fp[fp == 0] <- NA
  }

  graphics::layout(matrix(1:5, nrow = 1), widths = c(1, 1, 0.20, 1, 0.20))

  # ---- panel 1 — derived drainage topology + terminal-unit classification ----
  # Classify every terminal unit (NEXT_DOWN == 0) honestly rather than calling
  # them all "outlets":
  #   demand outlet — reaches the floodplain demand area (a true outlet)
  #   domain edge   — drains OUT of the study area (its real downstream was trimmed)
  #   lake barrier  — large lake; routing deliberately severed in step 5c
  is_lake  <- fwa$is_lake_barrier %in% TRUE
  sink_oth <- fwa$is_sink & !is_lake
  fp_touch <- rep(FALSE, nrow(fwa))
  if (!is.null(fp) && any(sink_oth)) {
    ex <- terra::extract(fp, terra::vect(fwa[sink_oth, ]),
                         fun = function(x) sum(x > 0, na.rm = TRUE))
    fp_touch[which(sink_oth)] <- !is.na(ex[, 2]) & ex[, 2] > 0
  }
  is_demand <- sink_oth & fp_touch
  is_edge   <- sink_oth & !fp_touch

  graphics::par(mar = c(2, 1, 4, 1))
  plot(geom, col = "grey95", border = "grey80", lwd = 0.1, reset = FALSE,
       main = "Drainage topology — NEXT_DOWN tree + terminal units")
  if (!is.null(fp))
    terra::plot(fp, add = TRUE, col = "#2171b5", alpha = 0.6, legend = FALSE)

  # NEXT_DOWN edges: each unit's centroid -> its downstream neighbour's centroid
  centroids <- sf::st_coordinates(suppressWarnings(sf::st_centroid(geom)))
  rownames(centroids) <- as.character(fwa$HYBAS_ID)
  edge_df <- data.frame(id = fwa$HYBAS_ID, nd = fwa$NEXT_DOWN)
  edge_df <- edge_df[edge_df$nd != 0 & edge_df$nd %in% fwa$HYBAS_ID, ]
  x0 <- centroids[as.character(edge_df$id), "X"]; y0 <- centroids[as.character(edge_df$id), "Y"]
  x1 <- centroids[as.character(edge_df$nd), "X"]; y1 <- centroids[as.character(edge_df$nd), "Y"]
  graphics::segments(x0, y0, x1, y1, col = "white",  lwd = 1.1)  # halo
  graphics::segments(x0, y0, x1, y1, col = "grey10", lwd = 0.45) # core

  # outline each terminal class in its own colour
  if (any(is_demand))
    plot(sf::st_geometry(fwa[is_demand, ]), add = TRUE, border = "red", lwd = 1.4)
  if (any(is_edge))
    plot(sf::st_geometry(fwa[is_edge, ]), add = TRUE, border = "darkorange", lwd = 1.0)
  if (any(is_lake))
    plot(sf::st_geometry(fwa[is_lake, ]), add = TRUE, border = "#08519c", lwd = 2.2)
  plot(sf::st_geometry(downstream_aoi), add = TRUE, border = "black", lwd = 1.8)

  leg_txt <- c(sprintf("%d demand outlet(s) → floodplain", sum(is_demand)),
               sprintf("%d domain-edge (drains out of study area)", sum(is_edge)),
               sprintf("%d lake barrier(s) (routing severed)", sum(is_lake)),
               "NEXT_DOWN edge", "downstream AOI (admin)")
  leg_fill <- rep(NA, 5)
  leg_col  <- c("red", "darkorange", "#08519c", "grey10", "black")
  leg_lty  <- rep(1, 5); leg_lwd <- c(1.4, 1.0, 2.2, 1, 1.8)
  if (!is.null(fp)) {  # floodplain swatch only when 07 floodplain exists
    leg_txt  <- c("floodplain = demand area", leg_txt)
    leg_fill <- c("#2171b5", leg_fill); leg_col <- c(NA, leg_col)
    leg_lty  <- c(NA, leg_lty);         leg_lwd <- c(NA, leg_lwd)
  }
  graphics::legend("topleft", legend = leg_txt, fill = leg_fill, border = NA,
                   col = leg_col, lty = leg_lty, lwd = leg_lwd, bty = "n", cex = 0.9)

  # ---- panel 2 — flow distance to the admin AOI (drives the step-5 trim) ------
  # NOTE: network distance to the *admin* downstream AOI; the benefit model (11)
  # decays by reach to the floodplain, not by this distance
  flow_col <- qa_contcol(fwa$flow_dist_km, flow_pal)
  qa_fill_panel(geom, flow_col$col, "Flow distance to AOI (network km)")
  qa_draw_bar(flow_col$rng, flow_pal, "flow dist. (km)")

  # ---- panel 3 — measured main-channel reach length (step 4b) -----------------
  reach_col <- qa_contcol(fwa$reach_km, reach_pal)
  qa_fill_panel(geom, reach_col$col, "Measured reach length (FWA streams)")
  qa_draw_bar(reach_col$rng, reach_pal, "reach (km)")
})

message("✓ 02_subbasins_fwa.R — wrote FWA sub-basins + topology + upstream AOI (",
        nrow(fwa), " assessment watersheds; ", sum(fwa$is_sink), " terminal sinks incl. ",
        sum(fwa$is_lake_barrier %in% TRUE), " lake barriers)")