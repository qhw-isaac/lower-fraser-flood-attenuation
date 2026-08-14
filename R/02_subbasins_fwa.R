# ==============================================================================
# 02_subbasins_fwa.R: sub-basin polygons + downstream topology + upstream AOI
# ------------------------------------------------------------------------------
# Replaces the national-study HydroBASINS L12 with the BC Freshwater Atlas
# assessment watersheds. The FWA has no field naming the next watershed
# downstream, so this script works that out from the FWA watershed codes.
#
# Inputs:
#   WHSE_BASEMAPPING.FWA_ASSESSMENT_WATERSHEDS_POLY (DataBC WFS, via bcdata)
#   fwa_stream_networks (FWA_STREAM_NETWORKS_SP.gdb)
#   US transboundary providers: WBD HUC12 + NHDPlus flowlines (USGS web
#     services, via nhdplusTools), replacing the coarser HydroBASINS L12 splice
#   data/processed/01_aoi_downstream.gpkg
#
# Outputs (data/processed/):
#   02_subbasins.gpkg   sub-basin polygons + tags
#   02_topology.csv   focal_id, ds_id, reach_km
#   01_aoi_upstream.gpkg
# ==============================================================================

source(here::here("R", "00_setup.R"))

# the Lower Mainland demand polygon (MVRD ∪ FVRD)
downstream_aoi <- read_aoi("downstream")

# ---- 1. FWA assessment watersheds --------------------------------------------
# expand the AOI by the flow-distance cap + edge buffer
fwa_buf <- sf::st_buffer(sf::st_union(downstream_aoi),
                         MAX_FLOW_DIST_KM * 1000 + EDGE_BUFFER_M) |> sf::st_geometry()

# pull every watershed intersecting that buffer from the DataBC WFS
message("pulling FWA assessment watersheds from DataBC WFS…")

fwa <- bcdc_query_geodata("WHSE_BASEMAPPING.FWA_ASSESSMENT_WATERSHEDS_POLY") |>
  filter(INTERSECTS(fwa_buf)) |>
  collect()

# drop unneeded column, re-project to BC Albers
fwa <- fwa[, setdiff(names(fwa), "SE_ANNO_CAD_DATA")] |>
  sf::st_transform(PROJECT_CRS)

# Keep a copy of the survey's own shapes before 4c and 4d start cutting and
# merging them. Step 4e needs it to tell apart a fragment the FWA itself mapped
# from one this script created. 41 of these watersheds arrive in several pieces.
fwa_source_geom <- stats::setNames(sf::st_geometry(fwa),
                                   as.character(fwa$WATERSHED_FEATURE_ID))

# ---- 2. Parse watershed code -------------------------------------------------
# FWA_WATERSHED_CODE "100-190442-244975-...-000000" is a base drainage (100)
# plus six-digit junction segments ordered DOWNSTREAM -> UPSTREAM.
# LOCAL_WATERSHED_CODE adds one more segment giving position along the
# watershed's own stream, where smaller is further downstream.

# an FWA code's non-zero tributary segments
non_zero <- function(code) {
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
segments <- lapply(base$fwa, non_zero)

# depth = number of levels in the hierarchy (higher = farther upstream)
depth <- lengths(segments)

# local measure = watershed's position along its own tributary
measure <- as.numeric(
  mapply(\(l, d)
         strsplit(l, "-", fixed = TRUE)[[1]][d + 1],
         base$loc, depth))

# no local subdivision -> 0
measure[is.na(measure)] <- 0

# lookup from watershed id to row number, for faster indexing
row_of_id <- setNames(seq_len(nrow(base)), base$id)

# ---- 3. Pairwise downstream test ---------------------------------------------
# is_downstream(a, b) is TRUE when b lies downstream of a. The shared part of
# the FWA address tells us whether the two sit on the same drainage path. The
# local measure tells us whether b is further down that path than a.

is_downstream <- function(a, b) {

  row_a <- row_of_id[[as.character(a)]]
  row_b <- row_of_id[[as.character(b)]]

  addr_a <- segments[[row_a]]
  addr_b <- segments[[row_b]]

  # a downstream watershed cannot sit deeper in the tributary hierarchy
  if (length(addr_b) > length(addr_a)) return(FALSE)

  # b must share the start of a's FWA address
  if (!all(addr_b == addr_a[seq_along(addr_b)])) return(FALSE)

  # two coastal catchments can share the ocean as their outlet without ever
  # flowing into one another (Sunshine Coast vs North Shore)
  if (base$magnitude[row_a] == 0 && base$magnitude[row_b] == 0) return(FALSE)

  # flow only accumulates downstream, so b cannot carry less network than a
  if (base$magnitude[row_b] < base$magnitude[row_a]) return(FALSE)

  if (length(addr_b) == length(addr_a)) {
    # on the same tributary, b is downstream if it sits lower along it
    measure[row_b] < measure[row_a]

  } else {
    # b is on a downstream trunk. The next segment of a's address marks where
    # a's tributary joins b's stream, and b must be at or below that confluence
    confluence <- as.numeric(addr_a[length(addr_b) + 1])
    measure[row_b] <= confluence
  }
}

# ---- 4. NEXT_DOWN = adjacent downstream neighbour carrying the most flow ------
touch <- sf::st_touches(fwa)

# one row per (watershed, sharing a boundary with) pair
pairs <- data.frame(a = rep(base$id, lengths(touch)),
                    b = base$id[unlist(touch)])

pairs$flows <- mapply(is_downstream, pairs$a, pairs$b)
pairs$b_magnitude <- base$magnitude[row_of_id[as.character(pairs$b)]]

# one downstream neighbour per watershed, taking the highest magnitude as the
# trunk
next_down <- pairs |>
  dplyr::filter(flows) |>
  dplyr::group_by(a) |>
  dplyr::slice_max(b_magnitude, n = 1, with_ties = FALSE) |>
  dplyr::ungroup() |>
  dplyr::select(id = a, next_down = b)

# no downstream neighbour means a domain outlet, recorded as 0
base <- base |>
  dplyr::left_join(next_down, by = "id") |>
  dplyr::mutate(next_down = dplyr::coalesce(next_down, 0))

# ---- 4b. Reach length = measured stream length from FWA stream networks -------
# Preferred measure is the clipped length of the watershed's own main stem.
# Where that is missing, fall back to the longest single stream inside it.
# The stream network is read once and shared. 4b measures reach length here,
# and 4c reads the stream names to settle who owns the ground at the border.
.fwa_streams <- NULL
fwa_streams <- function(aoi_buf) {
  if (!is.null(.fwa_streams)) return(.fwa_streams)
  gdb <- data_path("fwa_stream_networks")
  message("reading AOI streams from FWA_STREAM_NETWORKS gdb…")

  # read only stream segments intersecting the buffered AOI
  bbox_wkt <- sf::st_as_text(sf::st_as_sfc(sf::st_bbox(aoi_buf)))
  layers <- sf::st_layers(gdb)
  spatial_layers <- layers$name[!is.na(layers$geomtype)]

  # keep only fields needed for stream-length assignment
  keep_cols <- c("FWA_WATERSHED_CODE", "BLUE_LINE_KEY", "LENGTH_METRE",
                 "GNIS_NAME")
  stream_layers <- list()

  for (layer in spatial_layers) {
    layer_streams <- sf::st_read(gdb, layer = layer, quiet = TRUE, wkt_filter = bbox_wkt)

    if (nrow(layer_streams) > 0) {
      stream_layers[[layer]] <- sf::st_zm(layer_streams[, intersect(keep_cols, names(layer_streams))])
    }
  }

  .fwa_streams <<- do.call(rbind, stream_layers)
  .fwa_streams
}

reach_from_streams <- function(fwa_sf, code_by_id, aoi_buf) { # polygon, lookup table, study area
  streams <- fwa_streams(aoi_buf)

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

  # primary reach length. Sum the clipped stream pieces whose FWA code matches
  # the watershed's own code. which() drops the unmatched rows, since a segment
  # with no code compares as NA and would otherwise subscript in a row of NAs.
  own_reach <- clip_tbl[
    which(clip_tbl$FWA_WATERSHED_CODE == code_by_id[as.character(clip_tbl$id)]),
  ] |>
    dplyr::group_by(id) |>
    dplyr::summarise(km = sum(clip_km), .groups = "drop")

  # fallback reach length. Use the single stream with the most clipped length
  # inside the watershed.
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

# ---- 4c. Splice in US (Whatcom) providers draining into the BC AOI ------------
# FWA coverage stops at the 49th parallel, but several Whatcom County
# watersheds (Chilliwack, Sumas, Silesia, Tomyhoi) drain north into the BC
# demand area. Pull those from the US national hydrography with nhdplusTools.
# WBD HUC12 gives the shapes and the routing, NHDPlus flowlines the reach
# length. This step keeps the units that route toward Canada, and step 5 drops
# anything too far upstream, so south-draining basins fall out on their own.

# Pair a US unit with the FWA watershed it sits under, not the one it nicks.
# Absolute overlap lets a large downstream unit (Sumas) win a 7 km² tail over
# the small FWA unit the US catchment actually abuts on 49°N. Prefer the
# longest shared 49th-parallel edge; if no candidate meets the border, take
# the FWA unit whose own area is most covered (f_fwa). Used for outlet_fwa
# (routing) and xb_partner (the 4d merge), which have to agree.
BORDER_LAT <- 49
BORDER_TOL_DEG <- 5e-4      # ~55 m
BORDER_EDGE_MIN_M <- 200    # shorter than this is a corner nick, not a sit-under
PARTNER_MIN_KM2 <- 0.5      # ignore overlap slivers when choosing a partner

border_edge_m <- function(us_g, fwa_g) {
  edge <- suppressWarnings(sf::st_intersection(
    sf::st_boundary(us_g), sf::st_boundary(fwa_g)))
  if (!length(edge) || all(sf::st_is_empty(edge))) return(0)
  if (any(sf::st_geometry_type(edge) == "GEOMETRYCOLLECTION"))
    edge <- sf::st_collection_extract(edge, "LINESTRING", warn = FALSE)
  edge <- edge[as.character(sf::st_geometry_type(edge)) %in%
                 c("LINESTRING", "MULTILINESTRING")]
  if (!length(edge) || all(sf::st_is_empty(edge))) return(0)
  edge_ll <- sf::st_transform(edge, 4326)
  band <- sf::st_as_sfc(sf::st_bbox(
    c(xmin = -180, ymin = BORDER_LAT - BORDER_TOL_DEG,
      xmax = 180,  ymax = BORDER_LAT + BORDER_TOL_DEG), crs = 4326))
  on <- suppressWarnings(sf::st_intersection(edge_ll, band))
  if (!length(on) || all(sf::st_is_empty(on))) return(0)
  as.numeric(sum(sf::st_length(sf::st_transform(on, PROJECT_CRS))))
}

# us_g: one US geometry. fwa_sf: candidate FWA polygons. ov_area: overlap m²
# aligned with fwa_sf rows. Returns a row index into fwa_sf, or NA.
best_fwa_partner <- function(us_g, fwa_sf, ov_area) {
  ok <- which(is.finite(ov_area) & ov_area >= PARTNER_MIN_KM2 * 1e6)
  if (!length(ok)) {
    ok <- which(is.finite(ov_area) & ov_area > 0)
    if (!length(ok)) return(NA_integer_)
  }
  us_g <- sf::st_geometry(us_g)
  fwa_g <- sf::st_geometry(fwa_sf)
  edge_m <- vapply(ok, function(j) border_edge_m(us_g, fwa_g[j]), numeric(1))
  if (any(edge_m >= BORDER_EDGE_MIN_M))
    return(ok[which.max(edge_m)])
  fwa_a <- as.numeric(sf::st_area(fwa_sf))
  f_fwa <- ov_area[ok] / pmax(fwa_a[ok], 1)
  ok[which.max(f_fwa)]
}

# columns carried on the spatial layer (must match the FWA layer's schema below)
us_cols <- c("WATERSHED_FEATURE_ID", "WATERSHED_ORDER", "WATERSHED_MAGNITUDE",
             "FWA_WATERSHED_CODE", "GNIS_NAME_1", "AREA_HA")

us <- tryCatch({
  message("pulling US WBD HUC12 + NHDPlus flowlines via nhdplusTools…")

  # US portion of the buffered AOI (the strip south of the 49th parallel).
  # fwa_buf is already an sfc (built with st_geometry above), so transform it
  # directly.
  buf_ll <- sf::st_transform(fwa_buf, 4326)
  south  <- sf::st_as_sfc(sf::st_bbox(
    c(xmin = -130, ymin = 40, xmax = -110, ymax = 49), crs = 4326))
  bb <- sf::st_bbox(sf::st_intersection(buf_ll, south))

  # WBD HUC12 over the strip. The WFS limits query size and refuses large or
  # complex areas, so cut the wide border strip into 2-degree longitude tiles,
  # pull each one by bounding box, and combine them. This takes the Washington
  # units and the Canadian ones they drain into, since those are what show a
  # unit routing across the border. Everything is clipped back afterwards.
  xcuts <- unname(unique(c(seq(bb["xmin"], bb["xmax"], by = 2), bb["xmax"])))
  ymin <- unname(bb["ymin"]); ymax <- unname(bb["ymax"])
  tiles <- Map(function(x0, x1) sf::st_as_sfc(sf::st_bbox(
    c(xmin = x0, ymin = ymin, xmax = x1, ymax = ymax), crs = 4326)),
    xcuts[-length(xcuts)], xcuts[-1])
  wbd_list <- lapply(tiles, function(tl)
    tryCatch(nhdplusTools::get_huc(AOI = tl, type = "huc12"),
             error = function(e) NULL))
  wbd_list <- wbd_list[!vapply(wbd_list, is.null, logical(1))]
  if (length(wbd_list) == 0) stop("WBD HUC12 pull returned nothing")
  wbd <- do.call(rbind, wbd_list)
  wbd <- wbd[!is.na(wbd$huc12) & !duplicated(wbd$huc12), ] |>
    sf::st_transform(PROJECT_CRS) |> sf::st_make_valid() |>
    # the WBD boundaries are very detailed, about 80 MB raw. Simplify to ~100 m
    # so the geometry work below and the saved layer stay manageable.
    sf::st_simplify(dTolerance = 100, preserveTopology = TRUE)
  wbd <- wbd[lengths(sf::st_intersects(wbd, fwa_buf)) > 0, ]
  wbd$huc12 <- as.character(wbd$huc12)
  wbd$tohuc <- as.character(wbd$tohuc)

  # Drop the marine and open-water units (hutype F = frontal, W = water). They
  # cover the Strait of Georgia and Puget Sound and overlap the coastal FWA
  # units, which would make it look as though water routes across the ocean.
  if ("hutype" %in% names(wbd)) {
    n_marine <- sum(wbd$hutype %in% c("F", "W"))
    if (n_marine > 0)
      message("  dropped ", n_marine, " marine/frontal WBD unit(s) (hutype F/W)")
    wbd <- wbd[!(wbd$hutype %in% c("F", "W")), ]
  }

  # A border outlet is a Washington unit that drains out of the US set and also
  # overlaps an FWA watershed. The unclipped WBD shapes reach past 49°N onto
  # FWA ground, so st_intersects finds the units whose streams cross on land.
  usid  <- wbd$huc12[grepl("WA", wbd$states)]
  to_us <- wbd$tohuc[match(usid, wbd$huc12)]
  exit_idx <- match(usid[!(to_us %in% usid)], wbd$huc12)   # units leaving the set
  overlaps <- lengths(sf::st_intersects(wbd[exit_idx, ], fwa)) > 0
  outlet_fwa <- list()
  message("  · pairing border outlets with the FWA unit they sit under ",
          "(49th-parallel edge, else f_fwa), not largest overlap")
  for (i in exit_idx[overlaps]) {
    hits <- sf::st_intersects(wbd[i, ], fwa)[[1]]
    if (length(hits) == 0) next
    cands <- fwa[hits, ]
    ov <- suppressWarnings(sf::st_intersection(
      sf::st_geometry(wbd[i, ]), sf::st_geometry(cands)))
    ov_area <- vapply(ov, function(g) {
      if (sf::st_is_empty(g)) return(0)
      as.numeric(sf::st_area(g))
    }, numeric(1))
    best <- best_fwa_partner(sf::st_geometry(wbd)[i], cands, ov_area)
    if (!is.na(best) && ov_area[best] > 0) {
      outlet_fwa[[wbd$huc12[i]]] <- cands$WATERSHED_FEATURE_ID[best]
      area_best <- which.max(ov_area)
      fwa_lab <- function(j) {
        nm <- cands$GNIS_NAME_1[j]
        if (is.na(nm) || !nzchar(nm)) paste0("#", cands$WATERSHED_FEATURE_ID[j])
        else nm
      }
      message(sprintf("      %-36s -> %-20s%s",
                      substr(wbd$name[i], 1, 36),
                      substr(fwa_lab(best), 1, 20),
                      if (!is.na(area_best) && area_best != best)
                        paste0("  (was ", fwa_lab(area_best), " by area)")
                      else ""))
    }
  }
  # provider set = the outlets plus everything upstream of them in the US graph
  eh <- data.frame(from = usid, to = to_us)
  eh <- eh[eh$to %in% usid, ]
  gU <- igraph::graph_from_data_frame(eh, vertices = data.frame(name = usid),
                                      directed = TRUE)
  keep <- unique(c(names(outlet_fwa), unlist(lapply(names(outlet_fwa),
    function(o) names(igraph::subcomponent(gU, o, mode = "in"))))))
  wbd_us <- wbd[wbd$huc12 %in% keep, ]

  # Everything above has to run on the unclipped US shapes. A watershed that
  # crosses the border only shows up as a border outlet because its shape
  # reaches north past 49°N onto ground the FWA also maps.
  #
  # That shared ground now has to go back to the FWA. About 660 km² north of
  # the line is mapped by both surveys, and if it stays in both, steps 3, 6 and
  # 9 count the same forest and the same rainfall twice.
  #
  # Handing back all of it goes too far, because some watersheds really do
  # straddle the line. Johnson Creek is 54 km², with 47 south of the border and
  # a 7 km² tail north of it that sits inside the FWA's Sumas River watershed.
  # Cut that tail off and the creek stops dead at the border while its own
  # headwaters are credited to a different river.
  #
  # So look at each overlap on its own and ask whether the FWA maps a named
  # stream there. If it does, the FWA has mapped that ground as part of a
  # watershed it recognises, and keeps it. If it maps only unnamed ditches, the
  # US survey has drawn the better catchment and the US unit keeps its tail.
  # Comparing areas cannot tell the two cases apart, and neither can the FWA
  # stream codes, which label every ditch on the old Sumas prairie with the
  # river it eventually reaches. LOBE_FWA_MAX stops the rule moving more than
  # a tail.
  #
  # Whatever a US unit keeps is also taken off the FWA watershed it came from,
  # just below the splice where `fwa` is in scope. Doing only one side of that
  # brings the double count straight back.
  LOBE_MIN_KM2 <- 0.5    # ignore slivers this small, they are not watersheds
  LOBE_FWA_MAX <- 0.25   # never move more than this share of an FWA watershed

  wbd_us <- wbd_us[!duplicated(wbd_us$huc12), ]
  fwa_near <- fwa[lengths(sf::st_intersects(fwa, wbd_us)) > 0, ]
  fwa_trim <- list()
  if (nrow(fwa_near) > 0) {
    before_km2 <- as.numeric(sf::st_area(wbd_us)) / 1e6
    fwa_area_m2 <- as.numeric(sf::st_area(fwa_near))
    # BC's streams over the border zone, which is what the test reads
    str_near <- fwa_streams(fwa_buf)
    str_near <- str_near[lengths(sf::st_intersects(
      str_near, sf::st_union(sf::st_geometry(wbd_us)))) > 0, ]
    message("  · ", nrow(str_near),
            " FWA stream segment(s) over the contested border zone")

    # Record which FWA watershed this US unit is the other half of, while the
    # US shape still holds its Canadian portion. Sit-under pairing (same rule
    # as outlet_fwa), not largest overlap: a 7 km² tail inside Sumas is not
    # the other half of Sumas. Unnamed FWA units still work.
    hits <- sf::st_intersects(wbd_us, fwa_near)
    xb_partner <- rep(NA_real_, nrow(wbd_us)); xb_frac <- rep(0, nrow(wbd_us))
    lobe_geom <- vector("list", nrow(wbd_us))   # ground each US unit keeps
    lobe_log  <- list()                         # one row per contested pair
    for (k in seq_len(nrow(wbd_us))) {
      h <- hits[[k]]; if (!length(h)) next
      gi <- lapply(h, function(q) suppressWarnings(sf::st_intersection(
        sf::st_geometry(wbd_us)[k], sf::st_geometry(fwa_near)[q])))
      a <- vapply(gi, function(g)
        if (!length(g)) 0 else sum(as.numeric(sf::st_area(g))), numeric(1))
      if (!length(a) || max(a) <= 0) next
      best <- best_fwa_partner(sf::st_geometry(wbd_us)[k], fwa_near[h, ], a)
      if (is.na(best)) next
      xb_partner[k] <- fwa_near$WATERSHED_FEATURE_ID[h[best]]
      xb_frac[k] <- a[best] / (before_km2[k] * 1e6)

      # apply the named-stream rule described above. For each overlap, check the
      # size floor and LOBE_FWA_MAX first, then look for a named FWA stream.
      f_us  <- a / (before_km2[k] * 1e6)
      f_fwa <- a / fwa_area_m2[h]
      take <- rep(FALSE, length(h))
      why  <- rep("", length(h))
      for (m in seq_along(h)) {
        if (a[m] < LOBE_MIN_KM2 * 1e6) { why[m] <- "below the size floor"; next }
        if (f_fwa[m] > LOBE_FWA_MAX) {
          why[m] <- sprintf("%.0f%% of the FWA unit, too much to move",
                            100 * f_fwa[m]); next
        }
        seg <- str_near[lengths(sf::st_intersects(str_near, gi[[m]])) > 0, ]
        named <- unique(stats::na.omit(seg$GNIS_NAME))
        if (length(named)) {
          why[m] <- paste0("FWA maps ", paste(named, collapse = ", "), " here")
          next
        }
        take[m] <- TRUE
        why[m] <- sprintf("%d unnamed ditch(es), FWA does not resolve it", nrow(seg))
      }
      for (m in which(a > 0 & f_us >= 0.02))
        lobe_log[[length(lobe_log) + 1]] <- data.frame(
          us = wbd_us$name[k], fwa_name = fwa_near$GNIS_NAME_1[h[m]],
          fwa_id = fwa_near$WATERSHED_FEATURE_ID[h[m]], km2 = a[m] / 1e6,
          f_us = f_us[m], f_fwa = f_fwa[m], moved = take[m],
          why = paste0(if (take[m]) "-> US   " else "-> FWA  ", why[m]))
      if (!any(take)) next
      lobe_geom[[k]] <- sf::st_make_valid(
        sf::st_union(do.call(c, gi[which(take)])))
      # this US unit keeps its tail and stays separate, so clear the overlap
      # record, or 4d merges back what was just separated
      if (!is.na(xb_partner[k]) &&
          xb_partner[k] %in% fwa_near$WATERSHED_FEATURE_ID[h[which(take)]]) {
        xb_partner[k] <- NA_real_; xb_frac[k] <- 0
      }
      # the same ground has to leave the FWA watershed it came from. Carry the
      # US unit's id so the subtraction below knows where it is going.
      for (m in which(take))
        fwa_trim[[length(fwa_trim) + 1]] <- sf::st_sf(
          WATERSHED_FEATURE_ID = as.numeric(fwa_near$WATERSHED_FEATURE_ID[h[m]]),
          us_id = as.numeric(wbd_us$huc12[k]),
          geom = sf::st_union(gi[[m]]))
    }
    message("  \u00b7 ", sum(xb_frac >= 0.15),
            " US unit(s) overlap a single FWA watershed by >=15% of their area")

    if (length(lobe_log)) {
      lg <- do.call(rbind, lobe_log)
      lg <- lg[order(-lg$km2), ]
      message("  \u00b7 contested ground, US unit vs the FWA watershed claiming it:")
      for (r in seq_len(nrow(lg)))
        message(sprintf("      %-36s %-20s %6.1f km2  f_us %.2f  f_fwa %.2f  %s",
                        substr(lg$us[r], 1, 36),
                        substr(ifelse(is.na(lg$fwa_name[r]),
                                      paste0("#", lg$fwa_id[r]), lg$fwa_name[r]),
                               1, 20),
                        lg$km2[r], lg$f_us[r], lg$f_fwa[r], lg$why[r]))
    }

    us_geom <- suppressWarnings(sf::st_difference(
      sf::st_geometry(wbd_us), sf::st_union(sf::st_geometry(fwa_near))))
    # hand the transferred lobes back to the US units they drain to
    for (k in which(!vapply(lobe_geom, is.null, logical(1))))
      us_geom[k] <- sf::st_union(us_geom[k], lobe_geom[[k]])
    us_geom <- sf::st_make_valid(us_geom)
    # the same collection trap as on the FWA side below, so keep it polygonal
    for (k in which(sf::st_geometry_type(us_geom) == "GEOMETRYCOLLECTION")) {
      p <- sf::st_collection_extract(us_geom[k], "POLYGON", warn = FALSE)
      us_geom[k] <- if (length(p) > 1) sf::st_combine(p) else p
    }
    sf::st_geometry(wbd_us) <- us_geom
    # A unit swallowed whole by the FWA footprint is entirely double-counted
    # and has nothing of its own left to contribute. Dropping it breaks the
    # tohuc chain of anything upstream, which then falls to next_down = 0 and
    # is trimmed in step 5 rather than routed on a unit that no longer exists.
    gone <- sf::st_is_empty(sf::st_geometry(wbd_us))
    if (any(gone)) {
      message("  dropped ", sum(gone),
              " US unit(s) lying wholly inside the FWA footprint")
      wbd_us <- wbd_us[!gone, ]
      before_km2 <- before_km2[!gone]
      xb_partner <- xb_partner[!gone]; xb_frac <- xb_frac[!gone]
    }
    wbd_us$areasqkm <- as.numeric(sf::st_area(wbd_us)) / 1e6
    handed_km2 <- sum(before_km2) - sum(wbd_us$areasqkm)
    message("  handed ", round(handed_km2, 1),
            " km² of overlapping BC territory back to the FWA (",
            sum(before_km2 - wbd_us$areasqkm > 0.01), " of ", nrow(wbd_us),
            " US unit(s) trimmed)")

    # one row per (FWA watershed, US unit) transfer, to apply after the splice
    moved_km2 <- 0
    if (length(fwa_trim)) {
      fwa_trim <- do.call(rbind, fwa_trim)
      moved_km2 <- sum(as.numeric(sf::st_area(fwa_trim))) / 1e6
      message("  · moved ", round(moved_km2, 1),
              " km² of cross-border tail the other way, to the US unit it",
              " drains to (", nrow(fwa_trim), " FWA watershed(s) trimmed)")
    } else {
      fwa_trim <- NULL
      message("  · no contested ground moved: BC's streams put every lobe in ",
              "the drainage of the FWA watershed already holding it")
    }
  }
  keep_ids <- wbd_us$huc12

  # NHDPlus flowlines over the kept footprint -> measured main-channel reach per
  # unit (mainstem = the highest stream order present in the unit). Replaces the
  # Hack's-law area proxy used for the old HydroBASINS splice.
  fl <- nhdplusTools::get_nhdplus(
    AOI = sf::st_as_sfc(sf::st_bbox(sf::st_transform(wbd_us, 4326))),
    realization = "flowline")
  ord_col <- grep("streamorde", names(fl), ignore.case = TRUE, value = TRUE)[1]
  fl <- sf::st_zm(fl) |> sf::st_transform(PROJECT_CRS)
  fl$ord <- if (!is.na(ord_col)) as.integer(fl[[ord_col]]) else 1L
  unit_poly <- wbd_us["huc12"]; names(unit_poly)[1] <- "id"
  clip <- suppressWarnings(sf::st_intersection(fl[, "ord"], unit_poly))
  clip <- clip[sf::st_geometry_type(clip) %in% c("LINESTRING", "MULTILINESTRING"), ]
  clip_tbl <- sf::st_drop_geometry(clip)
  clip_tbl$km <- as.numeric(sf::st_length(clip)) / 1000
  reach_us <- clip_tbl |>
    dplyr::group_by(id) |>
    dplyr::summarise(km = sum(km[ord == max(ord)], na.rm = TRUE), .groups = "drop")
  reach_km_us <- reach_us$km[match(keep_ids, reach_us$id)]
  reach_km_us <- ifelse(is.na(reach_km_us) | reach_km_us <= 0,
                        1.4 * wbd_us$areasqkm^0.6, reach_km_us)

  # Set the downstream neighbour. A unit inside the US set points at the next US
  # unit, and a border outlet points at the FWA watershed recorded above. Step 5
  # then drops any whose chain never reaches the demand area.
  nd <- rep(0, nrow(wbd_us))
  internal <- wbd_us$tohuc %in% keep_ids
  nd[internal] <- as.numeric(wbd_us$tohuc[internal])
  out_match <- match(keep_ids, names(outlet_fwa))
  nd[!is.na(out_match)] <- as.numeric(unlist(outlet_fwa)[out_match[!is.na(out_match)]])

  list(
    attr = data.frame(
      id = as.numeric(keep_ids), magnitude = NA_real_,
      area_km2 = wbd_us$areasqkm, fwa = NA_character_, loc = NA_character_,
      next_down = nd, reach_km = reach_km_us,
      xb_partner = if (exists("xb_partner")) xb_partner else NA_real_,
      xb_frac = if (exists("xb_frac")) xb_frac else 0),
    fwa_trim = fwa_trim,
    stats = list(handed_km2 = if (exists("handed_km2")) handed_km2 else NA_real_,
                 moved_km2  = if (exists("moved_km2"))  moved_km2  else NA_real_),
    polys = wbd_us |>
      dplyr::transmute(
        WATERSHED_FEATURE_ID = as.numeric(huc12), WATERSHED_ORDER = NA_real_,
        WATERSHED_MAGNITUDE = NA_real_, FWA_WATERSHED_CODE = NA_character_,
        GNIS_NAME_1 = as.character(name), AREA_HA = areasqkm * 100))
}, error = function(e) {
  message("  ! US provider pull skipped: ", conditionMessage(e)); NULL
})

if (!is.null(us) && nrow(us$attr) > 0) {
  # The other half of the handover in 4c. Take each piece of ground a US unit
  # kept off the FWA watershed that also claimed it. This runs here rather than
  # in 4c, which builds the US layer and cannot reach `fwa`.
  if (inherits(us$fwa_trim, "sf") && nrow(us$fwa_trim) > 0) {
    # The two surveys drew the shared boundary with different vertices, so
    # subtracting a piece leaves hairline fragments along it. Sumas River had 27
    # of them, the largest 0.05 km². Left in, they turn the watershed into a
    # 28-part shape that GDAL's RFC7946 writer refuses, and it drops off the map
    # entirely. Send them along with the piece, so both watersheds stay single
    # shapes covering the same ground.
    SLIVER_M2 <- 1e5
    gcol <- sf::st_geometry(fwa)
    ucol <- sf::st_geometry(us$polys)
    moved <- 0; sliver_km2 <- 0
    for (r in seq_len(nrow(us$fwa_trim))) {
      j <- match(us$fwa_trim$WATERSHED_FEATURE_ID[r], fwa$WATERSHED_FEATURE_ID)
      k <- match(us$fwa_trim$us_id[r], us$polys$WATERSHED_FEATURE_ID)
      if (is.na(j)) next
      g <- suppressWarnings(sf::st_difference(
        gcol[j], sf::st_geometry(us$fwa_trim)[r]))
      # Order matters here. st_make_valid is what turns a difference into a
      # GEOMETRYCOLLECTION, so extract after it. Extract first and a collection
      # reaches the gpkg, where exactextractr refuses the whole layer in 10.
      g <- sf::st_make_valid(g)
      if (any(sf::st_geometry_type(g) == "GEOMETRYCOLLECTION"))
        g <- sf::st_collection_extract(g, "POLYGON", warn = FALSE)
      if (!length(g) || all(sf::st_is_empty(g))) next
      parts <- suppressWarnings(sf::st_cast(g, "POLYGON"))
      pa <- as.numeric(sf::st_area(parts))
      big <- pa >= SLIVER_M2
      if (!any(big)) next
      gcol[j] <- if (sum(big) > 1) sf::st_combine(parts[big]) else parts[big]
      if (any(!big) && !is.na(k)) {
        ucol[k] <- sf::st_make_valid(
          sf::st_union(ucol[k], sf::st_combine(parts[!big])))
        sliver_km2 <- sliver_km2 + sum(pa[!big]) / 1e6
      }
      # Subtract from AREA_HA rather than remeasuring the shape. AREA_HA is the
      # FWA's own surveyed figure and every other area here comes from it, so
      # remeasuring would fold a survey-versus-shape difference into the move.
      lost_km2 <- (as.numeric(sf::st_area(sf::st_geometry(us$fwa_trim)[r])) +
                     sum(pa[!big])) / 1e6
      fwa$AREA_HA[j] <- max(fwa$AREA_HA[j] - lost_km2 * 100, 0)
      moved <- moved + lost_km2
      bj <- match(fwa$WATERSHED_FEATURE_ID[j], base$id)
      if (!is.na(bj)) base$area_km2[bj] <- fwa$AREA_HA[j] / 100
      if (!is.na(k)) {
        us$polys$AREA_HA[k] <- us$polys$AREA_HA[k] + sum(pa[!big]) / 1e4
        us$attr$area_km2[k] <- us$attr$area_km2[k] + sum(pa[!big]) / 1e6
      }
    }
    sf::st_geometry(fwa) <- gcol
    sf::st_geometry(us$polys) <- ucol
    message("  · ", round(moved, 1), " km² subtracted from ", nrow(us$fwa_trim),
            " FWA watershed(s) that the transferred tails came from (",
            round(sliver_km2, 3), " km² of that as boundary slivers)")
  }
  base <- dplyr::bind_rows(base, us$attr)
  fwa  <- rbind(fwa[, us_cols], us$polys[, us_cols])
  message("spliced ", nrow(us$attr),
          " US (WBD HUC12) providers; flow-distance trim decides which survive")
}

# ---- 4d. Reunite watersheds that only the border divides ---------------------
# A stream crossing 49°N is mapped twice, once by each survey. Saar Creek is an
# FWA assessment watershed AND a WBD subwatershed, two shapes whose only shared
# edge is the survey boundary. Merge a pair like that back into one watershed.
#
# Both parts of the rule are needed.
#   1. Their shared boundary really is the border, measured as the share of it
#      sitting on 49°N. Below BORDER_SHARE_MIN the two meet at a real ridge.
#   2. They are the same named stream. Touching alone is not enough, since
#      several different US units sit against the same BC one, and merging on
#      that would swallow real headwater units into their downstream neighbour.
#
# norm_stream() drops WBD's "Bear Creek-Chilliwack River" compounds, the
# Upper/Lower/Headwater prefixes and the Creek/River suffix, and uses an alias
# table for the streams the two countries spell differently. Names then have to
# match exactly.
BORDER_SHARE_MIN <- 0.9
DIVIDE_MIN <- 0.5           # below this share of ridge, the edge is no divide
XB_FRAC_MIN <- 0.15         # share of a US unit that lay inside one FWA unit

dsm <- tryCatch(terra::rast(data_path("dem_glo30")), error = function(e) NULL)

CROSS_BORDER_ALIASES <- c(silesia = "slesse", tomyhoi = "tamihi",
                          ensawkwatch = "nesakwatch")

norm_stream <- function(x) {
  x <- tolower(trimws(ifelse(is.na(x), "", x)))
  x <- sub("^.*-", "", x)                                   # WBD compound names
  x <- gsub("^(upper|lower|middle|headwaters?|north|south|east|west) ", "", x)
  x <- trimws(gsub(" (creek|river|brook|slough)$", "", x))
  ifelse(x %in% names(CROSS_BORDER_ALIASES), CROSS_BORDER_ALIASES[x], x)
}

if (!is.null(us) && nrow(us$attr) > 0) {
  us_i <- which(is.na(fwa$FWA_WATERSHED_CODE))
  bc_i <- which(!is.na(fwa$FWA_WATERSHED_CODE))
  nm <- norm_stream(fwa$GNIS_NAME_1)
  ll <- sf::st_transform(fwa, 4326)
  touch <- sf::st_intersects(fwa[us_i, ], fwa[bc_i, ])

  merges <- list()
  for (a in seq_along(us_i)) {
    for (b in touch[[a]]) {
      i <- us_i[a]; j <- bc_i[b]

      # 1. their shared boundary is the border, measured on the lon/lat copy
      edge_ll <- suppressWarnings(sf::st_intersection(
        sf::st_boundary(sf::st_geometry(ll)[i]),
        sf::st_boundary(sf::st_geometry(ll)[j])))
      if (any(sf::st_geometry_type(edge_ll) == "GEOMETRYCOLLECTION"))
        edge_ll <- sf::st_collection_extract(edge_ll, "LINESTRING", warn = FALSE)
      edge_ll <- edge_ll[as.character(sf::st_geometry_type(edge_ll)) %in%
                           c("LINESTRING", "MULTILINESTRING")]
      if (!length(edge_ll) || all(sf::st_is_empty(edge_ll))) next
      co_ll <- sf::st_coordinates(edge_ll)
      on_border <- mean(abs(co_ll[, 2] - BORDER_LAT) < BORDER_TOL_DEG)

      # The strongest evidence, and it stands on its own. Before the clip, this
      # US unit covered part of the very FWA watershed now being tested, so both
      # surveys drew the same ground. It skips the border-share test the weaker
      # checks need.
      bi <- match(fwa$WATERSHED_FEATURE_ID[i], base$id)
      overlapped <- isTRUE(!is.na(base$xb_partner[bi]) &&
        base$xb_partner[bi] == fwa$WATERSHED_FEATURE_ID[j] &&
        base$xb_frac[bi] >= XB_FRAC_MIN)

      if (!overlapped && on_border < BORDER_SHARE_MIN) next

      # 2. they are the same watershed, on any one of three kinds of evidence
      same_name <- nzchar(nm[i]) && nm[i] == nm[j]
      drains_in <- isTRUE(base$next_down[bi] == fwa$WATERSHED_FEATURE_ID[j])

      # The third test asks whether the shared edge is really a drainage divide.
      # Sample the elevation along it and 300 m either side. Where the edge is
      # rarely the local high point it is a hillslope rather than a ridge, and
      # nothing but the survey line separates the two watersheds. Sweltzer and
      # Frosst measure as ridge over 6% of their edge and Kendall's over 29%,
      # against close to 100% for a real ridgeline. This reads the raw DSM,
      # because 05_dem.tif is built from this script's own AOI and using it here
      # would be circular.
      not_a_divide <- FALSE
      if (!overlapped && !same_name && !drains_in && !is.null(dsm)) {
        edge_m <- sf::st_transform(edge_ll, PROJECT_CRS)
        pts <- suppressWarnings(sf::st_cast(sf::st_line_sample(
          sf::st_cast(edge_m, "LINESTRING", warn = FALSE), density = 1 / 300), "POINT"))
        if (length(pts) >= 5) {
          co <- sf::st_coordinates(pts)[, 1:2, drop = FALSE]
          zz <- function(dy) terra::extract(dsm, terra::vect(
            cbind(co[, 1], co[, 2] + dy), crs = PROJECT_CRS, type = "points"),
            ID = FALSE)[, 1]
          z0 <- zz(0); zN <- zz(300); zS <- zz(-300)
          good <- is.finite(z0) & is.finite(zN) & is.finite(zS)
          if (sum(good) >= 5)
            not_a_divide <- mean(z0[good] >= pmax(zN[good], zS[good]) - 5) < DIVIDE_MIN
        }
      }

      if (!overlapped && !same_name && !drains_in && !not_a_divide) next
      merges[[length(merges) + 1]] <- list(
        us = i, bc = j,
        # rank the evidence, strongest first
        rank = if (overlapped) 0L else if (same_name) 1L else
          if (drains_in) 2L else 3L,
        km = sum(as.numeric(sf::st_length(sf::st_transform(edge_ll, PROJECT_CRS)))))
    }
  }

  # Only one US unit may join each BC unit. Without that limit, every US unit
  # draining into a BC one gets absorbed by it and separate tributaries collapse
  # together, which grew Sumas River to 205 km², 69% of it American. Keep the
  # partner with the strongest evidence, and break ties on the longest shared
  # border.
  if (length(merges) > 0) {
    ord <- order(vapply(merges, `[[`, integer(1), "rank"),
                 -vapply(merges, `[[`, numeric(1), "km"))
    merges <- merges[ord]
    seen_bc <- integer(0); pick <- list()
    for (m in merges) {
      if (m$bc %in% seen_bc) next
      seen_bc <- c(seen_bc, m$bc); pick[[length(pick) + 1]] <- m
    }
    dropped_n <- length(merges) - length(pick)
    if (dropped_n > 0)
      message("  held back ", dropped_n,
              " further US unit(s) that drain into an already-reunited watershed")
    merges <- pick
    keep <- rep(TRUE, nrow(fwa))
    g <- sf::st_geometry(fwa)
    for (m in merges) {
      i <- m$us; j <- m$bc
      if (!keep[i]) next
      g[j] <- sf::st_union(g[j], g[i])
      # many FWA units have no GNIS name, so take the US half's where there is
      # one rather than falling back to "Upstream area"
      if (is.na(fwa$GNIS_NAME_1[j]) || !nzchar(fwa$GNIS_NAME_1[j]))
        fwa$GNIS_NAME_1[j] <- fwa$GNIS_NAME_1[i]
      # the reunited watershed carries both halves' area and mainstem length
      fwa$AREA_HA[j] <- fwa$AREA_HA[j] + fwa$AREA_HA[i]
      bj <- match(fwa$WATERSHED_FEATURE_ID[j], base$id)
      bi <- match(fwa$WATERSHED_FEATURE_ID[i], base$id)
      base$area_km2[bj] <- base$area_km2[bj] + base$area_km2[bi]
      base$reach_km[bj] <- base$reach_km[bj] + base$reach_km[bi]
      # anything that drained into the US half now drains into the whole
      base$next_down[base$next_down == base$id[bi]] <- base$id[bj]
      keep[i] <- FALSE
    }
    sf::st_geometry(fwa) <- sf::st_make_valid(g)
    dropped <- fwa$WATERSHED_FEATURE_ID[!keep]
    named <- unique(norm_stream(fwa$GNIS_NAME_1[!keep]))
    fwa <- fwa[keep, ]
    base <- base[!(base$id %in% dropped), ]
    N_REUNITED <- length(dropped)
    message("reunited ", length(dropped),
            " watershed(s) split by the 49th parallel: ",
            paste(sort(named), collapse = ", "))
  }
}

# ---- 4e. Rehome fragments this script stranded --------------------------------
# The cutting in 4c and the merging in 4d can leave a watershed holding a piece
# of ground that no longer touches the rest of it. Fishtrap Creek ended up with
# a 2 km² parcel 178 m away across Johnson Creek. Fishtrap drains south to the
# Nooksack and stops there, so that parcel's retention went to a watershed that
# delivers nothing downstream, while everything around it drains north into the
# demand area.
#
# The rule is deliberately narrow, because watersheds in several pieces are
# normal here. 41 of the 613 FWA units arrive from the survey that way and must
# be left alone. Move a piece only when both of these hold.
#   · It does not touch the part of the watershed that carries the survey's own
#     shape.
#   · It does not overlap that shape at all, so it is ground this script
#     attached rather than ground the FWA put there.
# The piece then goes to whichever neighbour shares the longest boundary with
# it. If no neighbour will take it, leave it where it is, so the watersheds
# still cover the same ground between them.
rehomed <- 0; rehomed_km2 <- 0
gcol <- sf::st_geometry(fwa)
for (i in seq_len(nrow(fwa))) {
  # Single brackets keep the lookup an sfc with its CRS. `[[` would hand back a
  # bare sfg and every st_intersection below would refuse it.
  si <- match(as.character(fwa$WATERSHED_FEATURE_ID[i]), names(fwa_source_geom))
  if (is.na(si)) next                        # a US unit, no FWA shape to judge
  src <- fwa_source_geom[si]
  parts <- suppressWarnings(sf::st_cast(gcol[i], "POLYGON"))
  if (length(parts) < 2) next

  ov <- vapply(seq_along(parts), function(p) {
    i2 <- suppressWarnings(sf::st_intersection(parts[p], src))
    if (length(i2) == 0) 0 else sum(as.numeric(sf::st_area(i2)))
  }, numeric(1))
  if (all(ov <= 0)) next                     # nothing anchors it, so leave it
  anchor <- which.max(ov)
  loose <- which(ov <= 0 &
                 lengths(sf::st_intersects(parts, parts[anchor])) == 0)
  if (!length(loose)) next

  for (p in loose) {
    frag <- parts[p]
    # The longest shared boundary wins. st_touches would miss a neighbour that
    # sits across a hairline gap, so measure the shared boundary directly.
    cand <- setdiff(which(lengths(sf::st_intersects(gcol, frag)) > 0), i)
    if (!length(cand)) next
    shared <- vapply(cand, function(k) {
      b <- suppressWarnings(sf::st_intersection(sf::st_boundary(gcol[k]),
                                                sf::st_boundary(frag)))
      if (length(b) == 0) 0 else sum(as.numeric(sf::st_length(b)))
    }, numeric(1))
    if (max(shared) <= 0) next
    to <- cand[which.max(shared)]
    frag_km2 <- as.numeric(sf::st_area(frag)) / 1e6

    gcol[i]  <- sf::st_make_valid(sf::st_difference(gcol[i], frag))
    gcol[to] <- sf::st_make_valid(sf::st_union(gcol[to], frag))

    # move the surveyed area along with the ground, as 4c does, so AREA_HA stays
    # the figure every later area comes from
    fwa$AREA_HA[i]  <- max(fwa$AREA_HA[i] - frag_km2 * 100, 0)
    fwa$AREA_HA[to] <- fwa$AREA_HA[to] + frag_km2 * 100
    bi <- match(fwa$WATERSHED_FEATURE_ID[i], base$id)
    bt <- match(fwa$WATERSHED_FEATURE_ID[to], base$id)
    if (!is.na(bi)) base$area_km2[bi] <- fwa$AREA_HA[i] / 100
    if (!is.na(bt)) base$area_km2[bt] <- fwa$AREA_HA[to] / 100

    rehomed <- rehomed + 1; rehomed_km2 <- rehomed_km2 + frag_km2
    message("  · rehomed ", round(frag_km2, 2), " km² from ",
            fwa$WATERSHED_FEATURE_ID[i], " (",
            ifelse(is.na(fwa$GNIS_NAME_1[i]), "unnamed", fwa$GNIS_NAME_1[i]),
            ") to ", fwa$WATERSHED_FEATURE_ID[to], " (",
            ifelse(is.na(fwa$GNIS_NAME_1[to]), "unnamed", fwa$GNIS_NAME_1[to]),
            ")")
  }
}
sf::st_geometry(fwa) <- gcol
message("rehomed ", rehomed, " stranded fragment(s), ",
        round(rehomed_km2, 2), " km² total")

# ---- 4f. Repair NEXT_DOWN where the splice broke adjacency -------------------
# Step 4 picks the downstream neighbour from the uncut FWA shapes. 4c–4e then
# cut and merge along the 49th parallel, so a unit can still point at an outlet
# it no longer touches. Unnamed 8172 is a Sumas tributary in the FWA address
# but sits against Johnson Creek; the map drew a hop over the land in between.
# Point any such unit at a neighbour that still reaches that same outlet.
ws_lab <- function(id) {
  i <- match(id, fwa$WATERSHED_FEATURE_ID)
  if (is.na(i)) return(as.character(id))
  nm <- fwa$GNIS_NAME_1[i]
  if (is.na(nm) || !nzchar(nm)) paste0("#", id) else nm
}
reaches <- function(start, target, nd, ids) {
  seen <- integer(0)
  x <- start
  for (step in seq_len(length(ids) + 1L)) {
    if (is.na(x) || identical(x, 0) || identical(x, 0L)) return(FALSE)
    if (identical(x, target)) return(TRUE)
    if (x %in% seen) return(FALSE)
    seen <- c(seen, x)
    j <- match(x, ids)
    if (is.na(j)) return(FALSE)
    x <- nd[j]
  }
  FALSE
}
contact_score <- function(a, b) {
  xi <- suppressWarnings(sf::st_intersection(a, b))
  area <- if (!length(xi) || all(sf::st_is_empty(xi))) 0
    else sum(as.numeric(sf::st_area(xi)))
  be <- suppressWarnings(sf::st_intersection(sf::st_boundary(a),
                                            sf::st_boundary(b)))
  if (length(be) && any(sf::st_geometry_type(be) == "GEOMETRYCOLLECTION"))
    be <- sf::st_collection_extract(be, "LINESTRING", warn = FALSE)
  len <- if (!length(be) || all(sf::st_is_empty(be))) 0
    else sum(as.numeric(sf::st_length(be)))
  c(area = area, len = len)
}

ids <- base$id
nd <- base$next_down
g <- sf::st_geometry(fwa)
gi <- match(ids, fwa$WATERSHED_FEATURE_ID)
OVERLAP_M2 <- 1e3
rerouted <- 0L
for (k in which(!is.na(nd) & nd != 0)) {
  i <- gi[k]
  j <- gi[match(nd[k], ids)]
  if (is.na(i) || is.na(j)) next
  if (lengths(sf::st_intersects(g[i], g[j])) > 0) next

  hits <- setdiff(sf::st_intersects(g[i], g)[[1]], i)
  if (!length(hits)) next
  sc <- vapply(hits, function(o) contact_score(g[i], g[o]), numeric(2))
  keep <- sc["area", ] > 0 | sc["len", ] > 0
  hits <- hits[keep]
  sc <- sc[, keep, drop = FALSE]
  if (!length(hits)) next

  hit_ids <- fwa$WATERSHED_FEATURE_ID[hits]
  ok <- vapply(hit_ids, function(oid)
    identical(oid, nd[k]) || reaches(oid, nd[k], nd, ids), logical(1))
  hits <- hits[ok]
  sc <- sc[, ok, drop = FALSE]
  if (!length(hits)) next

  pick <- if (any(sc["area", ] > OVERLAP_M2))
    hits[which.max(sc["area", ])] else hits[which.max(sc["len", ])]
  new_nd <- fwa$WATERSHED_FEATURE_ID[pick]
  if (identical(new_nd, nd[k])) next
  if (reaches(new_nd, ids[k], nd, ids)) next

  message("  · rerouted ", ws_lab(ids[k]), " (", ids[k], ") ",
          ws_lab(nd[k]), " -> ", ws_lab(new_nd), " (", new_nd, ")")
  nd[k] <- new_nd
  rerouted <- rerouted + 1L
}
base$next_down <- nd
message("rerouted ", rerouted,
        " watershed(s) whose downstream neighbour no longer adjoins them")

# ---- 5. Trim to ≤ MAX_FLOW_DIST_KM upstream of the downstream AOI -------------
# copy NEXT_DOWN and reach length onto the spatial layer
fwa$NEXT_DOWN <- base$next_down[match(fwa$WATERSHED_FEATURE_ID, base$id)]
fwa$reach_km <- base$reach_km[match(fwa$WATERSHED_FEATURE_ID, base$id)]

# all watershed ids = graph vertices
all_ids <- base$id

# directed flow edges, each watershed to its downstream neighbour, weighted by
# reach length in km
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
# A true island or sliver such as Bowen or Gambier shares no boundary with any
# other watershed AND has no flow connection into the study area. Shape alone
# is not enough to tell, because the US units come from a different survey and
# share no edge with the BC ones even though they do drain in.
geom_isolated <- lengths(sf::st_touches(fwa)) == 0
flows_out <- fwa$NEXT_DOWN %in% fwa$WATERSHED_FEATURE_ID   # drains into a kept unit
flows_in  <- fwa$WATERSHED_FEATURE_ID %in% base$next_down  # a kept unit drains into it
island <- geom_isolated & !flows_out & !flows_in

if (any(island)) {
  message("dropped ", sum(island), " drainage-isolated island/sliver watershed(s)")
  # drop islands as they do not attenuate any downstream areas
  fwa <- fwa[!island, ]
  # trim the attribute table to match. Nothing routes into an island, so this
  # leaves no references pointing at a watershed that is gone.
  base <- base[base$id %in% fwa$WATERSHED_FEATURE_ID, ]
}

# ---- 5c. Lake barriers: sever routing through large lakes --------------------
# The FWA has no field marking a lake, so this restores the lake-as-sink rule
# that HydroBASINS supplied. A big lake or reservoir takes up the storm pulse in
# its own storage, so holding water back above it buys little downstream. Flag
# each lake that qualifies and cut its downstream connection.
#
# DEFERRED. The runoff held back above a lake should be credited to the LAKE
# rather than discarded. That needs the lake separated from the watershed land
# it sits in, and a reservoir's capacity handled differently from a natural
# lake's overflow. Until then lakes stay terminal, which at least avoids
# crediting a watershed with retention it does not provide.
base$is_lake_barrier <- FALSE
lakes <- tryCatch({
  lakes_all <- bcdc_query_geodata("WHSE_BASEMAPPING.FWA_LAKES_POLY") |>
    filter(INTERSECTS(fwa_buf)) |> collect()
  # a lake qualifies on surface area, or by being a named reservoir
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

  # A watershed becomes a barrier in either of two ways. A qualifying lake
  # covers at least LAKE_COVER_FRAC of it, which handles lakes spread across
  # several watersheds. Or the lake overlaps this watershed more than any other,
  # which catches small reservoirs sitting inside a larger mixed watershed, such
  # as Coquitlam, Alouette and Capilano.

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

# flag the outlets, meaning watersheds with no downstream neighbour in the model
fwa$is_sink <- fwa$NEXT_DOWN == 0

# tag each watershed by regional district
# load regional district polygons and project them to the working CRS
districts <- bcmaps::regional_districts(ask = FALSE) |>
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

# the US watersheds have no FWA code, so tag them separately
fwa$admin_district[is.na(fwa$FWA_WATERSHED_CODE)] <- "Whatcom (US)"

# ---- 7. Persist sub-basin layer ----------------------------------------------
# columns needed by downstream scripts, plus selected original FWA identifiers
keep_cols <- c("HYBAS_ID", "NEXT_DOWN", "SUB_AREA", "flow_dist_km", "reach_km",
               "in_downstream_aoi", "is_sink", "is_lake_barrier", "admin_district",
               "WATERSHED_FEATURE_ID", "WATERSHED_ORDER", "WATERSHED_MAGNITUDE",
               "FWA_WATERSHED_CODE", "GNIS_NAME_1")

# Nothing but polygons leaves this script. Every st_make_valid above, in the 4c
# handover, the 4d merges and the lake severing, can hand back a
# GEOMETRYCOLLECTION, and that fails quietly further on. exactextractr rejects
# the whole layer in step 10, and GDAL's RFC7946 writer drops the watershed from
# the map. Fix it once here, where every path out of the script passes, rather
# than guarding each step.
gcol <- sf::st_geometry(fwa)
gc_i <- which(sf::st_geometry_type(gcol) == "GEOMETRYCOLLECTION")
for (i in gc_i) {
  p <- sf::st_collection_extract(gcol[i], "POLYGON", warn = FALSE)
  gcol[i] <- if (length(p) > 1) sf::st_combine(p) else p
}
if (length(gc_i)) {
  sf::st_geometry(fwa) <- gcol
  message("normalised ", length(gc_i),
          " geometry collection(s) to polygons before writing")
}

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
  # a mitre join keeps the sharp watershed corners that the default buffer
  # would round off
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
# Two checks. The left panel asks whether the routing tree branches like a real
# river network and converges on the outlets. The right asks whether flow
# distance falls smoothly toward the demand area, since that field drives the
# 100 km trim in step 5.
QA_PAL <- c(bc = "#a6cbe6", bc_dark = "#4f83ad", us = "#e9d08a",
            us_dark = "#b8891f", ink = "#0f1722", slate = "#5b6b7e",
            hot = "#df744a", teal = "#2a9d8f", lake = "#2c7fb8")

# Use the coast-clipped demand boundary, and open the frame wide enough to hold
# both it and the sub-basins. Framed on the sub-basins alone, the boundary's
# western limb runs off the panel and appears to stop in mid-ocean.
aoi_disp <- aoi_display()
qa_view <- span_bbox(fwa, aoi_disp)

qa_png("02_subbasins_fwa.png", ncol = 2, panel_w = 1250, panel_h = 1250,
       plot_fn = function() {
  op <- graphics::par(mfrow = c(1, 2))
  on.exit(graphics::par(op), add = TRUE)

  # ---- panel 1: routing tree and terminal units ------------------------------
  # A watershed with no downstream neighbour is a real outlet if it sits inside
  # the demand area. Outside it, the real downstream was simply trimmed away.
  is_lake   <- fwa$is_lake_barrier %in% TRUE
  is_demand <- fwa$is_sink & !is_lake & (fwa$in_downstream_aoi %in% TRUE)
  is_edge   <- fwa$is_sink & !is_lake & !(fwa$in_downstream_aoi %in% TRUE)

  fwa$terminal <- ifelse(is_demand, "demand outlet",
                  ifelse(is_edge, "domain-edge sink",
                  ifelse(is_lake, "lake barrier", "provider")))

  # terra orders classes alphabetically and ignores factor levels, so the
  # colours are listed in that order rather than in the order above
  terra::plot(terra::vect(fwa), "terminal", type = "classes",
              col = c("demand outlet"    = "#f6c6b0",
                      "domain-edge sink" = "#f0dca8",
                      "lake barrier"     = QA_PAL[["bc"]],
                      "provider"         = "#f2f5f8"),
              border = "white", lwd = 0.15, axes = TRUE,
              xlim = qa_view$xlim, ylim = qa_view$ylim,
              main = "Routing tree and terminal units")

  # routing edges, drawn faintly so the terminal fills stay readable
  centroids <- sf::st_coordinates(suppressWarnings(
    sf::st_centroid(sf::st_geometry(fwa))))
  rownames(centroids) <- as.character(fwa$HYBAS_ID)
  edge_df <- data.frame(id = fwa$HYBAS_ID, nd = fwa$NEXT_DOWN)
  edge_df <- edge_df[edge_df$nd != 0 & edge_df$nd %in% fwa$HYBAS_ID, ]
  graphics::segments(centroids[as.character(edge_df$id), "X"],
                     centroids[as.character(edge_df$id), "Y"],
                     centroids[as.character(edge_df$nd), "X"],
                     centroids[as.character(edge_df$nd), "Y"],
                     col = "#0f172259", lwd = 0.55)
  plot(sf::st_geometry(aoi_disp), add = TRUE,
       border = QA_PAL[["ink"]], lwd = 2)

  # ---- panel 2: flow distance to the demand area -----------------------------
  terra::plot(terra::vect(fwa), "flow_dist_km", type = "interval",
              breaks = c(0, 20, 40, 60, 80, 100),
              col = c("#eef4f9", "#a6cbe6", "#4f83ad", "#1f4e73", "#0f1722"),
              border = "grey85", lwd = 0.1, axes = TRUE,
              xlim = qa_view$xlim, ylim = qa_view$ylim,
              main = "Flow distance to the demand area (km)")
  plot(sf::st_geometry(aoi_disp), add = TRUE, border = "black", lwd = 2)
})

# ---- QA: border stitching diagnostic ----------------------------------------
# Checks whether the two national surveys meet cleanly along the 49th parallel.
# FWA watersheds are blue, US ones amber, and the routing edges carrying water
# north across the line are orange. A gap between blue and amber means the join
# failed.
#
# Each edge runs from the US watershed's centre to the point where it meets the
# FWA watershed it drains into, rather than to that watershed's centre. A
# centre-to-centre line would track over whatever watersheds lie between and
# read as water crossing ground it never touches. The receiving watersheds are
# tinted a darker shade to show the pairing instead.
is_us <- is.na(fwa$FWA_WATERSHED_CODE)
if (any(is_us)) {
  qa_png("02_border_stitch.png", panel_w = 1500, panel_h = 1200,
         plot_fn = function() {
    centroids <- sf::st_coordinates(suppressWarnings(
      sf::st_point_on_surface(sf::st_geometry(fwa))))
    rownames(centroids) <- as.character(fwa$HYBAS_ID)
    edge_df <- data.frame(id = fwa$HYBAS_ID, nd = fwa$NEXT_DOWN, src_us = is_us)
    edge_df <- edge_df[edge_df$nd != 0 & edge_df$nd %in% fwa$HYBAS_ID, ]
    cross <- edge_df[edge_df$src_us &
                       !is_us[match(edge_df$nd, fwa$HYBAS_ID)], ]

    # Frame on the watersheds taking part in the join. The full US extent runs
    # south to Bellingham Bay and has nothing to do with the border.
    join_ids <- unique(c(cross$id, cross$nd))
    focus <- fwa[fwa$HYBAS_ID %in% join_ids, ]
    fb <- sf::st_bbox(if (nrow(focus) > 0) focus else fwa[is_us, ])
    pad <- 0.1 * (fb[["xmax"]] - fb[["xmin"]])

    role <- ifelse(!is_us & fwa$HYBAS_ID %in% cross$nd, "FWA (BC), receiving",
            ifelse(!is_us, "FWA (BC)",
            ifelse(fwa$HYBAS_ID %in% cross$id, "WBD (US), draining north",
                   "WBD (US)")))
    # Name the colours. terra drops any class nothing falls into and then reads
    # the colour vector by position, so with no cross-border edge the two
    # "draining" classes vanish and an unnamed vector would recolour the rest.
    survey_cols <- c("FWA (BC)"                 = QA_PAL[["bc"]],
                     "FWA (BC), receiving"      = QA_PAL[["bc_dark"]],
                     "WBD (US)"                 = QA_PAL[["us"]],
                     "WBD (US), draining north" = QA_PAL[["us_dark"]])
    fwa$survey <- factor(role, levels = names(survey_cols))
    terra::plot(terra::vect(fwa), "survey", type = "classes",
                col = survey_cols,
                border = "white", lwd = 0.3, axes = TRUE,
                mar = c(4.6, 3.1, 3.1, 9.5),
                main = "Cross-border routing: US units draining into the FWA",
                xlim = c(fb[["xmin"]] - pad, fb[["xmax"]] + pad),
                ylim = c(fb[["ymin"]] - pad, fb[["ymax"]] + pad))

    if (nrow(cross) > 0) {
      # Find where the pair touches. st_nearest_points returns a segment whose
      # far end sits on the receiving watershed, and two touching shapes give a
      # zero-length segment on their shared edge. It falls back to the receiving
      # centre only if the pair is apart, the failure this figure looks for.
      exit <- t(vapply(seq_len(nrow(cross)), function(k) {
        a <- sf::st_geometry(fwa)[match(cross$id[k], fwa$HYBAS_ID)]
        b <- sf::st_geometry(fwa)[match(cross$nd[k], fwa$HYBAS_ID)]
        p <- suppressWarnings(sf::st_coordinates(sf::st_nearest_points(a, b)))
        if (nrow(p) >= 2) p[2, c("X", "Y")]
        else centroids[as.character(cross$nd[k]), c("X", "Y")]
      }, numeric(2)))

      graphics::segments(centroids[as.character(cross$id), "X"],
                         centroids[as.character(cross$id), "Y"],
                         exit[, 1], exit[, 2],
                         col = QA_PAL[["hot"]], lwd = 2.5)
      graphics::points(exit[, 1], exit[, 2], pch = 21, cex = 0.9,
                       col = QA_PAL[["hot"]], bg = "white", lwd = 1.6)
      graphics::mtext(paste0(nrow(cross), " cross-border edge(s); the arrow ",
                             "runs to the point where the unit meets its ",
                             "receiving FWA watershed"),
                      side = 1, line = 3.2, cex = 0.75,
                      col = QA_PAL[["slate"]])
    }
    # Crop rather than rely on par(xpd). terra keeps the map's aspect ratio, so
    # the coordinate window is wider than the box it draws, and a boundary
    # running past the frame gets drawn across the legend. Take the boundary
    # first and crop that, since cropping the polygon would turn the frame's own
    # edges into part of the outline.
    aoi_win <- suppressWarnings(sf::st_crop(
      sf::st_boundary(sf::st_geometry(aoi_disp)),
      c(xmin = fb[["xmin"]] - pad, xmax = fb[["xmax"]] + pad,
        ymin = fb[["ymin"]] - pad, ymax = fb[["ymax"]] + pad)))
    if (length(aoi_win) > 0)
      plot(aoi_win, add = TRUE, border = QA_PAL[["ink"]], lwd = 1.5, lty = 2)
  })
}

message("✓ 02_subbasins_fwa.R: wrote FWA sub-basins + topology + upstream AOI (",
        nrow(fwa), " assessment watersheds; ", sum(fwa$is_sink), " terminal sinks incl. ",
        sum(fwa$is_lake_barrier %in% TRUE), " lake barriers)")