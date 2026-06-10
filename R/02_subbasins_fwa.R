# ==============================================================================
# 02_subbasins_fwa.R — Sub-basin polygons + topology + upstream AOI, FWA variant
# ------------------------------------------------------------------------------
# Scenario / granularity variant of 02_subbasins.R. Where 02_subbasins.R uses
# HydroBASINS L12 (the finest *global* unit, with ready-made NEXT_DOWN/DIST_SINK
# connectivity), this script uses the BC Freshwater Atlas **Assessment
# Watersheds** — ~3x finer in this study area (789 vs a few hundred polygons),
# which is what makes the interactive map more granular.
#
# FWA polygons do NOT ship with a NEXT_DOWN field, so the core of this script is
# *deriving* the downstream topology from the FWA watershed-code system, then
# producing the SAME output contract as 02_subbasins.R so 03–13 run unchanged.
#
# Drop-in compatible: writes 02_subbasins.gpkg, 02_topology.csv and
# 01_aoi_upstream.gpkg, with the analysis-unit id carried in a column literally
# named HYBAS_ID (= FWA WATERSHED_FEATURE_ID) so every downstream join still
# works. Run THIS script for the FWA scenario; run 02_subbasins.R to restore
# HydroBASINS. Both overwrite the same three files, so whichever you ran last is
# the active sub-basin layer for 03+.
#
# ---- How the FWA watershed code encodes connectivity -------------------------
# FWA_WATERSHED_CODE = "100-190442-244975-...-000000": a base drainage (100)
# followed by up to 20 six-digit junction segments, zero-padded. Read
# left→right it walks DOWNSTREAM→UPSTREAM. The non-zero segments are the
# "address" of the stream a watershed sits on; a shorter address that is a
# prefix of a longer one is that stream's PARENT (it flows into the prefix).
# LOCAL_WATERSHED_CODE adds one more segment giving the watershed's route-
# measure (position) along its own stream — smaller measure = further
# downstream (verified: measure falls as Shreve magnitude rises).
#
# We use the codes only for flow DIRECTION and require the immediate downstream
# watershed to be a physical neighbour, which keeps the topology hydrologically
# clean and turns clip-boundary cases into honest domain outlets.
#
# Inputs:
#   data_path("fwa_assessment_watersheds")          (DataBC WFS extract)
#   data/processed/01_aoi_outcome.gpkg, 01_aoi_downstream.gpkg
#
# Outputs (data/processed/):
#   02_subbasins.gpkg     sub-basin polygons + tags (HYBAS_ID = FWA id)
#   02_topology.csv       focal_id, ds_id, reach_km   (watershed-code-derived)
#   01_aoi_upstream.gpkg  hydrologically refined upstream AOI
#
# SCOPE — why we pull our own extract here (do not pre-clip to the upstream AOI):
#   The downstream topology is a directed flow graph, and the ≤ MAX_FLOW_DIST_KM
#   trim keeps only watersheds that can *reach the Lower Mainland through that
#   graph*. If the extract is clipped to a ragged outline (e.g. the HydroBASINS-
#   derived upstream AOI), intermediate mainstem watersheds go missing, the
#   chain dead-ends at the gap (next_down = 0), and everything above it is
#   orphaned (flow_dist = Inf) and wrongly dropped — the bug that shrank the
#   domain to ~the Lower Mainland. Pulling a *continuous* buffer around the
#   downstream AOI keeps the drainage network unbroken; the trim then removes
#   only genuinely-distant or non-contributing watersheds. Validated: this
#   recovers 28,305 km² of contributing area vs HydroBASINS' 27,829 km².
#
# KNOWN LIMITATIONS (flagged for production hardening):
#   • reach_km is a Hack's-law estimate of main-channel length, not a measured
#     stream length — swap in summed FWA_STREAM_NETWORKS lengths per watershed
#     for exact distances. (Validated to match HydroBASINS' contributing area.)
#   • is_sink only marks domain outlets; lake/reservoir routing breaks (which
#     02_subbasins.R gets from HydroBASINS LAKE/ENDO) are not represented.
# ==============================================================================

source(here::here("R", "00_setup.R"))

if (!requireNamespace("igraph", quietly = TRUE)) stop("install igraph")

downstream_aoi <- read_aoi("downstream")
outcome_aoi    <- read_aoi("outcome")

# ---- 1. FWA assessment watersheds, scoped to a continuous contributing buffer -
# Pull from the DataBC WFS scoped to the downstream AOI grown by the full flow-
# distance cap (+ edge buffer) — the same preclip 02_subbasins.R applies to
# HydroBASINS. The result is cached to the registry path; if the WFS is
# unreachable we fall back to that cache. Set FLOOD_FWA_REFRESH= to force a fresh
# pull (e.g. after the downstream AOI changes).
fwa_path  <- data_path("fwa_assessment_watersheds")
fwa_buf   <- sf::st_buffer(sf::st_union(downstream_aoi),
                           MAX_FLOW_DIST_KM * 1000 + EDGE_BUFFER_M) |> sf::st_geometry()
need_pull <- !file.exists(fwa_path) || nzchar(Sys.getenv("FLOOD_FWA_REFRESH"))

fwa <- NULL
if (need_pull && requireNamespace("bcdata", quietly = TRUE)) {
  fwa <- tryCatch({
    message("  · pulling FWA assessment watersheds from DataBC WFS ",
            "(downstream AOI + ", MAX_FLOW_DIST_KM, " km buffer)…")
    # attach bcdata so the spatial predicate INTERSECTS() and the bcdc filter
    # method resolve inside the (non-standard-evaluated) query expression.
    suppressPackageStartupMessages(library(bcdata))
    q <- bcdc_query_geodata("WHSE_BASEMAPPING.FWA_ASSESSMENT_WATERSHEDS_POLY") |>
      filter(INTERSECTS(fwa_buf)) |>
      collect()
    q <- q[, setdiff(names(q), "SE_ANNO_CAD_DATA")] |> sf::st_transform(PROJECT_CRS)
    dir.create(dirname(fwa_path), showWarnings = FALSE, recursive = TRUE)
    sf::st_write(q, fwa_path, delete_dsn = TRUE, quiet = TRUE)   # refresh cache
    q
  }, error = function(e) {
    message("  ! WFS pull failed (", conditionMessage(e), ")"); NULL
  })
}
if (is.null(fwa)) {
  if (!file.exists(fwa_path))
    stop("no FWA extract: WFS pull failed and no cache at ", fwa_path)
  fwa <- sf::st_read(fwa_path, quiet = TRUE) |> sf::st_transform(PROJECT_CRS)
}

# ---- 2. Parse each watershed code into stream address + route-measure --------
# nz(): the non-zero FWA-code segments = a watershed's stream address.
nz <- function(code) { s <- strsplit(code, "-", fixed = TRUE)[[1]]; s[s != "000000"] }

base <- fwa |>
  sf::st_drop_geometry() |>
  dplyr::transmute(id = WATERSHED_FEATURE_ID,
                   magnitude = WATERSHED_MAGNITUDE,
                   area_km2  = AREA_HA / 100,
                   fwa = FWA_WATERSHED_CODE, loc = LOCAL_WATERSHED_CODE)

segs  <- lapply(base$fwa, nz)          # list-column of stream addresses
depth <- lengths(segs)                 # how many junctions from the ocean
# measure = the LOCAL segment one past the stream's depth = position on stream
measure <- as.numeric(mapply(\(l, d) strsplit(l, "-", fixed = TRUE)[[1]][d + 1],
                             base$loc, depth))
measure[is.na(measure)] <- 0           # most-downstream watershed: measure 0

ix <- setNames(seq_len(nrow(base)), base$id)   # id -> row, for O(1) lookups

# ---- 3. Pairwise "does a flow into b?" predicate -----------------------------
# b is downstream of a when ALL hold:
#   (i)   b's stream address is a prefix of a's (b is on a's stream or an
#         ancestor stream it ultimately drains to) — otherwise different branch;
#   (ii)  magnitude(b) >= magnitude(a) — flow only accumulates downstream, which
#         also rejects small lateral catchments that merely share a's code;
#   (iii) b sits at or below the point a's flow reaches b's stream:
#           same stream  -> measure(b) < measure(a)
#           ancestor     -> measure(b) <= the confluence segment of a's address.
is_downstream <- function(a, b) {
  ia <- ix[[as.character(a)]]; ib <- ix[[as.character(b)]]
  sa <- segs[[ia]];           sb <- segs[[ib]]
  if (length(sb) > length(sa)) return(FALSE)
  if (!all(sb == sa[seq_along(sb)])) return(FALSE)
  if (base$magnitude[ib] < base$magnitude[ia]) return(FALSE)
  if (length(sb) == length(sa)) {
    measure[ib] < measure[ia]
  } else {
    conf <- as.numeric(sa[length(sb) + 1])   # measure where a's line joins b's stream
    measure[ib] <= conf
  }
}

# ---- 4. NEXT_DOWN = adjacent downstream neighbour carrying the most flow ------
# Water leaves a watershed into exactly one physically adjacent neighbour, so we
# only consider neighbours (shared boundary). Among those that flow the right
# way, the real downstream is the trunk = the one with the largest magnitude.
# No adjacent downstream neighbour ⇒ the watershed drains out of the modelled
# area ⇒ outlet (NEXT_DOWN = 0).
touch <- sf::st_touches(fwa)
pairs <- data.frame(a = rep(base$id, lengths(touch)),
                    b = base$id[unlist(touch)])
pairs$flows <- mapply(is_downstream, pairs$a, pairs$b)
pairs$b_mag <- base$magnitude[ix[as.character(pairs$b)]]

next_down <- pairs |>
  dplyr::filter(flows) |>
  dplyr::group_by(a) |>
  dplyr::slice_max(b_mag, n = 1, with_ties = FALSE) |>
  dplyr::ungroup() |>
  dplyr::select(id = a, next_down = b)

base <- base |>
  dplyr::left_join(next_down, by = "id") |>
  dplyr::mutate(next_down = dplyr::coalesce(next_down, 0))

# reach length crossed when leaving a watershed: Hack's law main-channel length
# L(km) ≈ 1.4 · area_km2^0.6. Cumulative down a chain ≈ along-channel flow
# distance — the same network-following distance 11_routing_decay.R decays over.
base$reach_km <- 1.4 * base$area_km2^0.6

# ---- 5. Trim to ≤ MAX_FLOW_DIST_KM upstream of the downstream AOI -------------
# Same idea as 02_subbasins.R: keep only watersheds whose water reaches the
# Lower Mainland within MAX_FLOW_DIST_KM of along-network distance. Build the
# directed flow graph (focal → next_down, weight = reach_km), then measure each
# watershed's shortest downstream distance to any outlet (a watershed touching
# the downstream AOI) and drop those beyond the cap.
fwa$NEXT_DOWN <- base$next_down[match(fwa$WATERSHED_FEATURE_ID, base$id)]
fwa$reach_km  <- base$reach_km[match(fwa$WATERSHED_FEATURE_ID, base$id)]

ids   <- base$id
edges <- base |>
  dplyr::filter(next_down != 0, next_down %in% ids) |>
  dplyr::transmute(from = as.character(id), to = as.character(next_down),
                   weight = reach_km)
g <- igraph::graph_from_data_frame(edges, vertices = data.frame(name = as.character(ids)))

outlet_idx <- as.logical(sf::st_intersects(fwa, downstream_aoi, sparse = FALSE)[, 1])
outlet_ids <- fwa$WATERSHED_FEATURE_ID[outlet_idx]

dmat <- igraph::distances(g, to = igraph::V(g)[name %in% as.character(outlet_ids)],
                          mode = "out")
flow_dist <- apply(dmat, 1, min, na.rm = TRUE)
base$flow_dist_km <- flow_dist[match(base$id, as.numeric(igraph::V(g)$name))]

keep <- is.finite(base$flow_dist_km) & base$flow_dist_km <= MAX_FLOW_DIST_KM
fwa  <- fwa[match(base$id[keep], fwa$WATERSHED_FEATURE_ID), ]
base <- base[keep, ]

# next_down references to dropped watersheds become domain outlets (0)
dangling <- !base$next_down %in% c(0, base$id)
if (any(dangling)) {
  message("  · rewrote ", sum(dangling), " out-of-domain NEXT_DOWN references to 0")
  base$next_down[dangling] <- 0
}
fwa$NEXT_DOWN <- base$next_down[match(fwa$WATERSHED_FEATURE_ID, base$id)]

# ---- 6. AOI flags + sink flag + admin district -------------------------------
# Column names mirror 02_subbasins.R exactly so 10–14 need no changes.
fwa$HYBAS_ID         <- fwa$WATERSHED_FEATURE_ID            # analysis-unit id
fwa$SUB_AREA         <- fwa$AREA_HA / 100                   # km²
fwa$flow_dist_km     <- base$flow_dist_km[match(fwa$WATERSHED_FEATURE_ID, base$id)]
fwa$in_downstream_aoi <- as.logical(sf::st_intersects(sf::st_centroid(fwa), downstream_aoi, sparse = FALSE)[, 1])
fwa$in_outcome_aoi    <- as.logical(sf::st_intersects(sf::st_centroid(fwa), outcome_aoi,    sparse = FALSE)[, 1])
# No FWA lake/endorheic attributes here, so a sink is simply a domain outlet.
fwa$is_sink <- fwa$NEXT_DOWN == 0

if (requireNamespace("bcmaps", quietly = TRUE)) {
  rds <- bcmaps::regional_districts() |> sf::st_transform(PROJECT_CRS)
  fwa$admin_district <- "Other"
  for (nm in c("Metro Vancouver Regional District",
               "Fraser Valley Regional District",
               "Squamish-Lillooet Regional District")) {
    poly <- dplyr::filter(rds, ADMIN_AREA_NAME == nm)
    if (nrow(poly) == 1) {
      hit <- as.logical(sf::st_intersects(sf::st_centroid(fwa), poly, sparse = FALSE)[, 1])
      fwa$admin_district[hit] <- gsub(" Regional District$", "", nm)
    }
  }
}

# ---- 7. Persist sub-basin layer ----------------------------------------------
keep_cols <- c("HYBAS_ID", "NEXT_DOWN", "SUB_AREA", "flow_dist_km",
               "in_downstream_aoi", "in_outcome_aoi", "is_sink", "admin_district",
               "WATERSHED_FEATURE_ID", "WATERSHED_ORDER", "WATERSHED_MAGNITUDE",
               "FWA_WATERSHED_CODE", "GNIS_NAME_1")
sf::st_write(fwa[, intersect(keep_cols, names(fwa))],
             file.path(paths()$processed, "02_subbasins.gpkg"),
             delete_dsn = TRUE, quiet = TRUE)

# ---- 8. Topology CSV (focal → next_down with reach length km) ----------------
topo <- base |>
  dplyr::transmute(focal_id = id, ds_id = next_down, reach_km = reach_km)
readr::write_csv(topo, file.path(paths()$processed, "02_topology.csv"))

# ---- 9. Refined upstream AOI -------------------------------------------------
refined <- fwa |>
  sf::st_union() |>
  sf::st_buffer(EDGE_BUFFER_M) |>
  sf::st_sf(role = "upstream", crs = PROJECT_CRS)
sf::st_write(refined, file.path(paths()$processed, "01_aoi_upstream.gpkg"),
             delete_dsn = TRUE, quiet = TRUE)

# ---- QA preview --------------------------------------------------------------
qa_png("02_subbasins_fwa.png", panel_w = 1200, function() {
  op <- graphics::par(mar = c(2, 2, 4, 1), oma = c(0, 0, 1, 0))
  on.exit(graphics::par(op), add = TRUE)
  plot(fwa["flow_dist_km"],
       main = paste0("FWA assessment watersheds by flow distance to downstream AOI (",
                     nrow(fwa), " polygons)"),
       border = "grey40", lwd = 0.2, key.pos = 4, reset = FALSE)
  plot(sf::st_geometry(downstream_aoi), add = TRUE, border = "black", lwd = 1.2)
  plot(sf::st_geometry(outcome_aoi),    add = TRUE, border = "red",   lwd = 1.2)
})

message("✓ 02_subbasins_fwa.R — wrote FWA sub-basins + topology + upstream AOI (",
        nrow(fwa), " assessment watersheds; ", sum(fwa$is_sink), " outlets)")
