# ==============================================================================
# 02_subbasins.R — Sub-basin polygons + topology + upstream AOI
# ------------------------------------------------------------------------------
# Walks the HydroBASINS L12 NEXT_DOWN graph upstream from the downstream AOI
# and keeps every sub-basin within ≤ MAX_FLOW_DIST_KM flow distance.
#
# Runs before the raster scripts (03–07) so they inherit the hydrologically
# refined upstream AOI rather than a regional-district placeholder, which
# would clip basins draining into the Lower Mainland from outside MVRD ∪ FVRD
# ∪ SLRD (e.g. Fraser/Thompson basins through TNRD).
#
# Sub-basin tags:
#   HYBAS_ID, NEXT_DOWN, SUB_AREA, DIST_SINK, LAKE, ENDO   (HydroBASINS)
#   in_downstream_aoi, in_outcome_aoi                      (centroid tests)
#   is_sink                                                (LAKE ≥ 1 & SUB_AREA
#                                                           ≥ LAKE_SINK_KM2
#                                                           OR ENDO == 2)
#   admin_district                                         (MVRD/FVRD/SLRD/Other)
#
# Inputs:
#   data_path("hydrobasins_l12")
#   data/processed/01_aoi_outcome.gpkg, 01_aoi_downstream.gpkg
#
# Outputs (data/processed/):
#   02_subbasins.gpkg     sub-basin polygons + tags
#   02_topology.csv       focal_id, ds_id, reach_km   (DIST_SINK-derived)
#   01_aoi_upstream.gpkg  hydrologically refined upstream AOI
# ==============================================================================

source(here::here("R", "00_setup.R"))

# ---- 1. Load HydroBASINS L12 -------------------------------------------------
hb <- sf::st_read(data_path("hydrobasins_l12"), quiet = TRUE) |>
  sf::st_transform(PROJECT_CRS)

hb$LAKE[is.na(hb$LAKE)] <- 0L
hb$ENDO[is.na(hb$ENDO)] <- 0L

# preclip to save resources: anything within MAX_FLOW_DIST_KM (+ edge buffer)
# of the downstream AOI is a candidate upstream provider.
downstream_aoi <- read_aoi("downstream")
prelim <- sf::st_buffer(downstream_aoi, MAX_FLOW_DIST_KM * 1000 + EDGE_BUFFER_M)
hb <- hb[as.logical(sf::st_intersects(hb, prelim, sparse = FALSE)), ]

# ---- 2. Trim to ≤ MAX_FLOW_DIST_KM upstream of downstream AOI ----------------
outlet_idx <- as.logical(sf::st_intersects(hb, downstream_aoi, sparse = FALSE))
outlet_ids <- hb$HYBAS_ID[outlet_idx]

# Drop edges whose downstream endpoint is outside our vertex set: those
# basins (the MVRD outlets that drain to coastal sub-basins, plus any
# interior basin whose NEXT_DOWN happened to fall outside the prelim
# buffer) are *natural sinks* in our subdomain — they have nowhere left
# to flow within the modeled area. Without this filter,
# igraph::graph_from_data_frame() refuses to build because edges
# reference vertices that don't exist.
hb_ids <- hb$HYBAS_ID

edges <- hb |>
  sf::st_drop_geometry() |>
  dplyr::transmute(from = HYBAS_ID, to = NEXT_DOWN, dist_sink = DIST_SINK) |>
  dplyr::filter(to != 0, to %in% hb_ids) |>
  dplyr::left_join(
    hb |> sf::st_drop_geometry() |> dplyr::select(HYBAS_ID, ds_dist = DIST_SINK),
    by = c("to" = "HYBAS_ID")
  ) |>
  dplyr::mutate(weight = pmax(dist_sink - ds_dist, 0)) |>
  dplyr::select(from, to, weight)

g <- igraph::graph_from_data_frame(edges, vertices = data.frame(name = hb_ids))

dmat <- igraph::distances(
  g,
  to = igraph::V(g)[name %in% as.character(outlet_ids)],
  mode = "out"
)
flow_dist <- apply(dmat, 1, min, na.rm = TRUE)

hb$flow_dist_km <- flow_dist[match(hb$HYBAS_ID, as.numeric(igraph::V(g)$name))]
hb <- hb[is.finite(hb$flow_dist_km) & hb$flow_dist_km <= MAX_FLOW_DIST_KM, ]

# After the flow-distance trim, some basins' NEXT_DOWN may now point to
# basins that survived the prelim buffer but failed the MAX_FLOW_DIST_KM
# trim. Rewrite those dangling references to 0 so the topology CSV
# (and 11_routing_decay.R's graph build) treat them as our domain's outlets.
dangling <- !hb$NEXT_DOWN %in% c(0L, hb$HYBAS_ID)
n_dangling <- sum(dangling)
if (n_dangling > 0L) {
  message("  · rewrote ", n_dangling, " out-of-domain NEXT_DOWN references to 0 ",
          "(natural drains-out-of-modeled-area; treated as sinks)")
  hb$NEXT_DOWN[dangling] <- 0L
}

# ---- 3. AOI flags + sink flag + admin district -------------------------------
outcome_aoi <- read_aoi("outcome")
hb$in_downstream_aoi <- as.logical(sf::st_intersects(sf::st_centroid(hb), downstream_aoi, sparse = FALSE)[, 1])
hb$in_outcome_aoi    <- as.logical(sf::st_intersects(sf::st_centroid(hb), outcome_aoi,    sparse = FALSE)[, 1])
hb$is_sink <- (hb$LAKE >= 1 & hb$SUB_AREA >= LAKE_SINK_KM2) | hb$ENDO == 2

if (requireNamespace("bcmaps", quietly = TRUE)) {
  rds <- bcmaps::regional_districts() |> sf::st_transform(PROJECT_CRS)
  hb$admin_district <- "Other"
  for (nm in c("Metro Vancouver Regional District",
               "Fraser Valley Regional District",
               "Squamish-Lillooet Regional District")) {
    poly <- dplyr::filter(rds, ADMIN_AREA_NAME == nm)
    if (nrow(poly) == 1) {
      hit <- as.logical(sf::st_intersects(sf::st_centroid(hb), poly, sparse = FALSE)[, 1])
      hb$admin_district[hit] <- gsub(" Regional District$", "", nm)
    }
  }
}

# ---- 4. Persist sub-basin layer ----------------------------------------------
sf::st_write(hb, file.path(paths()$processed, "02_subbasins.gpkg"),
             delete_dsn = TRUE, quiet = TRUE)

# ---- 5. Topology CSV (focal → next_down with reach length km) ----------------
topo <- hb |>
  sf::st_drop_geometry() |>
  dplyr::select(HYBAS_ID, NEXT_DOWN, DIST_SINK) |>
  dplyr::rename(focal_id = HYBAS_ID, ds_id = NEXT_DOWN) |>
  dplyr::left_join(
    hb |> sf::st_drop_geometry() |> dplyr::select(HYBAS_ID, ds_dist = DIST_SINK),
    by = c("ds_id" = "HYBAS_ID")
  ) |>
  dplyr::mutate(reach_km = pmax(DIST_SINK - ds_dist, 0)) |>
  dplyr::select(focal_id, ds_id, reach_km)
readr::write_csv(topo, file.path(paths()$processed, "02_topology.csv"))

# ---- 6. Refined upstream AOI -------------------------------------------------
refined <- hb |>
  sf::st_union() |>
  sf::st_buffer(EDGE_BUFFER_M) |>
  sf::st_sf(role = "upstream", crs = PROJECT_CRS)
sf::st_write(refined, file.path(paths()$processed, "01_aoi_upstream.gpkg"),
             delete_dsn = TRUE, quiet = TRUE)

# ---- QA preview --------------------------------------------------------------
qa_png("02_subbasins.png", panel_w = 1200, function() {
  op <- graphics::par(mar = c(2, 2, 4, 1), oma = c(0, 0, 1, 0))
  on.exit(graphics::par(op), add = TRUE)
  plot(hb["flow_dist_km"],
       main = paste0("Sub-basins by flow distance to downstream AOI outlet (",
                     nrow(hb), " polygons)"),
       border = "grey40", lwd = 0.2, key.pos = 4,
       reset = FALSE)
  plot(sf::st_geometry(downstream_aoi), add = TRUE, border = "black", lwd = 1.2)
  plot(sf::st_geometry(outcome_aoi),    add = TRUE, border = "red",   lwd = 1.2)
})

message("✓ 02_subbasins.R — wrote 02_subbasins.gpkg + 01_aoi_upstream.gpkg (",
        nrow(hb), " HydroBASINS L12 polygons)")