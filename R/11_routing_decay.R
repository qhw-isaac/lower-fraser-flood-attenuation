# ==============================================================================
# 11_routing_decay.R — Distance-decayed downstream demand per sub-basin
# ------------------------------------------------------------------------------
# Mirrors Duarte's `scripts_OSF/03_sba/04_upstream_selected_basins.R` and
# `05_downstream_beneficiaries_mod_distance_decay_onlyfloodplains.R`.
#
# For every sub-basin j (provider candidate), sum the demand of every
# *downstream* basin i, weighted by the flow distance from j → i:
#
#       TDA_j = Σ_{i ∈ downstream(j)} demand_i × 0.5^(d_ji / HALFLIFE_KM)
#
# stopping at sinks (LAKE ≥ 1 with SUB_AREA ≥ LAKE_SINK_KM2, or ENDO == 2)
# and at MAX_FLOW_DIST_KM. Implementation: bottom-up traversal of the
# directed graph from sinks upward; each node's TDA is computed in O(deg)
# from its NEXT_DOWN's TDA — same shape as Duarte's loop, but vectorised
# via topological order from igraph.
#
# Inputs (data/processed/):
#   02_subbasins.gpkg, 02_topology.csv
#   10_demand_subbasin.gpkg
#
# Outputs (data/processed/):
#   11_tda_subbasin.gpkg
#       columns:
#         HYBAS_ID, in_downstream_aoi, in_outcome_aoi, admin_district,
#         own_demand, tda_both, tda_built, tda_crops,
#         tda_crops_w (regional), tda_total_w (regional headline)
# ==============================================================================

source(here::here("R", "00_setup.R"))

if (!requireNamespace("igraph", quietly = TRUE)) stop("install igraph")

sb   <- sf::st_read(file.path(paths()$processed, "02_subbasins.gpkg"), quiet = TRUE)
topo <- readr::read_csv(file.path(paths()$processed, "02_topology.csv"),
                        show_col_types = FALSE)
dem  <- sf::st_read(file.path(paths()$processed, "10_demand_subbasin.gpkg"),
                    quiet = TRUE) |> sf::st_drop_geometry()

# Patch NEXT_DOWN at sinks to 0 so the graph terminates correctly.
sink_ids <- sb$HYBAS_ID[sb$is_sink]
topo$ds_id[topo$focal_id %in% sink_ids] <- 0

# Belt-and-suspenders: drop any edge whose ds_id isn't a known vertex —
# 02_subbasins.R already rewrites dangling NEXT_DOWN to 0, but be defensive
sb_ids <- sb$HYBAS_ID

# Build the directed graph (focal → downstream), edge weight = reach length.
edges <- topo |>
  dplyr::filter(ds_id != 0, !is.na(ds_id), ds_id %in% sb_ids) |>
  dplyr::transmute(from = as.character(focal_id),
                   to   = as.character(ds_id),
                   weight = reach_km)
g <- igraph::graph_from_data_frame(
  edges, directed = TRUE,
  vertices = data.frame(name = as.character(sb_ids))
)

# Topological order (sinks last) → process upstream-first.
order_ids <- igraph::topo_sort(g, mode = "out")
order_ids <- as.numeric(igraph::V(g)$name[order_ids])

# Pre-index demand vectors for fast lookup.
demand_cols <- c("ben_both", "ben_built", "ben_crops",
                 "ben_crops_w", "total_demand_w")

idx <- match(sb$HYBAS_ID, dem$HYBAS_ID)
demand_mat <- as.matrix(dem[idx, demand_cols, drop = FALSE])
demand_mat[is.na(demand_mat)] <- 0
rownames(demand_mat) <- as.character(sb$HYBAS_ID)

ds_id   <- setNames(topo$ds_id,   as.character(topo$focal_id))
reach_km <- setNames(topo$reach_km, as.character(topo$focal_id))

# Bottom-up TDA accumulation. order_ids is "sources first, sinks last"
# `tda_partial`: contribution to TDA from each downstream node, *as seen at*
# this node — built up by halving along each reach.
tda <- demand_mat            # initialize: own demand sits at distance 0
walk <- rev(order_ids)
for (id in walk) {
  k <- as.character(id)
  if (!k %in% rownames(tda)) next
  ds <- ds_id[[k]]
  if (is.null(ds) || is.na(ds) || ds == 0) next
  reach <- reach_km[[k]]
  if (is.null(reach) || is.na(reach) || reach < 0) reach <- 0
  decay <- 0.5 ^ (reach / HALFLIFE_KM)
  ds_k  <- as.character(ds)
  if (ds_k %in% rownames(tda)) {
    # Note: own demand of `id` should NOT be decayed (it's already at j=id),
    # but its *downstream contribution* should be: pull DS's TDA, decay it,
    # add it to id's accumulator.
    tda[k, ] <- tda[k, ] + decay * tda[ds_k, ]
  }
}

# Pack into the sub-basin layer.
sb$own_demand   <- demand_mat[, "ben_both"]
sb$tda_both     <- tda[, "ben_both"]
sb$tda_built    <- tda[, "ben_built"]
sb$tda_crops    <- tda[, "ben_crops"]
sb$tda_crops_w  <- tda[, "ben_crops_w"]
sb$tda_total_w  <- tda[, "total_demand_w"]

sf::st_write(sb, file.path(paths()$processed, "11_tda_subbasin.gpkg"),
             delete_dsn = TRUE, quiet = TRUE)

# ---- QA preview --------------------------------------------------------------
# Distance-decayed downstream demand: should be highest in upstream sub-basins
# whose flow paths feed dense MVRD/FVRD demand pixels, decaying with distance.
qa_png("11_tda_total_w.png", function() {
  op <- graphics::par(mar = c(2, 2, 3, 1))
  on.exit(graphics::par(op), add = TRUE)
  plot(sb[sb$tda_total_w > 0, "tda_total_w"],
       main = "Total downstream demand (decayed; tda_total_w)",
       border = "grey60", lwd = 0.2, key.pos = 4,
       logz = TRUE, reset = FALSE)
  downstream_aoi <- read_aoi("downstream")
  plot(sf::st_geometry(downstream_aoi), add = TRUE, border = "black", lwd = 1.2)
})

message("✓ 11_routing_decay.R — wrote TDA columns for ", nrow(sb), " sub-basins")