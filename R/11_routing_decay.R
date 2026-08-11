# ==============================================================================
# 11_routing_decay.R: Distance-decayed downstream demand per sub-basin
# ------------------------------------------------------------------------------
# For each sub-basin j, sum the demand of every downstream basin i, weighted by
# the flow distance from j to i. Demand decays by half every HALFLIFE_KM of
# travel (decay = 0.5 ^ (reach_km / HALFLIFE_KM)), so a sub-basin is credited
# most for the value at risk in the floodplains nearest below it.
#
# Inputs (data/processed/):
#   02_subbasins.gpkg, 02_topology.csv, 10_demand_subbasin.gpkg
#
# Outputs (data/processed/):
#   11_tda_subbasin.gpkg
# ==============================================================================

source(here::here("R", "00_setup.R"))

# ---- 1. Load sub-basins, topology, demand ------------------------------------
sb <- sf::st_read(file.path(paths()$processed, "02_subbasins.gpkg"), quiet = TRUE)
topo <- readr::read_csv(file.path(paths()$processed, "02_topology.csv"),
                        show_col_types = FALSE)
dem <- sf::st_read(file.path(paths()$processed, "10_demand_subbasin.gpkg"),
                    quiet = TRUE) |> sf::st_drop_geometry()

# ---- 2. Build the directed flow graph ----------------------------------------
# force sinks to ds_id 0 so the graph terminates
sink_ids <- sb$HYBAS_ID[sb$is_sink]
topo$ds_id[topo$focal_id %in% sink_ids] <- 0

sb_ids <- sb$HYBAS_ID

# keep only edges whose ds_id is a known vertex
edges <- topo |>
  dplyr::filter(ds_id != 0, !is.na(ds_id), ds_id %in% sb_ids) |>
  dplyr::transmute(from = as.character(focal_id),
                   to = as.character(ds_id),
                   weight = reach_km)
g <- igraph::graph_from_data_frame(
  edges, directed = TRUE,
  vertices = data.frame(name = as.character(sb_ids))
)

# ---- 3. Accumulate distance-decayed downstream demand ------------------------
if (!igraph::is_dag(g)) stop("flow graph is cyclic; check 02_topology.csv")
order_ids <- igraph::topo_sort(g, mode = "out")
order_ids <- as.numeric(igraph::V(g)$name[order_ids])

demand_cols <- c("ben_both", "ben_built", "ben_crops",
                 "ben_crops_w", "ben_facilities", "total_demand_w")

idx <- match(sb$HYBAS_ID, dem$HYBAS_ID)
demand_mat <- as.matrix(dem[idx, demand_cols, drop = FALSE])
demand_mat[is.na(demand_mat)] <- 0
rownames(demand_mat) <- as.character(sb$HYBAS_ID)

ds_id <- setNames(topo$ds_id, as.character(topo$focal_id))
reach_km <- setNames(topo$reach_km, as.character(topo$focal_id))

tda <- demand_mat
walk <- rev(order_ids)
for (id in walk) {
  k <- as.character(id)
  if (!k %in% rownames(tda)) next
  ds <- ds_id[k]
  if (is.na(ds) || ds == 0) next
  reach <- reach_km[k]
  if (is.na(reach) || reach < 0) reach <- 0
  decay <- 0.5 ^ (reach / HALFLIFE_KM)
  ds_k <- as.character(ds)
  if (ds_k %in% rownames(tda)) {
    tda[k, ] <- tda[k, ] + decay * tda[ds_k, ]
  }
}

# ---- 4. Attach columns and write ---------------------------------------------
sb$own_demand <- demand_mat[, "ben_both"]
sb$tda_both <- tda[, "ben_both"]
sb$tda_built <- tda[, "ben_built"]
sb$tda_crops <- tda[, "ben_crops"]
sb$tda_crops_w <- tda[, "ben_crops_w"]
sb$tda_facilities <- tda[, "ben_facilities"]
sb$tda_total_w <- tda[, "total_demand_w"]

sf::st_write(sb, file.path(paths()$processed, "11_tda_subbasin.gpkg"),
             delete_dsn = TRUE, quiet = TRUE)

# ---- QA preview --------------------------------------------------------------
# routed demand should peak in sub-basins whose flow paths feed dense MVRD/FVRD
# exposure and fade with distance. 
qa_png("11_tda_total_w.png", function() {
  terra::plot(terra::vect(sb), "tda_total_w", type = "interval",
              breaks = c(0, 0.01, 0.1, 1, 10, Inf),
              col = grDevices::hcl.colors(5, "YlGnBu", rev = TRUE),
              border = "grey70", lwd = 0.1, axes = TRUE, main = "")
  plot(sf::st_geometry(read_aoi("downstream")), add = TRUE,
       border = "black", lwd = 1.4)
})

message("✓ 11_routing_decay.R: wrote TDA columns for ", nrow(sb), " sub-basins")