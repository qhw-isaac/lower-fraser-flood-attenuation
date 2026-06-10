# ==============================================================================
# 13_population.R — Benefiting people downstream of priority sub-basins
# ------------------------------------------------------------------------------
#   Population unit = 2021 Census Dissemination Areas (DAs, ~400-700 people),
#   NOT StatCan Population Centres. The single 2.4M-person "Vancouver" centre
#   made floodplain population a uniform-density areal estimate and inflated
#   "indirect" to the whole metro; DAs fix both by resolving exposure at
#   small-area scale. DGUID now holds a DA GeoUID, COUNT_TOTAL its population.
#
#   1. Fetch BC DAs (geometry + population) via cancensus; cache to data/raw;
#      keep those meeting the downstream AOI.
#   2. Intersect each DA polygon with sub-basin polygons → area_fraction.
#   3. Per DA × sub-basin row, compute:
#        pop_in_basin   = COUNT_TOTAL × area_fraction
#        flood_km2      = DA ∩ floodplain ∩ sub-basin
#        flood_fraction = flood_km2 / DA area
#        pop_in_flood   = COUNT_TOTAL × flood_fraction      (≈ "directly affected")
#   4. For each RI interval (12_realised_benefit_*.gpkg), find every sub-basin
#      downstream of an interval member (igraph subcomponent, mode = "out"):
#        Direct   = Σ pop_in_flood for DAs in those basins
#        Indirect = (total population of those DAs) − Direct
#        n_popctr = unique DA count
#      This mirrors Duarte's 05_realized_benefits/03 downstream tally. Their
#      traversal applied a 500 km distance cutoff (downstream_focalbasins_thredist);
#      we apply none here because the whole domain is already trimmed to
#      ≤ MAX_FLOW_DIST_KM in 02_subbasins.R, so reachable == in-range.
#
# Inputs (data/processed/):
#   07_floodplain.tif, 02_subbasins.gpkg, 02_topology.csv
#   12_realised_benefit_<scenario>_all.gpkg
#
# Inputs (data/raw/population/ — StatCan-direct, no API key):
#   da_boundary_2021   DA cartographic boundary shapefile (lda_000b21a_e.shp)
#   da_population_2021  Geographic Attribute File CSV (2021_92-151_X.csv)
#   → joined on DAUID, clipped to AOI, cached as da_2021_lowermainland.gpkg
#
# Outputs (output/tables/):
#   13_benefiting_population_<scenario>.csv
#       columns: ri_interval, n_popctr, direct_pop, indirect_pop
# ==============================================================================

source(here::here("R", "00_setup.R"))

# ---- 1. Fine-geography population: 2021 Census Dissemination Areas -----------
# Repointed from StatCan Population Centres to Dissemination Areas (DAs, the
# ~400-700-person small-area unit). The single 2.4M-person "Vancouver" centre
# made floodplain population a uniform-density areal estimate (Richmond/Delta
# undercounted) and inflated "indirect" to the whole metro. DAs resolve
# floodplain exposure at small-area scale and make the direct/indirect split
# meaningful.
#
# Source = two StatCan-direct downloads (no API key / quota), joined on DAUID,
# clipped to the downstream AOI, then cached to data/raw for fast reruns:
#   da_boundary_2021   — DA cartographic boundary shapefile (geometry + DAUID)
#   da_population_2021  — Geographic Attribute File: one row per dissemination
#                         BLOCK with its population; summed by DAUID → DA pop.
# Column names are detected (StatCan suffixes them bilingually). See
# data_sources.csv for the exact download URLs.
da_path <- file.path(paths()$raw, "population", "da_2021_lowermainland.gpkg")
if (!file.exists(da_path)) {
  dir.create(dirname(da_path), showWarnings = FALSE, recursive = TRUE)

  # geometry + DA id, clipped to the downstream AOI (the boundary file is national)
  bnd <- sf::st_read(data_path("da_boundary_2021", must_exist = TRUE), quiet = TRUE) |>
    sf::st_transform(PROJECT_CRS)
  id_col <- grep("^DAUID", names(bnd), value = TRUE, ignore.case = TRUE)[1]
  if (is.na(id_col)) stop("no DAUID column in DA boundary file (cols: ",
                          paste(names(bnd), collapse = ", "), ")")
  bnd$DGUID <- as.character(bnd[[id_col]])
  bnd <- bnd[as.logical(sf::st_intersects(
    bnd, sf::st_union(read_aoi("downstream")), sparse = FALSE)[, 1]), ]

  # population: detect DAUID + block-population columns, aggregate to DA level
  gaf <- readr::read_csv(data_path("da_population_2021", must_exist = TRUE),
                         show_col_types = FALSE, guess_max = 50000)
  g_id  <- grep("DAUID", names(gaf), value = TRUE, ignore.case = TRUE)[1]
  g_pop <- grep("DBPOP|DB_POP|POP_?2021|POPULATION", names(gaf),
                value = TRUE, ignore.case = TRUE)[1]
  if (is.na(g_id) || is.na(g_pop))
    stop("could not find DAUID / population columns in da_population_2021 (cols: ",
         paste(names(gaf), collapse = ", "), ")")
  pop_tbl <- gaf |>
    dplyr::transmute(DGUID = as.character(.data[[g_id]]),
                     pop   = suppressWarnings(as.numeric(.data[[g_pop]]))) |>
    dplyr::group_by(DGUID) |>
    dplyr::summarise(COUNT_TOTAL = sum(pop, na.rm = TRUE), .groups = "drop")

  da <- bnd |>
    dplyr::select(DGUID) |>
    dplyr::left_join(pop_tbl, by = "DGUID") |>
    sf::st_make_valid()
  sf::st_write(da, da_path, delete_dsn = TRUE, quiet = TRUE)
  message("  · built ", nrow(da), " dissemination areas (downstream AOI) → ", da_path)
}
pop <- sf::st_read(da_path, quiet = TRUE) |> sf::st_transform(PROJECT_CRS)
if (any(is.na(pop$COUNT_TOTAL))) {
  warning(sum(is.na(pop$COUNT_TOTAL)), " DA(s) missing population — set to 0")
  pop$COUNT_TOTAL[is.na(pop$COUNT_TOTAL)] <- 0
}
pop$pop_area_km2 <- as.numeric(sf::st_area(pop)) * 1e-6
message("  · ", nrow(pop), " dissemination areas in downstream AOI")

# ---- 2. Intersect with sub-basins, build flood overlap -----------------------
# Duarte's 03_population_centers/02 dropped pop-centre slivers falling in
# big-lake sub-basins (LAKE < 1 & SUB_AREA < lake_size) before splitting
# population by area. We don't replicate that filter: there are no
# routing-breaking lakes inside the Lower Mainland pop-centre footprint, so it
# would remove nothing. Revisit if the AOI ever extends to large interior lakes.
sb <- sf::st_read(file.path(paths()$processed, "02_subbasins.gpkg"), quiet = TRUE)
fp <- terra::rast(file.path(paths()$processed, "07_floodplain.tif"))

inter <- sf::st_intersection(pop, sb |> dplyr::select(HYBAS_ID))
inter$area_km2 <- as.numeric(sf::st_area(inter)) * 1e-6
inter <- inter |>
  dplyr::group_by(DGUID) |>
  dplyr::mutate(area_fraction = area_km2 / sum(area_km2)) |>
  dplyr::ungroup()
inter$pop_in_basin <- inter$COUNT_TOTAL * inter$area_fraction

# Flood overlap per intersection polygon = flood-pixel count × pixel area (km²).
px_km2 <- (terra::xres(fp) * terra::yres(fp)) * 1e-6
inter$flood_km2 <- px_km2 * exactextractr::exact_extract(fp == 1, inter, "sum")
inter <- inter |>
  dplyr::group_by(DGUID) |>
  dplyr::mutate(flood_fraction = flood_km2 / sum(area_km2)) |>
  dplyr::ungroup()
inter$pop_in_flood <- inter$COUNT_TOTAL * inter$flood_fraction

# ---- 3. Downstream tally per RI interval -------------------------------------
topo <- readr::read_csv(file.path(paths()$processed, "02_topology.csv"),
                        show_col_types = FALSE)
sink_ids <- sb$HYBAS_ID[sb$is_sink]
topo$ds_id[topo$focal_id %in% sink_ids] <- 0

if (!requireNamespace("igraph", quietly = TRUE)) stop("install igraph")
g <- igraph::graph_from_data_frame(
  topo |> dplyr::filter(ds_id != 0, !is.na(ds_id)) |>
    dplyr::transmute(from = as.character(focal_id), to = as.character(ds_id)),
  vertices = data.frame(name = as.character(sb$HYBAS_ID))
)

scenarios <- list.files(paths()$processed,
                        pattern = "^12_realised_benefit_.*_all\\.gpkg$",
                        full.names = TRUE)

for (rb_path in scenarios) {
  scen <- sub("^12_realised_benefit_(.*)_all\\.gpkg$", "\\1", basename(rb_path))
  rb <- sf::st_read(rb_path, quiet = TRUE)

  out_rows <- list()
  for (lab in unique(stats::na.omit(rb$ri_interval))) {
    seeds <- rb$HYBAS_ID[rb$ri_interval == lab & !is.na(rb$ri_interval)]
    if (length(seeds) == 0) next
    # All basins reachable downstream from any seed (mode = "out").
    reach <- igraph::subcomponent(g, v = as.character(seeds), mode = "out")
    reach_ids <- as.numeric(igraph::V(g)$name[reach])

    members <- inter[inter$HYBAS_ID %in% reach_ids, ]
    direct  <- sum(members$pop_in_flood, na.rm = TRUE)
    centres <- unique(members$DGUID)
    total_in_centres <- sum(pop$COUNT_TOTAL[pop$DGUID %in% centres], na.rm = TRUE)
    indirect <- total_in_centres - direct

    out_rows[[length(out_rows) + 1]] <- dplyr::tibble(
      ri_interval = as.character(lab),
      n_popctr    = length(centres),
      direct_pop  = direct,
      indirect_pop = indirect
    )
  }

  out <- dplyr::bind_rows(out_rows) |>
    dplyr::arrange(dplyr::desc(ri_interval))
  f <- file.path(paths()$tables, glue::glue("13_benefiting_population_{scen}.csv"))
  readr::write_csv(out, f)
  message("  ✓ '", scen, "' → ", basename(f))
}

# ---- QA preview --------------------------------------------------------------
downstream_aoi <- read_aoi("downstream")

# Per-DA flood summary (DAs are atomic; one row per DA).
da_summary <- inter |>
  sf::st_drop_geometry() |>
  dplyr::group_by(DGUID) |>
  dplyr::summarise(
    pop_in_flood = sum(pop_in_flood, na.rm = TRUE),
    flood_km2    = sum(flood_km2,    na.rm = TRUE),
    area_km2     = sum(area_km2,     na.rm = TRUE),
    .groups = "drop"
  )

# Choropleth layer: each DA shaded by its own flood exposure %.
pop_map <- pop |>
  dplyr::left_join(da_summary |> dplyr::select(DGUID, flood_km2, area_km2),
                   by = "DGUID")
pop_map$flood_pct <- ifelse(!is.na(pop_map$area_km2) & pop_map$area_km2 > 0,
                            100 * pop_map$flood_km2 / pop_map$area_km2, 0)
pop_map$flood_pct[is.na(pop_map$flood_pct)] <- 0

# Aggregate DA floodplain population to municipalities (named) so the previews
# still read by community even though the analysis now runs at DA resolution.
muni <- NULL
if (requireNamespace("bcmaps", quietly = TRUE)) {
  muni <- bcmaps::municipalities() |> sf::st_transform(PROJECT_CRS) |>
    sf::st_filter(downstream_aoi)
  muni$nm <- tools::toTitleCase(tolower(muni$ADMIN_AREA_NAME))
}
centres_ranked <- NULL
if (!is.null(muni) && nrow(muni) > 0) {
  hit    <- sf::st_intersects(sf::st_centroid(sf::st_geometry(pop)), muni)
  da_nm  <- vapply(hit, function(h) if (length(h)) muni$nm[h[1]] else NA_character_,
                   character(1))
  centres_ranked <- da_summary |>
    dplyr::mutate(nm = da_nm[match(DGUID, pop$DGUID)]) |>
    dplyr::filter(!is.na(nm)) |>
    dplyr::group_by(nm) |>
    dplyr::summarise(pop_in_flood = sum(pop_in_flood, na.rm = TRUE),
                     flood_km2    = sum(flood_km2,    na.rm = TRUE),
                     area_km2     = sum(area_km2,     na.rm = TRUE),
                     .groups = "drop") |>
    dplyr::mutate(flood_pct = ifelse(area_km2 > 0, 100 * flood_km2 / area_km2, 0))
}

# (a) Map: DAs shaded by flood exposure %, with municipality labels
fp_poly <- sf::st_as_sf(terra::as.polygons(terra::ifel(fp == 1, 1L, NA))) |>
  sf::st_make_valid()

label_df <- if (!is.null(muni) && nrow(muni) > 0) {
  lc <- sf::st_coordinates(sf::st_centroid(sf::st_geometry(muni)))
  data.frame(nm = muni$nm, x = lc[, 1], y = lc[, 2])
} else data.frame(nm = character(0), x = numeric(0), y = numeric(0))

p_map <- ggplot2::ggplot() +
  ggplot2::geom_sf(data = downstream_aoi, fill = "grey96", colour = "grey50",
                   linewidth = 0.3) +
  ggplot2::geom_sf(data = fp_poly, fill = "#b3d4f7", colour = NA, alpha = 0.5) +
  ggplot2::geom_sf(data = pop_map, ggplot2::aes(fill = flood_pct),
                   colour = NA) +
  ggrepel::geom_label_repel(
    data = label_df,
    ggplot2::aes(x = x, y = y, label = nm),
    size = 2.6, colour = "grey10",
    fill = ggplot2::alpha("white", 0.75),
    label.size = 0,
    segment.colour = "grey40", segment.size = 0.3,
    max.overlaps = 30, min.segment.length = 0,
    box.padding = 0.5, point.padding = 0.3,
    force = 2, force_pull = 0.5,
    seed = 42
  ) +
  ggplot2::scale_fill_viridis_c(option = "inferno", direction = -1,
                                name = "Flood exposure\n(% of DA area)",
                                limits = c(0, NA)) +
  ggplot2::labs(title = "Dissemination-area flood exposure (MVRD \u222a FVRD)",
                subtitle = "Background: 100-yr floodplain (Mohanty)",
                x = "Longitude", y = "Latitude") +
  ggplot2::coord_sf(datum = sf::st_crs(4326)) +
  ggplot2::theme_minimal() +
  ggplot2::theme(legend.position = "right")

ggplot2::ggsave(file.path(qa_dir(), "13_population_map.png"),
                p_map, width = 10, height = 7, dpi = 150)
message("  \u00b7 qa preview: ", file.path(qa_dir(), "13_population_map.png"))

# (b) Paired bars by municipality: people in floodplain (left) + exposure % (right),
#     aggregated up from the DA-resolution counts.
if (!is.null(centres_ranked) && any(centres_ranked$pop_in_flood > 0)) {
  centres_ranked <- centres_ranked |>
    dplyr::filter(pop_in_flood > 0) |>
    dplyr::slice_max(pop_in_flood, n = 20) |>
    dplyr::arrange(pop_in_flood) |>
    dplyr::mutate(nm = factor(nm, levels = nm))

  p_abs <- ggplot2::ggplot(centres_ranked,
                           ggplot2::aes(x = pop_in_flood, y = nm)) +
    ggplot2::geom_col(fill = "#e6550d", width = 0.6) +
    ggplot2::geom_text(
      ggplot2::aes(label = scales::label_comma(accuracy = 1)(round(pop_in_flood))),
      hjust = -0.1, size = 3
    ) +
    ggplot2::scale_x_continuous(labels = scales::label_comma(),
                                expand = ggplot2::expansion(mult = c(0, 0.25))) +
    ggplot2::labs(title = "People in 100-yr floodplain (from DA counts)",
                  x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())

  p_pct <- ggplot2::ggplot(centres_ranked,
                           ggplot2::aes(x = flood_pct, y = nm)) +
    ggplot2::geom_col(fill = "#756bb1", width = 0.6) +
    ggplot2::geom_text(
      ggplot2::aes(label = paste0(round(flood_pct, 1), "%")),
      hjust = -0.1, size = 3
    ) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.3))) +
    ggplot2::labs(title = "Flood exposure (% of area)",
                  x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank()
    )

  p_paired <- patchwork::wrap_plots(p_abs, p_pct, ncol = 2, widths = c(1.2, 1))

  ggplot2::ggsave(file.path(qa_dir(), "13_population_exposure.png"),
                  p_paired, width = 11, height = 5, dpi = 150)
  message("  \u00b7 qa preview: ", file.path(qa_dir(), "13_population_exposure.png"))
}

message("✓ 13_population.R — completed ", length(scenarios), " scenario(s)")
