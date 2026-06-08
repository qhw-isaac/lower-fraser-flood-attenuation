# ==============================================================================
# 13_population.R — Benefiting people downstream of priority sub-basins
# ------------------------------------------------------------------------------
#   1. Dissolve StatCan Population Centres by DGUID; sum LANDAREA; join
#      Census 2021 'Population, 2021' (`COUNT_TOTAL`).
#   2. Intersect each pop-centre polygon with sub-basin polygons → area_fraction.
#   3. Per pop-centre × sub-basin row, compute:
#        pop_in_basin   = COUNT_TOTAL × area_fraction
#        flood_km2      = pop-centre ∩ floodplain ∩ sub-basin
#        flood_fraction = flood_km2 / pop-centre area
#        pop_in_flood   = COUNT_TOTAL × flood_fraction      (≈ "directly affected")
#   4. For each RI interval (12_realised_benefit_*.gpkg), find every sub-basin
#      downstream of an interval member (graph traversal, capped at MAX_FLOW_DIST_KM):
#        Direct   = Σ pop_in_flood for pop-centres in those basins
#        Indirect = (total population of those pop-centres) − Direct
#        n_pop_centres = unique DGUID count
#
# Inputs (data/processed/):
#   07_floodplain.tif, 02_subbasins.gpkg, 02_topology.csv
#   12_realised_benefit_<scenario>_all.gpkg
#
# Inputs (data/raw/):
#   data_path("pop_centres_2021")         — boundaries (shapefile)
#   data_path("pop_centres_2021_census")  — Census Profile CSV (Population, 2021)
#
# Outputs (output/tables/):
#   13_benefiting_population_<scenario>.csv
#       columns: ri_interval, n_popctr, direct_pop, indirect_pop
# ==============================================================================

source(here::here("R", "00_setup.R"))

# ---- 1. Dissolved + population-joined pop-centre layer -----------------------
pop_path <- data_path("pop_centres_2021")
pop_full <- sf::st_read(pop_path, quiet = TRUE) |> sf::st_transform(PROJECT_CRS)

if (!"DGUID" %in% names(pop_full)) {
  stop("Pop-centres layer missing DGUID column.")
}

pop <- pop_full |>
  dplyr::group_by(DGUID) |>
  # `summarise.sf()` auto-unions geometries per group. The earlier
  # `summarise(geometry = sf::st_union(sf::st_geometry(pop_full)))` was a
  # well-known dplyr×sf footgun: it unions *every* geometry once and assigns
  # the same blob to each group.
  dplyr::summarise(.groups = "drop")

if ("COUNT_TOTAL" %in% names(pop_full)) {
  pop_meta <- pop_full |>
    sf::st_drop_geometry() |>
    dplyr::select(DGUID, COUNT_TOTAL) |>
    dplyr::distinct(DGUID, .keep_all = TRUE)
} else {
  # StatCan boundary file + Census Profile CSV (DGUID must be read as character).
  census <- readr::read_csv(
    data_path("pop_centres_2021_census"),
    col_select = c(DGUID, CHARACTERISTIC_ID, C1_COUNT_TOTAL),
    col_types = readr::cols(
      DGUID = readr::col_character(),
      CHARACTERISTIC_ID = readr::col_double(),
      C1_COUNT_TOTAL = readr::col_double()
    ),
    show_col_types = FALSE
  )
  pop_meta <- census |>
    dplyr::filter(CHARACTERISTIC_ID == 1L) |>
    dplyr::transmute(DGUID, COUNT_TOTAL = as.numeric(C1_COUNT_TOTAL)) |>
    dplyr::distinct(DGUID, .keep_all = TRUE)
}
pop <- dplyr::left_join(pop, pop_meta, by = "DGUID")
if (any(is.na(pop$COUNT_TOTAL))) {
  n_miss <- sum(is.na(pop$COUNT_TOTAL))
  warning(n_miss, " pop-centre(s) missing COUNT_TOTAL after census join")
}
pop$pop_area_km2 <- as.numeric(sf::st_area(pop)) * 1e-6

# ---- 2. Intersect with sub-basins, build flood overlap -----------------------
sb <- sf::st_read(file.path(paths()$processed, "02_subbasins.gpkg"), quiet = TRUE)
fp <- terra::rast(file.path(paths()$processed, "07_floodplain.tif"))

inter <- sf::st_intersection(pop, sb |> dplyr::select(HYBAS_ID))
inter$area_km2 <- as.numeric(sf::st_area(inter)) * 1e-6
inter <- inter |>
  dplyr::group_by(DGUID) |>
  dplyr::mutate(area_fraction = area_km2 / sum(area_km2)) |>
  dplyr::ungroup()
inter$pop_in_basin <- inter$COUNT_TOTAL * inter$area_fraction

# Flood overlap area per intersection polygon.
fp_km2 <- terra::ifel(fp == 1, (terra::xres(fp) * terra::yres(fp)) * 1e-6, 0)
inter$flood_km2 <- exactextractr::exact_extract(fp_km2, inter, "sum")
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
qa_tables <- list()

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
  qa_tables[[scen]] <- out
  message("  ✓ '", scen, "' → ", basename(f))
}

# ---- QA preview --------------------------------------------------------------
downstream_aoi <- read_aoi("downstream")

# Per-centre summary: name, total pop, pop in floodplain, flood exposure %
pop_summary <- inter |>
  sf::st_drop_geometry() |>
  dplyr::group_by(DGUID) |>
  dplyr::summarise(
    pop_in_flood = sum(pop_in_flood, na.rm = TRUE),
    flood_km2    = sum(flood_km2,    na.rm = TRUE),
    area_km2     = sum(area_km2,     na.rm = TRUE),
    .groups = "drop"
  )
pcnames <- pop_full |>
  sf::st_drop_geometry() |>
  dplyr::distinct(DGUID, .keep_all = TRUE) |>
  dplyr::select(DGUID, PCNAME)

pop_map <- sf::st_filter(pop, downstream_aoi) |>
  dplyr::left_join(pop_summary, by = "DGUID") |>
  dplyr::left_join(pcnames,     by = "DGUID")
pop_map$pop_in_flood[is.na(pop_map$pop_in_flood)] <- 0
pop_map$flood_pct <- ifelse(pop_map$area_km2 > 0,
                            pop_map$flood_km2 / pop_map$area_km2 * 100, 0)

# (a) Map: pop centres shaded by flood exposure %, with community labels
fp_poly <- sf::st_as_sf(terra::as.polygons(terra::ifel(fp == 1, 1L, NA))) |>
  sf::st_make_valid()

label_pts <- sf::st_centroid(pop_map)
label_coords <- sf::st_coordinates(label_pts)
label_df <- data.frame(
  PCNAME = label_pts$PCNAME,
  x      = label_coords[, 1],
  y      = label_coords[, 2]
)

p_map <- ggplot2::ggplot() +
  ggplot2::geom_sf(data = downstream_aoi, fill = "grey96", colour = "grey50",
                   linewidth = 0.3) +
  ggplot2::geom_sf(data = fp_poly, fill = "#b3d4f7", colour = NA, alpha = 0.5) +
  ggplot2::geom_sf(data = pop_map, ggplot2::aes(fill = flood_pct),
                   colour = "grey30", linewidth = 0.15) +
  ggrepel::geom_label_repel(
    data = label_df,
    ggplot2::aes(x = x, y = y, label = PCNAME),
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
                                name = "Flood exposure\n(% of area)",
                                limits = c(0, NA)) +
  ggplot2::labs(title = "Population centre flood exposure (MVRD \u222a FVRD)",
                subtitle = "Background: 100-yr floodplain (Mohanty)",
                x = "Longitude", y = "Latitude") +
  ggplot2::coord_sf(datum = sf::st_crs(4326)) +
  ggplot2::theme_minimal() +
  ggplot2::theme(legend.position = "right")

ggplot2::ggsave(file.path(qa_dir(), "13_population_map.png"),
                p_map, width = 10, height = 7, dpi = 150)
message("  \u00b7 qa preview: ", file.path(qa_dir(), "13_population_map.png"))

# (b) Paired bars: people in floodplain (left) + flood exposure % (right)
centres_ranked <- pop_map |>
  sf::st_drop_geometry() |>
  dplyr::filter(pop_in_flood > 0) |>
  dplyr::arrange(pop_in_flood) |>
  dplyr::mutate(PCNAME = factor(PCNAME, levels = PCNAME))

if (nrow(centres_ranked) > 0) {
  p_abs <- ggplot2::ggplot(centres_ranked,
                           ggplot2::aes(x = pop_in_flood, y = PCNAME)) +
    ggplot2::geom_col(fill = "#e6550d", width = 0.6) +
    ggplot2::geom_text(
      ggplot2::aes(label = scales::label_comma(accuracy = 1)(round(pop_in_flood))),
      hjust = -0.1, size = 3
    ) +
    ggplot2::scale_x_continuous(labels = scales::label_comma(),
                                expand = ggplot2::expansion(mult = c(0, 0.25))) +
    ggplot2::labs(title = "People in 100-yr floodplain",
                  x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())

  p_pct <- ggplot2::ggplot(centres_ranked,
                           ggplot2::aes(x = flood_pct, y = PCNAME)) +
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
