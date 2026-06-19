# ==============================================================================
# 13_population.R — benefiting people downstream of priority sub-basins
# ------------------------------------------------------------------------------
# population unit = 2021 census dissemination areas
#
#   1. fetch BC DAs (geometry + population), clip to the downstream AOI
#   2. intersect each DA with sub-basins
#   3. per DA, split population by floodplain share
#   4. per RI interval, tally DAs in basins downstream of any interval member
#
# inputs (data/processed/): 07_floodplain.tif, 02_subbasins.gpkg, 02_topology.csv,
#   12_realised_benefit_<scenario>_all.gpkg
#
# inputs (data/raw/population/):
#   da_boundary_2021 (lda_000b21a_e.shp), 
#   da_population_2021 (2021_92-151_X.csv) cached as da_2021_lowermainland.gpkg
#
# outputs (output/tables/): 
#   13_benefiting_population_<scenario>.csv
#   (ri_interval, n_popctr, direct_pop, indirect_pop)
# ==============================================================================

source(here::here("R", "00_setup.R"))

# ---- 1. Fine-geography population: 2021 Census Dissemination Areas -----------
# two StatCan-direct downloads joined on DAUID, clipped to the downstream AOI
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
  g_id <- grep("DAUID", names(gaf), value = TRUE, ignore.case = TRUE)[1]
  g_pop <- grep("DBPOP|DB_POP|POP_?2021|POPULATION", names(gaf),
                value = TRUE, ignore.case = TRUE)[1]
  if (is.na(g_id) || is.na(g_pop))
    stop("could not find DAUID / population columns in da_population_2021 (cols: ",
         paste(names(gaf), collapse = ", "), ")")
  pop_tbl <- gaf |>
    dplyr::transmute(DGUID = as.character(.data[[g_id]]),
                     pop = suppressWarnings(as.numeric(.data[[g_pop]]))) |>
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

# ---- 2. intersect with sub-basins, build flood overlap -----------------------
sb <- sf::st_read(file.path(paths()$processed, "02_subbasins.gpkg"), quiet = TRUE)
fp <- terra::rast(file.path(paths()$processed, "07_floodplain.tif"))

inter <- sf::st_intersection(pop, sb |> dplyr::select(HYBAS_ID))
inter$area_km2 <- as.numeric(sf::st_area(inter)) * 1e-6

# flood overlap per piece = flood-pixel count * pixel area (km2); then split each
# DA's population by the share of its area in the floodplain.
px_km2 <- (terra::xres(fp) * terra::yres(fp)) * 1e-6
inter$flood_km2 <- px_km2 * exactextractr::exact_extract(fp == 1, inter, "sum")
inter <- inter |>
  dplyr::group_by(DGUID) |>
  dplyr::mutate(pop_in_flood = COUNT_TOTAL * flood_km2 / sum(area_km2)) |>
  dplyr::ungroup()

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
    # all basins reachable downstream from any seed (mode = "out").
    reach <- igraph::subcomponent(g, v = as.character(seeds), mode = "out")
    reach_ids <- as.numeric(igraph::V(g)$name[reach])

    members <- inter[inter$HYBAS_ID %in% reach_ids, ]
    direct <- sum(members$pop_in_flood, na.rm = TRUE)
    centres <- unique(members$DGUID)
    total_in_centres <- sum(pop$COUNT_TOTAL[pop$DGUID %in% centres], na.rm = TRUE)
    indirect <- total_in_centres - direct

    out_rows[[length(out_rows) + 1]] <- dplyr::tibble(
      ri_interval = as.character(lab),
      n_popctr = length(centres),
      direct_pop = direct,
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

# per-DA flood summary (one row per DA).
da_summary <- inter |>
  sf::st_drop_geometry() |>
  dplyr::group_by(DGUID) |>
  dplyr::summarise(
    pop_in_flood = sum(pop_in_flood, na.rm = TRUE),
    flood_km2 = sum(flood_km2, na.rm = TRUE),
    area_km2 = sum(area_km2, na.rm = TRUE),
    .groups = "drop"
  )

# choropleth layer: each DA shaded by its own flood exposure %.
pop_map <- pop |>
  dplyr::left_join(da_summary |> dplyr::select(DGUID, flood_km2, area_km2),
                   by = "DGUID")
pop_map$flood_pct <- ifelse(!is.na(pop_map$area_km2) & pop_map$area_km2 > 0,
                            100 * pop_map$flood_km2 / pop_map$area_km2, 0)
pop_map$flood_pct[is.na(pop_map$flood_pct)] <- 0

# clean community label from an official name: strip "The" + admin prefixes
# ("City/District/Township/... of") repeatedly until stable (some names nest
# them), then tag shared base names (City vs Township of Langley) to disambiguate.
tidy_muni_names <- function(raw_nm) {
  base <- tools::toTitleCase(tolower(raw_nm))
  repeat {
    prev <- base
    base <- sub("^[Tt]he ", "", base)
    base <- sub("^(City|District|Town|Township|Village|Corporation|Resort|Municipality)( (City|District|Municipality))? [Oo]f ",
                "", base)
    if (identical(base, prev)) break
  }
  base <- sub(" Island Municipality$", " Island", base)

  kind <- rep("", length(raw_nm))
  kind[grepl("Township", raw_nm, ignore.case = TRUE)] <- "Twp"
  kind[grepl("\\bCity\\b", raw_nm, ignore.case = TRUE)] <- "City"
  kind[grepl("District", raw_nm, ignore.case = TRUE) & kind == ""] <- "Dist"
  dup <- base %in% base[duplicated(base)]
  ifelse(dup & nzchar(kind), paste0(base, " (", kind, ")"), base)
}

# aggregate DA floodplain population up to named municipalities, so the previews
# still read by community even though the analysis runs at DA resolution.
muni <- NULL
if (requireNamespace("bcmaps", quietly = TRUE)) {
  muni <- bcmaps::municipalities() |> sf::st_transform(PROJECT_CRS) |>
    sf::st_filter(downstream_aoi)
  muni$nm <- tidy_muni_names(muni$ADMIN_AREA_NAME)
}
centres_ranked <- NULL
if (!is.null(muni) && nrow(muni) > 0) {
  hit <- sf::st_intersects(sf::st_centroid(sf::st_geometry(pop)), muni)
  da_nm <- vapply(hit, function(h) if (length(h)) muni$nm[h[1]] else NA_character_,
                   character(1))
  centres_ranked <- da_summary |>
    dplyr::mutate(nm = da_nm[match(DGUID, pop$DGUID)]) |>
    dplyr::filter(!is.na(nm)) |>
    dplyr::group_by(nm) |>
    dplyr::summarise(pop_in_flood = sum(pop_in_flood, na.rm = TRUE),
                     flood_km2 = sum(flood_km2, na.rm = TRUE),
                     area_km2 = sum(area_km2, na.rm = TRUE),
                     .groups = "drop") |>
    dplyr::mutate(flood_pct = ifelse(area_km2 > 0, 100 * flood_km2 / area_km2, 0))
}

# (a) map: DAs shaded by flood exposure %, with municipality labels
fp_poly <- sf::st_as_sf(terra::as.polygons(terra::ifel(fp == 1, 1L, NA))) |>
  sf::st_make_valid()

label_df <- if (!is.null(muni) && nrow(muni) > 0) {
  lc <- sf::st_coordinates(sf::st_centroid(sf::st_geometry(muni)))
  data.frame(nm = muni$nm, x = lc[, 1], y = lc[, 2])
} else data.frame(nm = character(0), x = numeric(0), y = numeric(0))

# frame on the municipal label points (Salish Sea -> Hope corridor); sparse
# upstream DAs up the Pitt/Stave/Harrison valleys fall outside this view.
view <- if (nrow(label_df) > 0) {
  c(xmin = min(label_df$x), ymin = min(label_df$y),
    xmax = max(label_df$x), ymax = max(label_df$y))
} else sf::st_bbox(fp_poly)
view_padx <- 0.08 * (view["xmax"] - view["xmin"])
view_pady <- 0.10 * (view["ymax"] - view["ymin"])
# tight east cap at Hope — FVRD continues up the Fraser Canyon, but those sparse
# canyon DAs are cropped from this view (display only; still in the analysis).
view_pad_e <- 0.015 * (view["xmax"] - view["xmin"])

# white municipal boundaries read as clean dividers over the choropleth.
muni_layer <- if (!is.null(muni) && nrow(muni) > 0) {
  ggplot2::geom_sf(data = muni, fill = NA, colour = "white", linewidth = 0.45)
} else NULL

# land fill = low end of the exposure ramp, so unexposed DAs blend into land and
# only flood-exposed areas read as colour; water is a bluer tone so the coast reads.
land_fill <- "#edeae2"
water_fill <- "#b7cfe0"

# land base = the DAs dissolved. the MVRD/FVRD boundary covers ocean, so using it
# as the mask would paint water as land; the DA footprint is land only, and
# dissolving fills sliver gaps between DAs.
land_poly <- sf::st_make_valid(sf::st_union(sf::st_make_valid(pop_map)))

p_map <- ggplot2::ggplot() +
  ggplot2::geom_sf(data = land_poly, fill = land_fill, colour = NA) +
  ggplot2::geom_sf(data = pop_map, ggplot2::aes(fill = flood_pct), colour = NA) +
  muni_layer +
  ggrepel::geom_label_repel(
    data = label_df,
    ggplot2::aes(x = x, y = y, label = nm),
    size = 2.5, colour = "#1f2d3a", fontface = "bold",
    fill = ggplot2::alpha("white", 0.82),
    label.size = 0, label.r = grid::unit(0.12, "lines"),
    label.padding = grid::unit(0.12, "lines"),
    segment.colour = "grey55", segment.size = 0.25,
    max.overlaps = 40, min.segment.length = 0,
    box.padding = 0.4, point.padding = 0.2,
    force = 3, force_pull = 0.5,
    seed = 42
  ) +
  ggplot2::scale_fill_gradientn(
    colours = c(land_fill, "#fee391", "#fe9929", "#e34a33", "#b10026"),
    values = scales::rescale(c(0, 3, 20, 50, 100)),
    limits = c(0, 100),
    name = "Flood\nexposure\n(% of DA area)") +
  ggplot2::labs(title = "Dissemination-area flood exposure (MVRD \u222a FVRD)",
                subtitle = "Census dissemination areas shaded by share within the 100-yr floodplain (Mohanty); white lines are municipal boundaries") +
  ggplot2::coord_sf(
    xlim = c(view["xmin"] - view_padx, view["xmax"] + view_pad_e),
    ylim = c(view["ymin"] - view_pady, view["ymax"] + view_pady),
    datum = sf::st_crs(4326), expand = FALSE) +
  ggplot2::theme_void(base_size = 12) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 14),
    plot.subtitle = ggplot2::element_text(colour = "grey35", size = 9,
                                          margin = ggplot2::margin(b = 6)),
    panel.background = ggplot2::element_rect(fill = water_fill, colour = NA),
    panel.border = ggplot2::element_rect(fill = NA, colour = "grey75", linewidth = 0.4),
    legend.position = "right",
    legend.title = ggplot2::element_text(size = 9),
    plot.margin = ggplot2::margin(10, 10, 8, 10)
  )

ggplot2::ggsave(file.path(qa_dir(), "13_population_map.png"),
                p_map, width = 11, height = 5.5, dpi = 150)
message("  \u00b7 qa preview: ", file.path(qa_dir(), "13_population_map.png"))

# (b) paired bars by municipality (from DA counts): people in floodplain (left)
#     + exposure % (right).
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
                                expand = ggplot2::expansion(mult = c(0, 0.30))) +
    ggplot2::scale_y_discrete(expand = ggplot2::expansion(add = 0.8)) +
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
    ggplot2::scale_x_continuous(limits = c(0, 100),
                                expand = ggplot2::expansion(mult = c(0, 0.08))) +
    ggplot2::scale_y_discrete(expand = ggplot2::expansion(add = 0.8)) +
    ggplot2::labs(title = "Flood exposure (% of area)",
                  x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank()
    )

  p_paired <- patchwork::wrap_plots(p_abs, p_pct, ncol = 2, widths = c(1.2, 1))

  ggplot2::ggsave(file.path(qa_dir(), "13_population_exposure.png"),
                  p_paired, width = 11, height = 6, dpi = 150)
  message("  \u00b7 qa preview: ", file.path(qa_dir(), "13_population_exposure.png"))
}

message("✓ 13_population.R — completed ", length(scenarios), " scenario(s)")