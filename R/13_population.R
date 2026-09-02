# ==============================================================================
# 13_population.R: Benefiting people downstream of priority sub-basins
# ------------------------------------------------------------------------------
# Population unit = whichever geography POP_UNIT names in 00_setup.R: 2021
# Census dissemination areas (DAs, the default) or dissemination blocks (DBs).
# Units are addressed by POP_UID throughout; DAUID is carried for labelling.
#
# Steps:
#   1. Fetch BC DAs (geometry + population), clip to the downstream AOI
#   2. Intersect each DA with sub-basins
#   3. Per DA, split population by floodplain share
#   4. Per RI interval, tally DAs in basins downstream of any interval member
#
# Inputs (data/processed/):
#   07_floodplain.tif, 02_subbasins.gpkg, 02_topology.csv,
#   12_realised_benefit_<scenario>.gpkg
#
# Inputs (via read_das() in 00_setup.R):
#   da_boundary_2021 (lda_000b21a_e.shp)
#   da_population_2021 (2021_92-151_X.csv)
#
# Outputs (output/tables/):
#   13_benefiting_population_<scenario>.csv
#   (ri_interval, n_popctr, direct_pop, indirect_pop)
# ==============================================================================

source(here::here("R", "00_setup.R"))

# ---- 1. Population by dissemination area -------------------------------------
# assembled by read_das() in 00_setup.R
pop <- read_das()
pop$pop_area_km2 <- as.numeric(sf::st_area(pop)) * 1e-6
message("  · ", nrow(pop), " dissemination areas in downstream AOI")

# ---- 2. Intersect with sub-basins, build flood overlap -----------------------
sb <- sf::st_read(file.path(paths()$processed, "02_subbasins.gpkg"), quiet = TRUE)
fp <- terra::rast(file.path(paths()$processed, "07_floodplain.tif"))

inter <- sf::st_intersection(pop, sb |> dplyr::select(HYBAS_ID))
inter$area_km2 <- as.numeric(sf::st_area(inter)) * 1e-6

# flood area per piece = flood-pixel count x pixel area, then split each DA's
# population by the share of its area sitting in the floodplain
px_km2 <- (terra::xres(fp) * terra::yres(fp)) * 1e-6
inter$flood_km2 <- px_km2 * exactextractr::exact_extract(fp == 1, inter, "sum")
inter <- inter |>
  dplyr::group_by(POP_UID) |>
  dplyr::mutate(pop_in_flood = COUNT_TOTAL * flood_km2 / sum(area_km2)) |>
  dplyr::ungroup()

# ---- 3. Downstream tally per RI interval -------------------------------------
topo <- readr::read_csv(file.path(paths()$processed, "02_topology.csv"),
                        show_col_types = FALSE)
sink_ids <- sb$HYBAS_ID[sb$is_sink]
topo$ds_id[topo$focal_id %in% sink_ids] <- 0

g <- igraph::graph_from_data_frame(
  topo |> dplyr::filter(ds_id != 0, !is.na(ds_id)) |>
    dplyr::transmute(from = as.character(focal_id), to = as.character(ds_id)),
  vertices = data.frame(name = as.character(sb$HYBAS_ID))
)

scenarios <- list.files(paths()$processed,
                        pattern = "^12_realised_benefit_.*\\.gpkg$",
                        full.names = TRUE)

for (rb_path in scenarios) {
  scen <- sub("^12_realised_benefit_(.*)\\.gpkg$", "\\1", basename(rb_path))
  rb <- sf::st_read(rb_path, quiet = TRUE)

  out_rows <- list()
  for (lab in unique(stats::na.omit(rb$ri_interval))) {
    seeds <- rb$HYBAS_ID[rb$ri_interval == lab & !is.na(rb$ri_interval)]
    if (length(seeds) == 0) next
    reach_ids <- unique(unlist(lapply(
      as.character(seeds),
      function(s) igraph::V(g)$name[igraph::subcomponent(g, s, mode = "out")])))
    reach_ids <- as.numeric(reach_ids)

    members <- inter[inter$HYBAS_ID %in% reach_ids, ]
    direct <- sum(members$pop_in_flood, na.rm = TRUE)
    centres <- unique(members$POP_UID)
    total_in_centres <- sum(pop$COUNT_TOTAL[pop$POP_UID %in% centres], na.rm = TRUE)
    indirect <- total_in_centres - direct

    out_rows[[length(out_rows) + 1]] <- dplyr::tibble(
      ri_interval = as.character(lab),
      n_popctr = length(centres),
      direct_pop = direct,
      indirect_pop = indirect
    )
  }

  # highest percentile tier first; the non-percentile "lake_buffer" row last
  # (a plain desc() would sort it above "95_100" alphabetically)
  out <- dplyr::bind_rows(out_rows) |>
    dplyr::arrange(ri_interval == "lake_buffer", dplyr::desc(ri_interval))
  f <- file.path(paths()$tables, glue::glue("13_benefiting_population_{scen}.csv"))
  readr::write_csv(out, f)
  message("  ✓ '", scen, "' → ", basename(f))
}

# ---- QA preview --------------------------------------------------------------
downstream_aoi <- read_aoi("downstream")

# per-unit flood summary (one row per DA, or per DB when POP_UNIT = "DB").
da_summary <- inter |>
  sf::st_drop_geometry() |>
  dplyr::group_by(POP_UID) |>
  dplyr::summarise(
    pop_in_flood = sum(pop_in_flood, na.rm = TRUE),
    flood_km2 = sum(flood_km2, na.rm = TRUE),
    area_km2 = sum(area_km2, na.rm = TRUE),
    .groups = "drop"
  )

# choropleth layer: each unit shaded by its own flood exposure %.
pop_map <- pop |>
  dplyr::left_join(da_summary |> dplyr::select(POP_UID, flood_km2, area_km2),
                   by = "POP_UID")
pop_map$flood_pct <- ifelse(!is.na(pop_map$area_km2) & pop_map$area_km2 > 0,
                            100 * pop_map$flood_km2 / pop_map$area_km2, 0)
pop_map$flood_pct[is.na(pop_map$flood_pct)] <- 0

pop_map$flood_pct <- pmin(pmax(pop_map$flood_pct, 0), 100)

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

# roll DA floodplain population up to municipalities, so the previews read by
# community even though the analysis runs at DA resolution
muni <- bcmaps::municipalities(ask = FALSE) |> sf::st_transform(PROJECT_CRS) |>
  sf::st_filter(downstream_aoi)
muni$nm <- tidy_muni_names(muni$ADMIN_AREA_NAME)

hit <- sf::st_intersects(sf::st_centroid(sf::st_geometry(pop)), muni)
da_nm <- vapply(hit, function(h) if (length(h)) muni$nm[h[1]] else NA_character_,
                 character(1))
centres_ranked <- da_summary |>
  dplyr::mutate(nm = da_nm[match(POP_UID, pop$POP_UID)]) |>
  dplyr::filter(!is.na(nm)) |>
  dplyr::group_by(nm) |>
  dplyr::summarise(pop_in_flood = sum(pop_in_flood, na.rm = TRUE),
                   flood_km2 = sum(flood_km2, na.rm = TRUE),
                   area_km2 = sum(area_km2, na.rm = TRUE),
                   .groups = "drop") |>
  dplyr::mutate(flood_pct = ifelse(area_km2 > 0, 100 * flood_km2 / area_km2, 0))

lc <- sf::st_coordinates(sf::st_centroid(sf::st_geometry(muni)))

view <- if (nrow(lc) > 0) {
  c(xmin = min(lc[, 1]), ymin = min(lc[, 2]),
    xmax = max(lc[, 1]), ymax = max(lc[, 2]))
} else sf::st_bbox(pop_map)
view_padx <- 0.08 * (view["xmax"] - view["xmin"])
view_pady <- 0.10 * (view["ymax"] - view["ymin"])
view_pad_e <- 0.015 * (view["xmax"] - view["xmin"])

land_fill <- "#f5f3ed"
water_fill <- "#d6e4ef"

# land base = DAs dissolved; land-only footprint (no ocean)
land_poly <- sf::st_make_valid(sf::st_union(sf::st_make_valid(pop_map)))

# floodplain outline for geographic context
fp_outline <- sf::st_as_sf(
  terra::as.polygons(terra::ifel(fp == 1, 1L, NA))) |>
  sf::st_make_valid() |>
  sf::st_union() |>
  sf::st_make_valid()

muni_layer <- ggplot2::geom_sf(data = muni, fill = NA, colour = "white",
                               linewidth = 0.35, alpha = 0.8)

outline_cols <- c("Modelled floodplain" = "#3182bd")

p_map <- ggplot2::ggplot() +
  ggplot2::geom_sf(data = land_poly, fill = land_fill, colour = NA) +
  ggplot2::geom_sf(data = pop_map, ggplot2::aes(fill = flood_pct),
                   colour = NA, na.rm = TRUE) +
  ggplot2::geom_sf(data = fp_outline,
                   ggplot2::aes(colour = "Modelled floodplain"),
                   fill = NA, linewidth = 0.25, alpha = 0.5) +
  muni_layer +
  ggplot2::scale_colour_manual(values = outline_cols, name = NULL,
                               guide = ggplot2::guide_legend(
                                 order = 2, override.aes = list(
                                   linewidth = 0.7, alpha = 1))) +
  ggplot2::scale_fill_gradientn(
    colours = c(land_fill, "#fef0d9", "#fdcc8a", "#fc8d59",
                "#e34a33", "#b30000"),
    values = c(0, 0.02, 0.10, 0.30, 0.60, 1),
    limits = c(0, 100),
    na.value = land_fill,
    name = "Flood exposure\n(% of area)",
    guide = ggplot2::guide_colourbar(
      barwidth = grid::unit(0.8, "lines"),
      barheight = grid::unit(8, "lines"),
      order = 1,
      ticks.linewidth = 0.3)) +
  ggplot2::coord_sf(
    xlim = c(view["xmin"] - view_padx, view["xmax"] + view_pad_e),
    ylim = c(view["ymin"] - view_pady, view["ymax"] + view_pady),
    datum = sf::st_crs(4326), expand = FALSE) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(
    axis.title = ggplot2::element_blank(),
    axis.text = ggplot2::element_text(size = 7, colour = "grey40"),
    panel.grid = ggplot2::element_line(colour = "grey85", linewidth = 0.25,
                                       linetype = "dotted"),
    panel.background = ggplot2::element_rect(fill = water_fill, colour = NA),
    legend.position = "right",
    legend.title = ggplot2::element_text(size = 8.5, face = "bold"),
    legend.text = ggplot2::element_text(size = 7.5),
    plot.margin = ggplot2::margin(8, 14, 6, 8)
  )

ggplot2::ggsave(file.path(qa_dir(), "13_population_map.png"),
                p_map, width = 13, height = 5.5, dpi = 180, bg = "white")
message("  \u00b7 qa preview: ", file.path(qa_dir(), "13_population_map.png"))

message("✓ 13_population.R: completed ", length(scenarios), " scenario(s)")