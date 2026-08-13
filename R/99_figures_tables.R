# ==============================================================================
# 99_figures_tables.R: Publication-ready figures + summary table
# ------------------------------------------------------------------------------
# Reads finished artefacts and produces report-quality outputs. Every map shows
# the full upstream provider extent but anchors the Lower Mainland demand area
# (MVRD ∪ FVRD), so the reader can see what is being protected.
#
# These do not restate the interactive map. The map answers "who benefits from
# this watershed", which moves with the demand weights, the floodplain extent
# and the population vintage. These cover what holds still: what makes a
# watershed retain runoff, which ecosystems do it, and which places rank highly
# whichever storm is modelled.
#
# Outputs:
#   Fig 2a    PRR (mm) per pixel, demand area outlined
#   Fig 2b    RI percentile per sub-basin (two-panel: overview + demand zoom)
#   Fig 3     Retention against its physical drivers (precipitation, slope, cover)
#   Fig 4     Retention by ecosystem type: total volume and depth per unit area
#   Fig 5     Scenario robustness: rank agreement + the consistently high set
#   Table 1   RI interval summary × downstream beneficiaries
#   Table 2   Retention by ecosystem × scenario
#   Table 3   Watersheds in the top tenth under both scenarios
# ==============================================================================

source(here::here("R", "00_setup.R"))

# ==============================================================================
# 0. SCENARIO + COMMON LAYERS
# ==============================================================================

scen <- DEFAULT_SCENARIO

# readable scenario names, defined here because the first figure already needs
# them for its title (the robustness figures further down use them too)
scen_lab <- c(ar2021 = "November 2021 storm",
              wettest_month = "Wettest month (climatology)")
pretty_scen <- function(s) unname(ifelse(s %in% names(scen_lab), scen_lab[s], s))

rb_path <- file.path(paths()$processed,
                     glue::glue("12_realised_benefit_{scen}.gpkg"))
rb  <- sf::st_read(rb_path, quiet = TRUE)
prr <- terra::rast(file.path(paths()$processed, "runoff", scen, "09_prr_mm.tif"))
aoi <- aoi_display()

# RI interval factor levels (used in multiple figures)
ri_lvls <- c("0_25", "25_50", "50_60", "60_70",
             "70_80", "80_85", "85_90", "90_95", "95_100")
ri_labs <- c("0\u201325", "25\u201350", "50\u201360", "60\u201370",
             "70\u201380", "80\u201385", "85\u201390", "90\u201395", "95\u2013100")

# drop tiny disconnected raster specks, for display only
drop_small_patches <- function(r, min_px = 100L) {
  bin <- terra::ifel(!is.na(r), 1L, NA_integer_)
  ids <- terra::patches(bin, directions = 8)
  pf  <- terra::freq(ids)
  if (is.null(pf) || nrow(pf) == 0L) return(r)
  keep <- pf$value[pf$count >= min_px]
  if (length(keep) == 0L) return(r)
  terra::mask(r, terra::ifel(ids %in% keep, 1L, NA_integer_))
}

# ==============================================================================
# FIG 2a: POTENTIAL RUNOFF RETENTION (mm) PER PIXEL
# ==============================================================================
# every modelled pixel, with the demand-area boundary drawn over it so the
# upstream-provision to downstream-benefit relationship reads directly

# --- raster masks: model domain, water bodies, snow/ice ---------------------
sb_mask_path <- file.path(paths()$tmp, glue::glue("99_sb_mask_{scen}.tif"))
sb_mask <- terra::rasterize(terra::vect(rb), prr, field = 1,
                            filename = sb_mask_path, overwrite = TRUE)

src_cov <- terra::rast(file.path(paths()$processed, "03_lulc_source_coverage.tif"))

nalcms_path <- file.path(paths()$tmp, glue::glue("99_nalcms_{scen}.tif"))
nalcms <- terra::project(terra::rast(data_path("nalcms_2020")), prr,
                         method = "near", filename = nalcms_path,
                         overwrite = TRUE)

prec <- terra::rast(file.path(paths()$processed, "precip",
                             glue::glue("06_p_{scen}.tif")))

in_domain  <- !is.na(sb_mask) & !is.na(src_cov)
model_mask <- in_domain & !is.na(prec)
water      <- in_domain & nalcms == 18L
snow_ice   <- in_domain & nalcms == 19L

prr_display <- terra::mask(prr, model_mask, maskvalues = c(NA, FALSE))
prr_display <- terra::ifel(is.na(prr_display) & snow_ice, 0, prr_display)
water_overlay <- terra::ifel(water, 1, NA)

# --- Crop extent to the downstream AOI's southern edge ----------------------
down_bb <- sf::st_bbox(aoi)
prr_bb  <- terra::ext(prr_display)
data_ext <- terra::ext(prr_bb[1], prr_bb[2], down_bb["ymin"], prr_bb[4])

prr_plot   <- terra::crop(prr_display, data_ext)
water_plot <- terra::crop(water_overlay, data_ext)

# clean speckle on the smaller extent
display_mask <- terra::ifel(!is.na(prr_plot) | !is.na(water_plot),
                            1L, NA_integer_)
display_mask <- drop_small_patches(display_mask, min_px = 100L)
display_mask <- terra::sieve(display_mask, threshold = 25L, directions = 8)
prr_plot   <- terra::mask(prr_plot, display_mask)
water_plot <- terra::mask(water_plot, display_mask)

# --- Sub-basin outlines (only over the coloured footprint) ------------------
color_poly <- terra::as.polygons(display_mask, dissolve = TRUE, values = FALSE)
rb_lines <- sf::st_intersection(sf::st_crop(rb, data_ext),
                                sf::st_as_sf(color_poly))
if (nrow(rb_lines) > 0L)
  rb_lines <- rb_lines[as.numeric(sf::st_length(rb_lines)) > 120, ]

# --- Padded extent so data doesn't touch the axis frame ---------------------
pad <- 0.04
xw  <- data_ext[2] - data_ext[1]
yw  <- data_ext[4] - data_ext[3]
fig_ext <- terra::ext(data_ext[1] - pad * xw, data_ext[2] + pad * xw,
                      data_ext[3] - pad * yw, data_ext[4] + pad * yw)

# --- Draw the figure --------------------------------------------------------
png(file.path(paths()$figures, glue::glue("99_fig2_prr_{scen}.png")),
    width = 1600, height = 1250, res = 200)

# short title on the raster itself: the report supplies the full caption, but
# these PNGs are also read on their own, where unlabelled they say nothing
terra::plot(prr_plot,
            background = "white",
            mar = c(3.1, 3.1, 4.1, 7.1),
            main = glue::glue("Potential Runoff Retention (mm): ",
                              "{pretty_scen(scen)}"),
            ext = fig_ext)

terra::plot(water_plot, add = TRUE, legend = FALSE,
            col = "#74add1", ext = fig_ext)

if (nrow(rb_lines) > 0L)
  plot(sf::st_geometry(rb_lines), add = TRUE, border = "#3d3d3d", lwd = 0.4)

# demand area, dashed and heavy, so the reader sees what is being protected
plot(sf::st_geometry(sf::st_transform(aoi, terra::crs(prr_plot))),
     add = TRUE, border = "#d62728", lwd = 2.0, lty = 2)

dev.off()

message("  \u2713 Fig 2a: PRR pixel map")

# ==============================================================================
# FIG 2b: RETENTION INDEX (RI) INTERVAL PER SUB-BASIN  (two-panel)
# ==============================================================================
# Panel A: the full study area, upstream providers plus the demand corridor.
# Panel B: a demand-area zoom, which is the reason for the split. RI is
# retention x demand, so natural patches inside the demand zone score high
# simply because nothing decays over zero distance.

sb_flag <- sf::st_read(file.path(paths()$processed, "02_subbasins.gpkg"),
                       quiet = TRUE) |>
  sf::st_drop_geometry() |>
  dplyr::select(HYBAS_ID, in_downstream_aoi)
rb <- dplyr::left_join(rb, sb_flag, by = "HYBAS_ID")

rb$ri_interval <- ifelse(is.na(rb$ri_interval), "nr", rb$ri_interval)
all_lvls <- c(ri_lvls, "lake_buffer", "nr")
all_labs <- c(ri_labs, "Lake buffer", "Not ranked")
rb$ri_interval <- factor(rb$ri_interval, levels = all_lvls, labels = all_labs)

ri_pal <- c(
  grDevices::colorRampPalette(
    c("#f7fcf5", "#c7e9c0", "#74c476", "#238b45", "#00441b")
  )(length(ri_lvls)),
  "#6baed6",
  "#f0e6d2"
)
# "Lake buffer" and "Not ranked" are exclusions, not low scores, so name them
# that way in the key. This goes through the scale's `labels`, since Table 1
# joins the population CSV on the factor labels of rb$ri_interval itself.
key_labs <- all_labs
key_labs[key_labs == "Lake buffer"] <- "Lake buffer (excluded)"
key_labs[key_labs == "Not ranked"]  <- "Outside the ranked set"

ri_scale <- ggplot2::scale_fill_manual(
  values = stats::setNames(ri_pal, all_labs),
  labels = stats::setNames(key_labs, all_labs),
  name = "Realized benefit\npercentile (RI)")
ri_theme <- ggplot2::theme_void(base_size = 10) +
  ggplot2::theme(
    plot.title    = ggplot2::element_text(face = "bold", size = 11),
    plot.subtitle = ggplot2::element_text(colour = "grey40", size = 8,
                                          margin = ggplot2::margin(t = 2, b = 6)),
    legend.title  = ggplot2::element_text(face = "bold", size = 9),
    plot.margin = ggplot2::margin(4, 4, 4, 4))

# --- Panel A: full study area -------------------------------------------------
p_full <- ggplot2::ggplot() +
  ggplot2::geom_sf(data = rb,
                   ggplot2::aes(fill = ri_interval),
                   colour = "white", linewidth = 0.05) +
  ggplot2::geom_sf(data = aoi, fill = NA,
                   colour = "#d62728", linewidth = 0.7, linetype = "dashed") +
  ri_scale + ri_theme +
  ggplot2::labs(
    title = "A. Full study area",
    subtitle = paste0("Every contributing watershed. The dashed red line is ",
                      "the demand area."))

# --- Panel B: demand-area zoom -----------------------------------------------
# the watersheds scoring high on zero decay distance sit in a cluster a few
# pixels wide at panel A's scale, so the effect is invisible there
aoi_bb <- sf::st_bbox(sf::st_buffer(aoi, 5000))
p_demand <- ggplot2::ggplot() +
  ggplot2::geom_sf(data = rb,
                   ggplot2::aes(fill = ri_interval),
                   colour = "white", linewidth = 0.12) +
  ggplot2::geom_sf(data = aoi, fill = NA,
                   colour = "#d62728", linewidth = 0.8, linetype = "dashed") +
  ggplot2::coord_sf(xlim = c(aoi_bb["xmin"], aoi_bb["xmax"]),
                    ylim = c(aoi_bb["ymin"], aoi_bb["ymax"]),
                    expand = FALSE) +
  ri_scale + ri_theme +
  ggplot2::labs(
    title = "B. The demand area, enlarged",
    subtitle = paste0("The same map inside the dashed outline. Watersheds ",
                      "here rank high partly because\nnothing decays over ",
                      "zero travel distance."))

# one key for two panels. `guides = "collect"` alone leaves a legend on each
# plot; the position also has to be set on the patch, which is what `&` does.
p_ri <- patchwork::wrap_plots(p_full, p_demand, ncol = 2, guides = "collect") &
  ggplot2::theme(legend.position = "right")

# bg: theme_void() blanks the plot rectangle, so ggsave would write a
# transparent PNG, which renders black on the dark slide decks these end up in
ggplot2::ggsave(file.path(paths()$figures,
                          glue::glue("99_fig2_ri_intervals_{scen}.png")),
                p_ri, width = 14, height = 7, dpi = 200, bg = "white")

message("  \u2713 Fig 2b: RI interval map (two-panel)")

# ==============================================================================
# TABLE 1: RI INTERVAL SUMMARY × DOWNSTREAM BENEFICIARIES
# ==============================================================================

pop_tbl <- readr::read_csv(
  file.path(paths()$tables,
            glue::glue("13_benefiting_population_{scen}.csv")),
  show_col_types = FALSE)
# include the lake_buffer tier so its beneficiaries aren't dropped from Table 1.
# Same levels as rb$ri_interval below, or the join coerces both to character.
pop_tbl$ri_interval <- factor(pop_tbl$ri_interval,
                              levels = c(ri_lvls, "lake_buffer", "nr"),
                              labels = c(ri_labs, "Lake buffer", "Not ranked"))

rb_summary <- rb |>
  sf::st_drop_geometry() |>
  dplyr::group_by(ri_interval) |>
  dplyr::summarise(
    n_subbasins     = dplyr::n(),
    sum_natural_km2 = sum(natural_km2, na.rm = TRUE),
    sum_prr_total   = sum(prr_total_mm_km2, na.rm = TRUE),
    sum_tda_built   = sum(tda_built, na.rm = TRUE),
    sum_tda_crops   = sum(tda_crops, na.rm = TRUE),
    .groups = "drop") |>
  dplyr::left_join(pop_tbl, by = "ri_interval") |>
  dplyr::arrange(dplyr::desc(ri_interval))

readr::write_csv(
  rb_summary,
  file.path(paths()$tables,
            glue::glue("99_table1_realised_benefit_summary_{scen}.csv")))

message("  \u2713 Table 1: RI summary")

# ==============================================================================
# SHARED VISUAL LANGUAGE
# ==============================================================================
# one palette across the static figures and the web map, so a reader moving
# between them is not relearning the colours: slate for structure, blues for
# provision, teal and green for ecosystems, orange for the thing being argued
PAL <- c(ink = "#0f1722", slate = "#5b6b7e", mute = "#c9d4de", pale = "#e8edf2",
         bc = "#a6cbe6", blue = "#4f83ad", deep = "#1f4e73",
         teal = "#2a9d8f", green = "#2a9d57", amber = "#e9b730",
         orange = "#df744a", plum = "#7d3f9e")

theme_report <- function(base = 11) {
  ggplot2::theme_minimal(base_size = base) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold", colour = PAL[["ink"]]),
      plot.subtitle = ggplot2::element_text(colour = PAL[["slate"]], size = base - 2),
      plot.caption  = ggplot2::element_text(colour = PAL[["slate"]], size = base - 3,
                                            hjust = 0),
      axis.title    = ggplot2::element_text(colour = PAL[["slate"]], size = base - 2),
      axis.text     = ggplot2::element_text(colour = PAL[["slate"]]),
      strip.text    = ggplot2::element_text(face = "bold", colour = PAL[["ink"]]),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = PAL[["pale"]]),
      legend.position = "top")
}

# both scenarios, since the robustness figures compare them
scen_files <- list.files(paths()$processed,
                         pattern = "^12_realised_benefit_.*\\.gpkg$",
                         full.names = TRUE)
scen_ids <- sub("^12_realised_benefit_(.*)\\.gpkg$", "\\1", basename(scen_files))

rb_all <- lapply(stats::setNames(scen_files, scen_ids), function(f)
  sf::st_read(f, quiet = TRUE) |> sf::st_drop_geometry())

# ==============================================================================
# FIG 3: WHAT DRIVES RETENTION
# ==============================================================================
# retention depth against each physical driver, one panel each, with the rank
# correlation printed: tests the claim that retention concentrates in the
# wettest, steepest, most densely forested ground. Terrain and climate only, so
# these do not move when the demand weights are revised.

sb_geom <- sf::st_read(file.path(paths()$processed, "02_subbasins.gpkg"),
                       quiet = TRUE)

slope_r <- terra::rast(file.path(paths()$processed, "05_slope_deg.tif"))
precip_r <- terra::rast(file.path(paths()$processed, "precip",
                                  glue::glue("06_p_{scen}.tif")))
nat_r <- terra::rast(file.path(paths()$processed, "03_lulc_natural_mask.tif"))

drivers <- data.frame(
  HYBAS_ID = sb_geom$HYBAS_ID,
  slope_deg = exactextractr::exact_extract(slope_r, sb_geom, "mean",
                                           progress = FALSE),
  precip_mm = exactextractr::exact_extract(precip_r, sb_geom, "mean",
                                           progress = FALSE),
  nat_frac  = exactextractr::exact_extract(nat_r, sb_geom, "mean",
                                           progress = FALSE))

drv <- merge(rb_all[[scen]][, c("HYBAS_ID", "prr_per_nat_mm_km2",
                                "prr_total_mm_km2", "natural_km2")],
             drivers, by = "HYBAS_ID")
drv <- merge(drv, sf::st_drop_geometry(sb_geom)[, c("HYBAS_ID", "SUB_AREA")],
             by = "HYBAS_ID")

# retention averaged over the whole watershed, not per km² of natural land. The
# per-natural-land figure divides out the extent of the ecosystem, which is half
# of what makes a watershed good at holding runoff back, and barely moves with
# natural cover anyway (Spearman +0.05). prr_total_mm_km2 is a mm·km² total, so
# dividing by area gives mean depth in mm.
drv$retain_mm <- drv$prr_total_mm_km2 / drv$SUB_AREA
drv <- drv[is.finite(drv$retain_mm) & drv$retain_mm > 0, ]

drv_long <- drv |>
  tidyr::pivot_longer(c(precip_mm, slope_deg, nat_frac),
                      names_to = "driver", values_to = "x") |>
  dplyr::mutate(
    x = ifelse(driver == "nat_frac", x * 100, x),
    driver = factor(driver,
      levels = c("precip_mm", "slope_deg", "nat_frac"),
      labels = c("Storm precipitation (mm)", "Mean slope (degrees)",
                 # "of modelled land", not "of watershed": 03 masks the natural
                 # layer to the mapped land-cover footprint, so open water and
                 # unmapped ground are out of the denominator, as they are in 09
                 "Natural land cover (% of modelled land)")))

rho_lab <- drv_long |>
  dplyr::group_by(driver) |>
  dplyr::summarise(
    rho = stats::cor(x, retain_mm, method = "spearman",
                     use = "complete.obs"),
    x = min(x, na.rm = TRUE),
    y = max(retain_mm, na.rm = TRUE), .groups = "drop") |>
  dplyr::mutate(lab = sprintf("Spearman rho = %+.2f", rho))

p_drivers <- ggplot2::ggplot(drv_long, ggplot2::aes(x, retain_mm)) +
  ggplot2::geom_point(colour = PAL[["blue"]], alpha = 0.45, size = 1.3) +
  ggplot2::geom_smooth(method = "loess", formula = y ~ x, se = FALSE,
                       colour = PAL[["orange"]], linewidth = 0.9) +
  ggplot2::geom_text(data = rho_lab,
                     ggplot2::aes(x = x, y = y, label = lab),
                     hjust = 0, vjust = 1, size = 3.4, fontface = "bold",
                     colour = PAL[["ink"]], inherit.aes = FALSE) +
  ggplot2::facet_wrap(~driver, scales = "free_x") +
  ggplot2::labs(x = NULL, y = "Retention (mm over the watershed)") +
  theme_report()

ggplot2::ggsave(
  file.path(paths()$figures, glue::glue("99_fig3_drivers_{scen}.png")),
  p_drivers, width = 11, height = 4.6, dpi = 200, bg = "white")

message("  ✓ Fig 3: physical drivers of retention")

# ==============================================================================
# FIG 4: WHICH ECOSYSTEMS DO THE RETAINING
# ==============================================================================
# retention summed by land-cover type, as total volume and as depth per km².
# The total is what the region gets, dominated by whatever covers the most
# ground; the depth is how hard each hectare works. The conservation argument
# differs depending on which of the two puts an ecosystem near the top.

legend_tbl <- readr::read_csv(
  file.path(paths()$processed, "03_lulc_class_legend.csv"),
  show_col_types = FALSE)
lulc_r <- terra::rast(file.path(paths()$processed, "03_lulc_values.tif"))

# group the NALCMS classes into the cover types a planner would recognise
cover_group <- function(id) {
  dplyr::case_when(
    id %in% c(1, 2, 5, 6) ~ "Forest",
    id %in% c(8, 11)      ~ "Shrubland",
    id %in% c(10, 12)     ~ "Grassland",
    id == 14              ~ "Wetland",
    id == 13              ~ "Alpine / lichen-moss",
    TRUE                  ~ NA_character_)
}

eco <- lapply(stats::setNames(scen_ids, scen_ids), function(s) {
  prr_s <- terra::rast(file.path(paths()$processed, "runoff", s, "09_prr_mm.tif"))
  # mask to the watershed footprint before summing: the rasters carry a 5 km
  # working buffer, so an unmasked total overstates regional retention by about
  # a seventh and no longer reconciles with the watershed totals from 12
  prr_s <- terra::mask(prr_s, sb_mask)
  st <- terra::zonal(prr_s, lulc_r, fun = "sum", na.rm = TRUE)
  ct <- terra::zonal(terra::ifel(is.na(prr_s), NA, 1), lulc_r,
                     fun = "sum", na.rm = TRUE)
  names(st) <- c("class_id", "prr_mm_sum"); names(ct) <- c("class_id", "px")
  m <- merge(st, ct, by = "class_id")
  m$cover <- cover_group(m$class_id)
  m <- m[!is.na(m$cover), ]
  # a 30 m pixel is 900 m2, so 1 mm over one pixel is 0.9 m3
  d <- m |>
    dplyr::group_by(cover) |>
    dplyr::summarise(vol_Mm3 = sum(prr_mm_sum * 0.9, na.rm = TRUE) / 1e6,
                     km2 = sum(px, na.rm = TRUE) * 900 / 1e6,
                     .groups = "drop") |>
    dplyr::mutate(depth_mm = vol_Mm3 * 1e6 / (km2 * 1e6) * 1e3,
                  scenario = pretty_scen(s))
  d
}) |> dplyr::bind_rows()

eco$cover <- factor(eco$cover, levels = eco |>
  dplyr::filter(scenario == pretty_scen(scen)) |>
  dplyr::arrange(vol_Mm3) |> dplyr::pull(cover))

cover_pal <- c(Forest = PAL[["green"]], Wetland = PAL[["teal"]],
               Shrubland = PAL[["amber"]], Grassland = PAL[["bc"]],
               `Alpine / lichen-moss` = PAL[["slate"]])

p_eco_a <- ggplot2::ggplot(eco,
    ggplot2::aes(vol_Mm3, cover, fill = cover, alpha = scenario)) +
  ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.72),
                    width = 0.68) +
  ggplot2::scale_fill_manual(values = cover_pal, guide = "none") +
  ggplot2::scale_alpha_manual(values = c(1, 0.55), name = NULL) +
  ggplot2::labs(x = "Total runoff retained (million m³)", y = NULL) +
  theme_report()

p_eco_b <- ggplot2::ggplot(eco,
    ggplot2::aes(depth_mm, cover, fill = cover, alpha = scenario)) +
  ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.72),
                    width = 0.68) +
  ggplot2::scale_fill_manual(values = cover_pal, guide = "none") +
  ggplot2::scale_alpha_manual(values = c(1, 0.55), name = NULL) +
  ggplot2::labs(x = "Retention per unit area (mm)", y = NULL) +
  theme_report() +
  ggplot2::theme(axis.text.y = ggplot2::element_blank())

p_eco <- patchwork::wrap_plots(p_eco_a, p_eco_b, nrow = 1, widths = c(1.25, 1),
                               guides = "collect") +
  patchwork::plot_annotation(theme = theme_report())

ggplot2::ggsave(file.path(paths()$figures, "99_fig4_ecosystems.png"),
                p_eco, width = 11, height = 5, dpi = 200, bg = "white")

readr::write_csv(eco[order(eco$scenario, -eco$vol_Mm3), ],
                 file.path(paths()$tables, "99_table2_retention_by_ecosystem.csv"))

message("  ✓ Fig 4 + Table 2: retention by ecosystem")

# ==============================================================================
# FIG 5: WHICH WATERSHEDS HELP UNDER EITHER SCENARIO
# ==============================================================================
# each watershed's rank under one scenario against the other, since the ranking
# a planner acts on should not depend on which storm was modelled. The top tenth
# under both scenarios is the robust set, mapped alongside.

if (length(scen_ids) >= 2) {
  a_id <- scen_ids[1]; b_id <- scen_ids[2]
  grp_pal <- c("Top tenth under both" = PAL[["orange"]],
               "Top tenth under one only" = PAL[["amber"]],
               "Neither" = PAL[["mute"]])

  # run once per ranking measure, since the two answer different questions:
  #   prr_total_mm_km2    whole-watershed retention, so a big watershed
  #                       outranks a small one holding back more per hectare
  #   prr_per_nat_mm_km2  retention per km² of natural land, the intensity term
  #                       the realized-benefit index uses
  # Neither carries demand: both are supply, ranked before exposure or decay.
  robustness <- function(metric, file, table_file, measure, note) {
    key <- function(s) {
      d <- rb_all[[s]][, c("HYBAS_ID", metric)]
      d$rank <- rank(-d[[metric]], ties.method = "average")
      d
    }
    cmp <- merge(key(a_id), key(b_id), by = "HYBAS_ID",
                 suffixes = c(".a", ".b"))
    rho <- stats::cor(cmp$rank.a, cmp$rank.b, method = "spearman")
    n <- nrow(cmp); cut10 <- ceiling(n * 0.10)
    cmp$grp <- dplyr::case_when(
      cmp$rank.a <= cut10 & cmp$rank.b <= cut10 ~ "Top tenth under both",
      cmp$rank.a <= cut10 | cmp$rank.b <= cut10 ~ "Top tenth under one only",
      TRUE ~ "Neither")
    cmp$grp <- factor(cmp$grp, levels = names(grp_pal))

    p_scatter <- ggplot2::ggplot(cmp,
                                 ggplot2::aes(rank.a, rank.b, colour = grp)) +
      ggplot2::geom_abline(slope = 1, intercept = 0, colour = PAL[["slate"]],
                           linetype = 2, linewidth = 0.4) +
      ggplot2::geom_point(size = 1.5, alpha = 0.8) +
      ggplot2::scale_colour_manual(values = grp_pal, name = NULL) +
      ggplot2::scale_y_reverse() + ggplot2::scale_x_reverse() +
      ggplot2::annotate("text", x = n, y = 1,
                        label = sprintf("Spearman rho = %.3f", rho),
                        hjust = 0, vjust = 1, fontface = "bold",
                        colour = PAL[["ink"]], size = 3.6) +
      ggplot2::labs(
        x = glue::glue("Rank by {measure}, {pretty_scen(a_id)}"),
        y = glue::glue("Rank by {measure}, {pretty_scen(b_id)}")) +
      theme_report()

    robust_ids <- cmp$HYBAS_ID[cmp$grp == "Top tenth under both"]
    sb_map <- sb_geom
    sb_map$grp <- cmp$grp[match(sb_map$HYBAS_ID, cmp$HYBAS_ID)]
    sb_map$grp[is.na(sb_map$grp)] <- "Neither"

    p_map <- ggplot2::ggplot() +
      ggplot2::geom_sf(data = sb_map, ggplot2::aes(fill = grp),
                       colour = "white", linewidth = 0.05) +
      ggplot2::geom_sf(data = sf::st_union(aoi), fill = NA,
                       colour = PAL[["ink"]], linewidth = 0.5) +
      ggplot2::scale_fill_manual(values = grp_pal, name = NULL) +
      ggplot2::theme_void(base_size = 11) +
      ggplot2::theme(legend.position = "top",
                     plot.subtitle = ggplot2::element_text(
                       colour = PAL[["slate"]], size = 9))

    p_rob <- patchwork::wrap_plots(p_scatter, p_map, nrow = 1,
                                   widths = c(1, 1.15), guides = "collect") +
      patchwork::plot_annotation(
        title = glue::glue("Ranked by {measure}"),
        subtitle = note,
        theme = theme_report())

    ggplot2::ggsave(file.path(paths()$figures, file),
                    p_rob, width = 12, height = 6, dpi = 200, bg = "white")

    robust_out <- sf::st_drop_geometry(sb_geom)[
      match(robust_ids, sb_geom$HYBAS_ID),
      c("HYBAS_ID", "GNIS_NAME_1", "SUB_AREA", "admin_district")]
    robust_out <- merge(robust_out, cmp[, c("HYBAS_ID", "rank.a", "rank.b")],
                        by = "HYBAS_ID")
    robust_out <- robust_out[order(robust_out$rank.a), ]
    readr::write_csv(robust_out, file.path(paths()$tables, table_file))

    message("  ✓ ", file, " (rho = ", sprintf("%.3f", rho), ", ",
            length(robust_ids), " robust watersheds)")
    list(rho = rho, robust_ids = robust_ids, cmp = cmp)
  }

  rob_total <- robustness(
    "prr_total_mm_km2",
    "99_fig5_scenario_robustness.png",
    "99_table3_robust_contributors.csv",
    "total retention",
    paste0("Retention summed over the whole watershed, so watershed size is ",
           "part of the rank. Supply only: no downstream demand applied."))

  # same figure on the intensity measure. 12 caps this column at mean + 3 SD and
  # zeroes watersheds with no natural land, both of which would pile up tied
  # ranks; neither bites on the current data. Re-check if land cover changes.
  rob_per_nat <- robustness(
    "prr_per_nat_mm_km2",
    "99_fig5b_scenario_robustness_per_nat.png",
    "99_table3b_robust_contributors_per_nat.csv",
    "retention per km² of natural land",
    paste0("Retention divided by the watershed's natural area, the intensity ",
           "term the realized-benefit index uses. Supply only: no downstream ",
           "demand applied."))

  # how far the two measures disagree about who the top tenth is
  overlap <- length(intersect(rob_total$robust_ids, rob_per_nat$robust_ids))
  message("  · top-tenth sets agree on ", overlap, " of ",
          length(rob_total$robust_ids), " (total) and ",
          length(rob_per_nat$robust_ids), " (per natural km²) watersheds")
}

# ==============================================================================
message("✓ 99_figures_tables.R: figures + tables written for '", scen, "'")
