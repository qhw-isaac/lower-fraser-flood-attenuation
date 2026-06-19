# ==============================================================================
# 99_figures_tables.R — Publication-ready figures + tables
# ------------------------------------------------------------------------------
# Reads finished artefacts from data/processed/ and output/tables/:
#   - Fig 2: PRR (mm) per pixel + RI interval per sub-basin
#   - Fig 3: Population exposure by RI interval (stacked + direct-only)
#   - Table 1: RI interval summary × downstream beneficiaries
# ==============================================================================

source(here::here("R", "00_setup.R"))

suppressPackageStartupMessages({
  library(ggplot2)
})

# ---- Scenario selection ------------------------------------------------------
scen <- Sys.getenv("FLOOD_SCENARIO", unset = "wettest_month")
rb_path <- file.path(paths()$processed, glue::glue("12_realised_benefit_{scen}_all.gpkg"))

if (!file.exists(rb_path)) {
  available <- list.files(paths()$processed,
                          pattern = "^12_realised_benefit_.*_all\\.gpkg$")
  stop("scenario '", scen, "' not found. Available: ",
       paste(available, collapse = ", "))
}

rb <- sf::st_read(rb_path, quiet = TRUE)
prr <- terra::rast(file.path(paths()$processed, "runoff", scen, "09_prr_mm.tif"))

# Drop disconnected display patches smaller than min_px (figure-only cleanup).
drop_small_patches <- function(r, min_px = 100L) {
  bin <- terra::ifel(!is.na(r), 1L, NA_integer_)
  ids <- terra::patches(bin, directions = 8)
  pf <- terra::freq(ids)
  if (is.null(pf) || nrow(pf) == 0L) return(r)
  keep <- pf$value[pf$count >= min_px]
  if (length(keep) == 0L) return(r)
  terra::mask(r, terra::ifel(ids %in% keep, 1L, NA_integer_))
}

# ---- Fig 2a: PRR per pixel ---------------------------------------------------
sb_mask <- terra::rasterize(terra::vect(rb), prr, field = 1)
src_cov <- terra::rast(file.path(paths()$processed, "03_lulc_source_coverage.tif"))
nalcms <- terra::rast(data_path("nalcms_2020")) |>
  terra::project(prr, method = "near")

prec_path <- file.path(paths()$processed, "precip", glue::glue("06_p_{scen}.tif"))
if (!file.exists(prec_path)) {
  prec_path <- file.path(paths()$processed, "precip", scen, glue::glue("06_p_{scen}.tif"))
}
prec <- terra::rast(prec_path)

in_domain <- !is.na(sb_mask) & !is.na(src_cov)
model_mask <- in_domain & !is.na(prec)
water <- in_domain & nalcms == 18L
snow_ice <- in_domain & nalcms == 19L

prr_display <- terra::mask(prr, model_mask, maskvalues = c(NA, FALSE))
prr_display <- terra::ifel(is.na(prr_display) & snow_ice, 0, prr_display)
water_overlay <- terra::ifel(water, 1, NA)

# Crop to the downstream AOI's southern edge, then clean specks on the smaller extent.
down_bb <- sf::st_bbox(read_aoi("downstream"))
prr_bb <- terra::ext(prr_display)
data_ext <- terra::ext(prr_bb[1], prr_bb[2], down_bb["ymin"], prr_bb[4])

prr_plot <- terra::crop(prr_display, data_ext)
water_plot <- terra::crop(water_overlay, data_ext)

display_mask <- terra::ifel(!is.na(prr_plot) | !is.na(water_plot), 1L, NA_integer_)
display_mask <- drop_small_patches(display_mask, min_px = 100L)
display_mask <- terra::sieve(display_mask, threshold = 25L, directions = 8)
prr_plot <- terra::mask(prr_plot, display_mask)
water_plot <- terra::mask(water_plot, display_mask)

# Sub-basin outlines only on the final coloured footprint.
has_color <- display_mask
color_poly <- terra::as.polygons(has_color, dissolve = TRUE, values = FALSE)
rb_lines <- sf::st_intersection(sf::st_crop(rb, data_ext), sf::st_as_sf(color_poly))
if (nrow(rb_lines) > 0L) {
  rb_lines <- rb_lines[as.numeric(sf::st_length(rb_lines)) > 120, ]
}

# Pad extent so data doesn't sit flush against the axis frame.
pad_frac <- 0.04
xw <- data_ext[2] - data_ext[1]
yw <- data_ext[4] - data_ext[3]
fig_ext <- terra::ext(
  data_ext[1] - pad_frac * xw,
  data_ext[2] + pad_frac * xw,
  data_ext[3] - pad_frac * yw,
  data_ext[4] + pad_frac * yw
)

png(file.path(paths()$figures, glue::glue("99_fig2_prr_{scen}.png")),
    width = 1600, height = 1250, res = 200)
terra::plot(prr_plot,
            main = glue::glue(
              "Potential Runoff Retention (mm) \u2014 {scen}\n",
              "water bodies (light blue) not modelled; snow/ice = 0"),
            background = "white",
            mar = c(3.1, 3.1, 5.1, 7.1),
            ext = fig_ext)
terra::plot(water_plot, add = TRUE, legend = FALSE, col = "#74add1", ext = fig_ext)
if (nrow(rb_lines) > 0L) {
  plot(sf::st_geometry(rb_lines), add = TRUE, border = "#3d3d3d", lwd = 0.55)
}
dev.off()

# ---- Fig 2b: RI interval map -------------------------------------------------
ri_lvls <- c("0_25", "25_50", "50_60", "60_70",
             "70_80", "80_85", "85_90", "90_95", "95_100")
ri_labs <- c("0\u201325", "25\u201350", "50\u201360", "60\u201370",
             "70\u201380", "80\u201385", "85\u201390", "90\u201395", "95\u2013100")

nr_label <- "Not ranked"
rb$ri_interval <- ifelse(is.na(rb$ri_interval), nr_label, rb$ri_interval)
all_lvls <- c(ri_lvls, nr_label)
all_labs <- c(ri_labs, nr_label)
rb$ri_interval <- factor(rb$ri_interval, levels = all_lvls, labels = all_labs)

ri_pal <- c(
  colorRampPalette(c("#e8f5e0", "#b2df8a", "#33a02c", "#006d2c"))(length(ri_lvls)),
  "#f0e6d2"
)

p_ri <- ggplot(rb) +
  geom_sf(aes(fill = ri_interval), colour = "white", linewidth = 0.05) +
  scale_fill_manual(values = stats::setNames(ri_pal, all_labs),
                    name = "RI (%)") +
  theme_minimal() +
  labs(title = glue::glue("Retention Index interval per sub-basin \u2014 {scen}"),
       subtitle = "Higher RI = greater natural runoff retention")

ggsave(file.path(paths()$figures, glue::glue("99_fig2_ri_intervals_{scen}.png")),
       p_ri, width = 9, height = 7, dpi = 200)

# ---- Table 1: interval × downstream beneficiaries ----------------------------
pop_csv_path <- file.path(paths()$tables,
                          glue::glue("13_benefiting_population_{scen}.csv"))

if (file.exists(pop_csv_path)) {
  pop_tbl <- readr::read_csv(pop_csv_path, show_col_types = FALSE)
  pop_tbl$ri_interval <- factor(pop_tbl$ri_interval,
                                levels = ri_lvls, labels = ri_labs)

  rb_summary <- rb |>
    sf::st_drop_geometry() |>
    dplyr::group_by(ri_interval) |>
    dplyr::summarise(
      n_subbasins = dplyr::n(),
      sum_natural_km2 = sum(natural_km2, na.rm = TRUE),
      sum_prr_total = sum(prr_total_mm_km2, na.rm = TRUE),
      sum_tda_built = sum(tda_built, na.rm = TRUE),
      sum_tda_crops = sum(tda_crops, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::left_join(pop_tbl, by = "ri_interval") |>
    dplyr::arrange(dplyr::desc(ri_interval))

  readr::write_csv(
    rb_summary,
    file.path(paths()$tables,
              glue::glue("99_table1_realised_benefit_summary_{scen}.csv"))
  )
}

# ---- Fig 3a: Population exposure (stacked) ------------------------------------
if (file.exists(pop_csv_path)) {
  pop_long <- pop_tbl |>
    dplyr::filter(direct_pop > 0 | indirect_pop > 0) |>
    tidyr::pivot_longer(c(direct_pop, indirect_pop),
                        names_to = "type", values_to = "people") |>
    dplyr::mutate(
      type = dplyr::recode(type,
        direct_pop = "Direct (in floodplain)",
        indirect_pop = "Indirect (in pop centre)"
      ),
      type = factor(type, levels = c("Indirect (in pop centre)",
                                     "Direct (in floodplain)"))
    )

  p_pop <- ggplot(pop_long, aes(ri_interval, people, fill = type)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = c("Indirect (in pop centre)" = "#d9d9d9",
                                 "Direct (in floodplain)" = "#e6550d")) +
    scale_y_continuous(labels = scales::label_comma()) +
    labs(title = glue::glue(
           "Downstream benefiting population by RI interval \u2014 {scen}"),
         x = "RI (%)", y = "People", fill = NULL) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "top",
      panel.grid.major.x = element_blank()
    )

  ggsave(file.path(paths()$figures,
                   glue::glue("99_fig3_population_{scen}.png")),
         p_pop, width = 9, height = 5.5, dpi = 200)

  # ---- Fig 3b: Direct exposure only -------------------------------------------
  pop_direct <- pop_tbl |>
    dplyr::filter(direct_pop > 0 | indirect_pop > 0)

  p_direct <- ggplot(pop_direct, aes(ri_interval, direct_pop)) +
    geom_col(fill = "#e6550d", width = 0.6) +
    geom_text(aes(label = scales::label_comma(accuracy = 1)(round(direct_pop))),
              vjust = -0.4, size = 3.2) +
    scale_y_continuous(labels = scales::label_comma(),
                       expand = expansion(mult = c(0, 0.15))) +
    labs(title = glue::glue(
           "Direct flood-exposed population by RI interval \u2014 {scen}"),
         subtitle = "People in the 100-yr floodplain downstream of each priority tier",
         x = "RI (%)", y = "People") +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major.x = element_blank()
    )

  ggsave(file.path(paths()$figures,
                   glue::glue("99_fig3b_direct_pop_{scen}.png")),
         p_direct, width = 9, height = 5, dpi = 200)
}

message("\u2713 99_figures_tables.R \u2014 figures + tables written for scenario '",
        scen, "'")