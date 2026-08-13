# ==============================================================================
# 12_realized_benefit.R: Realized benefit per sub-basin (provision * demand)
# ------------------------------------------------------------------------------
# For each scenario (one folder under data/processed/runoff/):
#
#   1. Aggregate PRR (mm*km²) per sub-basin
#   2. natural_km2_j = Σ natural-mask pixels × pixel area
#   3. PRR_per_NAT_j = PRR_total_j / natural_km2_j (mm)
#   4. RI_j = PRR_per_NAT_j × TDA_j (provision × demand)
#   5. Eligibility filter, outlier cap, percentile intervals
#
# Lake-barrier watersheds are excluded from the percentile ranking: their
# forests reduce inflow to engineered storage rather than acting as the primary
# flood defence. Their PRR is kept and labelled "lake_buffer" so the
# dam-overflow-prevention service stays visible.
#
# Inputs (data/processed/):
#   runoff/<scenario>/09_prr_mm.tif, 09_q_baseline_mm.tif
#   03_lulc_natural_mask.tif
#   02_subbasins.gpkg
#   11_tda_subbasin.gpkg
#
# Outputs (data/processed/):
#   12_realised_benefit_<scenario>.gpkg
# ==============================================================================

source(here::here("R", "00_setup.R"))

INTERVAL_QUANTILES <- c(0, 0.25, 0.50, 0.60, 0.70, 0.80, 0.85, 0.90, 0.95, 1.00)
INTERVAL_LABELS <- c("0_25","25_50","50_60","60_70","70_80","80_85","85_90","90_95","95_100")

sb <- sf::st_read(file.path(paths()$processed, "02_subbasins.gpkg"), quiet = TRUE)
tda <- sf::st_read(file.path(paths()$processed, "11_tda_subbasin.gpkg"), quiet = TRUE) |>
  sf::st_drop_geometry()
nat_mask <- terra::rast(file.path(paths()$processed, "03_lulc_natural_mask.tif"))

px_km2 <- (terra::xres(nat_mask) * terra::yres(nat_mask)) * 1e-6
sb$natural_km2 <- px_km2 * exactextractr::exact_extract(nat_mask, sb, "sum")

is_barrier <- sb$is_lake_barrier %in% TRUE

# is_barrier and the eligibility vectors below are positional, matched to sb's
# row order, so a duplicate id would make left_join return extra rows and shift
# them all silently
stopifnot(!anyDuplicated(sb$HYBAS_ID), !anyDuplicated(tda$HYBAS_ID))

scenarios <- list.dirs(file.path(paths()$processed, "runoff"),
                       recursive = FALSE, full.names = TRUE)
scenarios <- scenarios[basename(scenarios) != "runoff"]

for (scen_dir in scenarios) {
  scen <- basename(scen_dir)
  prr <- terra::rast(file.path(scen_dir, "09_prr_mm.tif"))
  q_b <- terra::rast(file.path(scen_dir, "09_q_baseline_mm.tif"))
  prr_total <- px_km2 * exactextractr::exact_extract(prr, sb, "sum")
  q_b_total <- px_km2 * exactextractr::exact_extract(q_b, sb, "sum")

  prr_per_nat <- ifelse(sb$natural_km2 > 0, prr_total / sb$natural_km2, 0)

  # outlier cap at mean + 3*SD over non-zero, non-barrier sub-basins (barriers
  # are out of the ranking below, so out of its cap too). Too few basins to fit
  # a cap means no cap, rather than an all-NaN column.
  cap_vals <- prr_per_nat[prr_per_nat > 0 & !is_barrier]
  cap <- if (sum(is.finite(cap_vals)) >= 2) {
    mean(cap_vals, na.rm = TRUE) + 3 * sd(cap_vals, na.rm = TRUE)
  } else Inf
  prr_per_nat_capped <- pmin(prr_per_nat, cap)

  joined <- dplyr::left_join(
    dplyr::tibble(HYBAS_ID = sb$HYBAS_ID,
                  natural_km2 = sb$natural_km2,
                  q_baseline_mm_km2 = q_b_total,
                  prr_total_mm_km2 = prr_total,
                  prr_per_nat_mm_km2 = prr_per_nat_capped),
    tda |> dplyr::select(HYBAS_ID, dplyr::starts_with("tda_")),
    by = "HYBAS_ID")
  joined$ri_index <- joined$prr_per_nat_mm_km2 * joined$tda_total_w

  # ---- Eligibility + lake-buffer labelling -----------------------------------
  # a basin is rankable if it has natural area, current runoff and downstream
  # demand, unless it holds a barrier lake, where the retention buffers dam
  # inflow (see the header)
  has_provision <- joined$natural_km2 > 0 &
    joined$q_baseline_mm_km2 > 0 &
    is.finite(joined$ri_index)

  eligible <- has_provision & joined$tda_total_w > 0 & !is_barrier

  # bin by percent rank (0-1) rather than by RI value: many basins sharing an RI
  # give duplicate quantile breaks, which cut() rejects
  ri_for_cut <- ifelse(eligible, joined$ri_index, NA_real_)
  joined$ri_interval <- as.character(
    cut(dplyr::percent_rank(ri_for_cut),
        breaks = INTERVAL_QUANTILES, labels = INTERVAL_LABELS,
        include.lowest = TRUE))

  # label barrier basins that do retain, so their service (less inflow, a
  # buffer against dam overflow) stays visible in the output
  joined$ri_interval[is_barrier & has_provision] <- "lake_buffer"
  joined$rank_eligible <- eligible

  # ---- Persist ---------------------------------------------------------------
  out_all <- sb |>
    dplyr::select(HYBAS_ID) |>
    dplyr::left_join(joined, by = "HYBAS_ID")

  f_all <- file.path(paths()$processed,
                     glue::glue("12_realised_benefit_{scen}.gpkg"))
  sf::st_write(out_all, f_all, delete_dsn = TRUE, quiet = TRUE)

  # ---- Summary: lake-buffer contribution ------------------------------------
  buf <- out_all[is_barrier & has_provision, ]
  buf_vol_Mm3 <- sum(buf$prr_total_mm_km2, na.rm = TRUE) / 1000
  lake_names <- stats::na.omit(
    sb$GNIS_NAME_1[match(buf$HYBAS_ID, sb$HYBAS_ID)])

  if (nrow(buf) > 0) {
    message(sprintf("  lake-buffer: %d watershed(s) retain %.1f Mm\u00B3 (%.1f km\u00B2 natural) buffering %s",
                    nrow(buf), buf_vol_Mm3,
                    sum(buf$natural_km2, na.rm = TRUE),
                    paste(unique(lake_names), collapse = ", ")))
  }

  # ---- QA preview ------------------------------------------------------------
  # the ranked tiers, with lake-buffer basins shown apart from the percentile scale
  qa_png(paste0("12_ri_interval_", scen, ".png"), function() {
    labs <- c("0-25%", "25-50%", "50-60%", "60-70%", "70-80%", "80-85%",
              "85-90%", "90-95%", "95-100%", "lake buffer", "not ranked")
    v <- terra::vect(out_all)
    v$tier <- factor(ifelse(is.na(out_all$ri_interval), "not ranked",
                            out_all$ri_interval),
                     levels = c(INTERVAL_LABELS, "lake_buffer", "not ranked"),
                     labels = labs)
    terra::plot(v, "tier", type = "classes",
                col = c(grDevices::hcl.colors(length(INTERVAL_LABELS),
                                              "YlGnBu", rev = TRUE),
                        "#6baed6", "grey93"),
                border = "grey70", lwd = 0.1, axes = TRUE, main = "")
    plot(sf::st_geometry(read_aoi("downstream")), add = TRUE,
         border = "red", lwd = 1.4)
  })

  message("  \u2713 scenario '", scen, "' \u2192 realised benefit + intervals written")
}

message("\u2713 12_realized_benefit.R \u2014 completed ", length(scenarios), " scenario(s)")
