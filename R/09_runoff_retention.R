# ==============================================================================
# 09_runoff_retention.R: SCS-CN runoff and potential runoff retention
# ------------------------------------------------------------------------------
# Runs the SCS Curve Number equation for each precipitation scenario produced by
# 06_precipitation.R
#
# For each scenario:
#   S = (1000 / CN) − 10                     potential retention (in)
#   Q = (P − 0.2S)^2 / (P + 0.8S)            runoff (in), if P > 0.2S
#     = 0                                    otherwise
#
#   Q_baseline       = SCS-CN(P, baseline CN)
#   Q_counterfactual = SCS-CN(P, counterfactual CN)
#   PRR (in)         = Q_counterfactual − Q_baseline
#   PRR (mm)         = PRR (in) × 25.4
#
# Precipitation is supplied in millimetres, converted to inches for the SCS
# equation, then converted back to millimetres for all output rasters.
#
# Inputs (data/processed/):
#   precip/06_p_<scenario>.tif
#   08_cn_baseline_slopeadj.tif
#   08_cn_counterfactual_slopeadj.tif
#
# Outputs (data/processed/runoff/<scenario>/):
#   09_q_baseline_mm.tif
#   09_prr_mm.tif
# ==============================================================================

source(here::here("R", "00_setup.R"))

scs_q_in <- function(P_in, CN) {
  S <- terra::ifel(CN > 0, 1000 / CN - 10, NA)
  # clamp negative S (CN just over 100 after the slope adjustment) to 0 while
  # keeping NA where CN is unmapped, so pixels with no curve number (snow/ice,
  # and cranberry, which cn_lookup treats as managed water) drop out of the
  # runoff surface instead of being modelled
  S <- terra::clamp(S, lower = 0)
  terra::ifel(P_in > 0.2 * S, (P_in - 0.2 * S)^2 / (P_in + 0.8 * S), 0)
}

# slope-adjusted AMC II curve numbers from 08: one pair, used by every scenario
cn <- list(
  b = terra::rast(file.path(paths()$processed, "08_cn_baseline_slopeadj.tif")),
  c = terra::rast(file.path(paths()$processed,
                            "08_cn_counterfactual_slopeadj.tif")))

precip_dir <- file.path(paths()$processed, "precip")
runoff_dir <- file.path(paths()$processed, "runoff")
dir.create(runoff_dir, showWarnings = FALSE, recursive = TRUE)

scenarios <- list.files(precip_dir, pattern = "^06_p_.*\\.tif$", full.names = TRUE)

blues <- grDevices::hcl.colors(50, palette = "Blues 3", rev = TRUE)
greens <- grDevices::hcl.colors(50, palette = "Greens 3", rev = TRUE)

# run one scenario: baseline runoff, counterfactual runoff, and retention
run_runoff <- function(P_mm, out, label) {
  dir.create(out, showWarnings = FALSE, recursive = TRUE)
  P_in <- P_mm / 25.4
  q_b <- scs_q_in(P_in, cn$b)
  q_c <- scs_q_in(P_in, cn$c)

  q_b_mm <- q_b * 25.4
  prr <- (q_c - q_b) * 25.4 # back to mm

  safe_writeRaster(q_b_mm, file.path(out, "09_q_baseline_mm.tif"))
  safe_writeRaster(prr, file.path(out, "09_prr_mm.tif"))

  # ---- QA preview ------------------------------------------------------------
  # how much runoff this storm produces, against how much natural cover holds back
  qa_png(paste0("09_runoff_", basename(out), ".png"), ncol = 2, function() {
    op <- graphics::par(mfrow = c(1, 2))
    on.exit(graphics::par(op), add = TRUE)
    terra::plot(q_b_mm, col = blues, axes = TRUE, main = "")
    terra::plot(prr, col = greens, axes = TRUE, main = "")
  })

  message("  ✓ runoff '", basename(out), "' written")
}

for (p_tif in scenarios) {
  scen <- sub("^06_p_", "", tools::file_path_sans_ext(basename(p_tif)))
  run_runoff(terra::rast(p_tif), file.path(runoff_dir, scen), scen)
}

message("✓ 09_runoff_retention.R: completed ", length(scenarios),
        " runoff run(s)")