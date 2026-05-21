# ==============================================================================
# 09_runoff_retention.R — SCS-CN runoff and "potential runoff retention"
# ------------------------------------------------------------------------------
# Additional scenarios can be added as a raster in data/processed/precip/
#
#   S = (1000 / CN) − 10                                  (potential retention, in)
#   Q = (P − 0.2 S)^2 / (P + 0.8 S)     if P > 0.2 S      (runoff, in)
#     = 0                               otherwise
#
#   Q_baseline       =  SCS-CN(P, CN_baseline_slopeadj)
#   Q_counterfactual =  SCS-CN(P, CN_counterfactual_slopeadj)        (natural → barren)
#   PRR (in)         =  Q_counterfactual − Q_baseline
#   PRR (mm)         =  PRR (in) × 25.4
#
# P is supplied in mm by 06_precipitation.R and converted to inches before
# applying the SCS equation, then PRR is converted back to mm; downstream rasters are in metric.
#
# Inputs (data/processed/):
#   precip/06_p_<scenario>.tif
#   08_cn_baseline_slopeadj.tif
#   08_cn_counterfactual_slopeadj.tif
#
# Outputs (data/processed/runoff/<scenario>/):
#   09_q_baseline_mm.tif
#   09_q_counterfactual_mm.tif
#   09_prr_mm.tif
#   09_runoff_pct_increase.tif
# ==============================================================================

source(here::here("R", "00_setup.R"))

scs_q_in <- function(P_in, CN) {
  # USDA TR-55 (1986), with P and S in inches; matches Duarte
  # `02_spa/04_runoff.R` exactly: S = (1000/CN) - 10.
  S <- terra::ifel(CN > 0, 1000 / CN - 10, NA)
  S <- terra::ifel(S < 0, 0, S)
  terra::ifel(P_in > 0.2 * S, (P_in - 0.2 * S)^2 / (P_in + 0.8 * S), 0)
}

cn_b <- terra::rast(file.path(paths()$processed, "08_cn_baseline_slopeadj.tif"))
cn_c <- terra::rast(file.path(paths()$processed, "08_cn_counterfactual_slopeadj.tif"))

precip_dir <- file.path(paths()$processed, "precip")
runoff_dir <- file.path(paths()$processed, "runoff")
dir.create(runoff_dir, showWarnings = FALSE, recursive = TRUE)

scenarios <- list.files(precip_dir, pattern = "^06_p_.*\\.tif$", full.names = TRUE)
if (length(scenarios) == 0) {
  stop("no precipitation scenarios in ", precip_dir, " — run 06_precipitation.R")
}

for (p_tif in scenarios) {
  scen <- sub("^06_p_", "", tools::file_path_sans_ext(basename(p_tif)))
  out  <- file.path(runoff_dir, scen)
  dir.create(out, showWarnings = FALSE, recursive = TRUE)

  P_mm <- terra::rast(p_tif)
  P_in <- P_mm / 25.4

  q_b <- scs_q_in(P_in, cn_b)
  q_c <- scs_q_in(P_in, cn_c)

  prr <- (q_c - q_b) * 25.4  # back to mm
  q_b_mm <- q_b * 25.4
  q_c_mm <- q_c * 25.4

  pct_inc <- terra::ifel(q_b_mm > 0, 100 * (q_c_mm - q_b_mm) / q_b_mm, NA)
  pct_inc <- terra::ifel(pct_inc > 1000, 1000, pct_inc)  # cap for plotting

  safe_writeRaster(q_b_mm,  file.path(out, "09_q_baseline_mm.tif"))
  safe_writeRaster(q_c_mm,  file.path(out, "09_q_counterfactual_mm.tif"))
  safe_writeRaster(prr,     file.path(out, "09_prr_mm.tif"))
  safe_writeRaster(pct_inc, file.path(out, "09_runoff_pct_increase.tif"))

  # ---- QA preview ------------------------------------------------------------
  # "what runs off today vs. what natural ecosystems hold back" comparison
  blues  <- grDevices::hcl.colors(50, palette = "Blues 3",  rev = TRUE)
  greens <- grDevices::hcl.colors(50, palette = "Greens 3", rev = TRUE)
  
  qa_png(paste0("09_runoff_", scen, ".png"), ncol = 2, function() {
    op <- graphics::par(mfrow = c(1, 2), mar = c(2, 2, 3, 6))
    on.exit(graphics::par(op), add = TRUE)
    terra::plot(q_b_mm, main = paste0("Q baseline (mm) — ", scen),  col = blues)
    terra::plot(prr,    main = paste0("PRR (mm) — ",        scen),  col = greens)
  })

  message("  ✓ runoff scenario '", scen, "' written to ", basename(out))
}

message("✓ 09_runoff_retention.R — completed ", length(scenarios), " scenario(s)")