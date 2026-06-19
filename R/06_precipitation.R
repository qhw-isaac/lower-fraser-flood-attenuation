# ==============================================================================
# 06_precipitation.R — Precipitation scenarios (P, in mm)
# ------------------------------------------------------------------------------
# Builds precipitation rasters P (mm) on a scenario basis. 
# Each scenario flows through the rest of the pipeline unchanged:
# 09_runoff_retention.R builds runoff/<scenario>/ for every 06_p_*.tif,
# 12_realized_benefit.R ranks each, and 14_interactive_map.R ships whichever
# FLOOD_SCENARIO selects.
#
# SCENARIO TYPES
#   1. Derived recipes (`recipes` below)
#   2. Drop-in grids (`precip/inputs/<name>.tif`, depth in mm)
#
# Restrict which scenarios build with FLOOD_PRECIP_SCENARIOS="a,b" (default all).
#
# Inputs:
#   data_path("prism_pcic") PCIC PRISM 800 m monthly normals (12 bands)
#   data/processed/03_lulc_values.tif
#   data/processed/precip/inputs/*.tif
#
# Outputs (data/processed/precip/):
#   06_p_<scenario>.tif P (mm) per scenario
# ==============================================================================
source(here::here("R", "00_setup.R"))

precip_dir <- file.path(paths()$processed, "precip")
inputs_dir <- file.path(precip_dir, "inputs")
dir.create(inputs_dir, showWarnings = FALSE, recursive = TRUE)

prism <- terra::rast(data_path("prism_pcic")) |> 
  align_to_grid(method = "bilinear")

lulc <- terra::rast(file.path(paths()$processed, "03_lulc_values.tif"))

# ---- Scenario registry -------------------------------------------------------
# name -> function() returning a precipitation raster P (mm) on the working grid
prism_max <- function() terra::app(prism, fun = max, na.rm = TRUE)

recipes <- list(
  wettest_month = prism_max, 
  wettest_month_plus10 = function() prism_max() * 1.10,
  wettest_month_plus20 = function() prism_max() * 1.20
)

# register any external drop-in grids (depth in mm) as scenarios
for (f in list.files(inputs_dir, pattern = "\\.tif$", full.names = TRUE)) {
  nm <- tools::file_path_sans_ext(basename(f))
  local({ ff <- f
    recipes[[nm]] <<- function() align_to_grid(terra::rast(ff), method = "bilinear") })
}

# Optional subset via env (comma-separated)
want <- Sys.getenv("FLOOD_PRECIP_SCENARIOS", "")
if (nzchar(want)) {
  keep <- trimws(strsplit(want, ",", fixed = TRUE)[[1]])
  recipes <- recipes[intersect(keep, names(recipes))]
  if (!length(recipes))
    stop("FLOOD_PRECIP_SCENARIOS matched no known scenario(s): ", want)
}

# ---- Build each scenario -----------------------------------------------------
for (nm in names(recipes)) {
  p <- terra::mask(recipes[[nm]](), lulc) # same footprint as the model grid
  out <- file.path(precip_dir, paste0("06_p_", nm, ".tif"))
  safe_writeRaster(p, out)

  qa_png(paste0("06_p_", nm, ".png"), function() {
    op <- graphics::par(mar = c(2, 2, 3, 6)); on.exit(graphics::par(op), add = TRUE)
    terra::plot(p, main = paste0("Precipitation scenario: ", nm, " (mm)"),
                col = grDevices::hcl.colors(50, palette = "Blues 3", rev = TRUE))
  })
  message("scenario '", nm, "' → ", basename(out))
}

message("✓ 06_precipitation.R — built ", length(recipes),
        " precip scenario(s): ", paste(names(recipes), collapse = ", "))