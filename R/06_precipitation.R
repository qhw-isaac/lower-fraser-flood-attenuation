# ==============================================================================
# 06_precipitation.R — Precipitation scenarios (P, in mm)
# ------------------------------------------------------------------------------
# Builds precipitation-depth rasters for one or more flood scenarios.
#
# Each output is a raster of storm/event precipitation depth (in mm) aligned to 
# the project grid. The rest of the pipeline treats each scenario the same way:
#
#   09_runoff_retention.R   builds runoff/<scenario>/ for each 06_p_*.tif
#   12_realized_benefit.R   ranks realized benefits for each scenario
#   14_interactive_map.R    serves the scenario selected by FLOOD_SCENARIO
#
# Scenario families:
#   1. wettest_month — PCIC PRISM climatological wettest month
#   2. ar2021        — November 2021 atmospheric-river storm total
#   3. drop-in grids — any precip/inputs/<name>.tif file, in mm
#
# To build only selected scenarios, set:
#   FLOOD_PRECIP_SCENARIOS="wettest_month,ar2021"
#
# Inputs:
#   data_path("prism_pcic")      PCIC PRISM 800 m monthly normals
#   data_path("rdpa_ar2021")     RDPA Nov. 13–16, 2021 storm total
#   03_lulc_values.tif           model footprint
#   precip/inputs/*.tif          optional external precipitation grids
#
# Outputs:
#   data/processed/precip/06_p_<scenario>.tif
# ==============================================================================

source(here::here("R", "00_setup.R"))

precip_dir <- file.path(paths()$processed, "precip")
inputs_dir <- file.path(precip_dir, "inputs")
dir.create(inputs_dir, showWarnings = FALSE, recursive = TRUE)

prism <- terra::rast(data_path("prism_pcic")) |> 
  align_to_grid(method = "bilinear")

lulc <- terra::rast(file.path(paths()$processed, "03_lulc_values.tif"))

# ---- Scenario registry -------------------------------------------------------

# climatological scenario: pixel-wise wettest PRISM normal month
prism_max <- function() {
  terra::app(prism, fun = max, na.rm = TRUE)
}

# historical validation/event scenario:November 2021 atmospheric-river event
ar2021_storm <- function() {
  r <- terra::rast(data_path("rdpa_ar2021", must_exist = TRUE))
  
  # fill one-cell edge gaps before alignment
  r <- terra::focal(
    terra::extend(r, 1),
    w = 3,
    fun = mean,
    na.policy = "only",
    na.rm = TRUE
  )
  
  align_to_grid(r, method = "bilinear")
}

# convenience helper for ad-hoc climate uplift tests.
uplift <- function(base, factor) {
  function() base() * factor
}

recipes <- list(
  wettest_month = prism_max,
  ar2021        = ar2021_storm
)

# add external drop-in precipitation grids
for (f in list.files(inputs_dir, pattern = "\\.tif$", full.names = TRUE)) {
  nm <- tools::file_path_sans_ext(basename(f))
  
  local({
    ff <- f
    recipes[[nm]] <<- function() {
      align_to_grid(terra::rast(ff), method = "bilinear")
    }
  })
}

# ---- scenario subset ---------------------------------------------------------

want <- Sys.getenv("FLOOD_PRECIP_SCENARIOS", "")

if (nzchar(want)) {
  keep <- trimws(strsplit(want, ",", fixed = TRUE)[[1]])
  recipes <- recipes[intersect(keep, names(recipes))]
  
  if (!length(recipes)) {
    stop("FLOOD_PRECIP_SCENARIOS matched no known scenario(s): ", want)
  }
}

# ---- build scenarios ---------------------------------------------------------

for (nm in names(recipes)) {
  p <- terra::mask(recipes[[nm]](), lulc)
  
  out <- file.path(precip_dir, paste0("06_p_", nm, ".tif"))
  safe_writeRaster(p, out)
  
  qa_png(paste0("06_p_", nm, ".png"), function() {
    op <- graphics::par(mar = c(2, 2, 3, 6))
    on.exit(graphics::par(op), add = TRUE)
    
    terra::plot(
      p,
      main = paste0("Precipitation scenario: ", nm, " (mm)"),
      col = grDevices::hcl.colors(50, palette = "Blues 3", rev = TRUE)
    )
  })
  
  message("scenario '", nm, "' → ", basename(out))
}

message(
  "✓ 06_precipitation.R — built ",
  length(recipes),
  " precip scenario(s): ",
  paste(names(recipes), collapse = ", ")
)