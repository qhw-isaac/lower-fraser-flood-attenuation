# ==============================================================================
# 06_precipitation.R: Precipitation scenarios (P, in mm)
# ------------------------------------------------------------------------------
# Precipitation depth rasters (mm), one per scenario. Everything downstream
# treats them identically: 09 builds runoff/<scenario>/ for each 06_p_*.tif,
# and the interactive map uses whichever FLOOD_SCENARIO selects.
#
# Scenario families (a baseline depth plus its climate-uplift variants):
#   wettest_month   PCIC PRISM climatological wettest month
#   ar2021          November 2021 atmospheric-river storm total
#
# Uplift variants (_plus10 / _plus20, roughly 7%/°C Clausius-Clapeyron) are
# opt-in via FLOOD_PRECIP_UPLIFTS=1; only the two base scenarios build by
# default. Drop-in grids (precip/inputs/<name>.tif, mm) build as standalone
# scenarios. FLOOD_PRECIP_SCENARIOS="wettest_month,ar2021_plus20" narrows the
# build to a named subset.
#
# Inputs:
#   data_path("prism_pcic")      PCIC PRISM 800 m monthly normals
#   data_path("rdpa_ar2021")     RDPA Nov. 13-16, 2021 storm total
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

# climatological scenario: the pixel-wise wettest PRISM normal month.
#
# PCIC's BC PRISM stops at the 49th parallel, leaving the Whatcom County
# watersheds that drain north into the Fraser Valley with no rainfall. The US
# side of the same PRISM family fills the gap: same variable, same 1991-2020
# normal period, same 30 arcsec NAD83 grid.
#
# PCIC wins wherever it has a value, so no BC cell moves and the US run only
# fills what PCIC leaves empty. They are independent runs and can disagree where
# they meet, which R/extensions/prep_prism_us.R measures over the overlap band.
# A missing US grid is not an error, it just returns the BC-only scenario.
prism_max <- function() {
  bc <- terra::app(prism, fun = max, na.rm = TRUE)

  us_f <- tryCatch(data_path("prism_us"), error = function(e) NA_character_)
  if (is.na(us_f) || !file.exists(us_f)) {
    message("  ! US PRISM normals not found. wettest_month stays BC-only ",
            "(run R/extensions/prep_prism_us.R)")
    return(bc)
  }

  # max first, then align: the wettest month is a per-cell choice among twelve,
  # so taking it on the native grid keeps interpolation to one step. The BC side
  # keeps its existing align-then-max order.
  us <- terra::rast(us_f) |>
    terra::app(fun = max, na.rm = TRUE) |>
    align_to_grid(method = "bilinear")

  filled <- terra::cover(bc, us) # PCIC first, US only where PCIC is NA
  gained <- terra::global(!is.na(filled) & is.na(bc), "sum", na.rm = TRUE)[1, 1]
  message("  · US PRISM filled ", format(gained, big.mark = ","),
          " cells south of the 49th parallel")

  # ---- QA: what the US grid contributed, and whether the two runs agree ------
  # each panel is trimmed to its own data rather than to the whole AOI, which
  # would leave two thirds of it empty north of 49°N and read as a failed render
  src <- terra::ifel(!is.na(bc), 1L, terra::ifel(!is.na(us), 2L, NA))
  # arithmetic already propagates NA, so the difference exists only on the band
  # where both runs have a value: no explicit overlap mask is needed
  overlap <- us - bc
  med <- terra::global(overlap, stats::median, na.rm = TRUE)[1, 1]
  n_ov <- terra::global(!is.na(overlap), "sum", na.rm = TRUE)[1, 1]

  qa_png("06_prism_us_fill.png", ncol = 2, panel_w = 1250, function() {
    op <- graphics::par(mfrow = c(1, 2))
    on.exit(graphics::par(op), add = TRUE)
    # levels/colours are built from the values actually present: terra matches
    # them positionally, so a hard-coded pair of labels would mislabel the map
    # if the US grid turned out to fill nothing
    src_key <- c("1" = "PCIC BC PRISM", "2" = "US PRISM (fill)")
    src_col <- c("1" = "#4f83ad", "2" = "#e9b730")
    present <- as.character(sort(terra::unique(src)[, 1]))
    terra::plot(terra::trim(src), type = "classes",
                levels = unname(src_key[present]),
                col = unname(src_col[present]), axes = TRUE,
                main = paste0("A. Source of the wettest-month depth\n",
                              "US PRISM fills ", format(gained, big.mark = ","),
                              " cells the BC grid leaves empty"))
    if (is.finite(n_ov) && n_ov > 0) {
      lim <- max(abs(terra::minmax(overlap)), na.rm = TRUE)
      terra::plot(terra::trim(overlap), range = c(-lim, lim),
                  col = grDevices::hcl.colors(50, "Blue-Red 3"), axes = TRUE,
                  main = paste0("B. US minus PCIC where both runs cover the\n",
                                "same ground (mm), median ",
                                sprintf("%+.0f", med), " over ",
                                format(n_ov, big.mark = ","), " cells"))
    } else {
      graphics::plot.new()
      graphics::text(0.5, 0.5, paste0("B. The two grids do not overlap,\n",
                                      "so there is no seam to compare"),
                     cex = 1.1)
    }
  })

  filled
}

# historical validation/event scenario: November 2021 atmospheric-river event
ar2021_storm <- function() {
  r <- terra::rast(data_path("rdpa_ar2021"))

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

# climate-uplift variant: scale a base depth by a Clausius-Clapeyron factor.
# force() pins base/factor so the closures survive the build loop
uplift <- function(base, factor) {
  force(base)
  force(factor)
  function() base() * factor
}

UPLIFTS <- c(plus10 = 1.10, plus20 = 1.20)

# families: each baseline depth plus its +10% / +20% climate-uplift variants
families <- list(
  wettest_month = list(fn = prism_max,    title = "Climatological wettest month"),
  ar2021        = list(fn = ar2021_storm, title = "Nov 2021 atmospheric river")
)

# one recipe per scenario. Uplifts are opt-in; by default only each family's
# base scenario builds and flows through 09/12
build_uplifts <- nzchar(Sys.getenv("FLOOD_PRECIP_UPLIFTS", ""))

recipes <- list()
for (fam in names(families)) {
  recipes[[fam]] <- families[[fam]]$fn
  if (build_uplifts) {
    for (u in names(UPLIFTS)) {
      recipes[[paste0(fam, "_", u)]] <- uplift(families[[fam]]$fn, UPLIFTS[[u]])
    }
  }
}

# add external drop-in precipitation grids as standalone scenarios
dropins <- character(0)
for (f in list.files(inputs_dir, pattern = "\\.tif$", full.names = TRUE)) {
  nm <- tools::file_path_sans_ext(basename(f))
  dropins <- c(dropins, nm)

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
}

# ---- build scenarios ---------------------------------------------------------

for (nm in names(recipes)) {
  p <- terra::mask(recipes[[nm]](), lulc)
  safe_writeRaster(p, file.path(precip_dir, paste0("06_p_", nm, ".tif")))
  message("scenario '", nm, "' → 06_p_", nm, ".tif")
}

# ---- QA previews -------------------------------------------------------------
# one panel per family member, one panel for each drop-in. Rasters are re-read
# from disk to keep memory light.

blues <- grDevices::hcl.colors(50, palette = "Blues 3", rev = TRUE)
read_p <- function(nm) terra::rast(file.path(precip_dir, paste0("06_p_", nm, ".tif")))

uplift_labs <- paste0("+", round((UPLIFTS - 1) * 100), "%")
names(uplift_labs) <- names(UPLIFTS)

for (fam in names(families)) {
  members <- c(fam, paste0(fam, "_", names(UPLIFTS)))
  members <- members[members %in% names(recipes)]   # respect scenario subset
  if (!length(members)) next

  qa_png(paste0("06_p_", fam, ".png"), ncol = length(members), function() {
    op <- graphics::par(mfrow = c(1, length(members)))
    on.exit(graphics::par(op), add = TRUE)
    for (m in members) terra::plot(read_p(m), col = blues, axes = TRUE, main = "")
  })
}

for (nm in dropins) {
  if (!nm %in% names(recipes)) next
  qa_png(paste0("06_p_", nm, ".png"), function() {
    terra::plot(read_p(nm), col = blues, axes = TRUE, main = "")
  })
}

message(
  "✓ 06_precipitation.R: built ",
  length(recipes),
  " precip scenario(s): ",
  paste(names(recipes), collapse = ", ")
)
