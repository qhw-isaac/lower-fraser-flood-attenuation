# ==============================================================================
# 03_soils.R — Hydrologic Soil Group (HSG) raster
# ------------------------------------------------------------------------------
# BC Soil Survey polygons mapped to USDA-HSG
# Priority: drain lookup → texture lookup → HYSOGs250m → D-fill
#
# Inputs:
#   data_path("bc_soil_survey"), data_path("hysogs250m"),
#   here("lookup/bc_drainage_to_hsg.csv")
#   here("lookup/bc_texture_to_hsg.csv")
#
# Outputs (data/processed/):
#   03_hsg.tif          — integer HSG raster (1=A, 2=B, 3=C, 4=D), working grid
#   03_hsg_source.tif   — source data set (1=BC SSI drain, 2=BC SSI texture,
#                                          3=HYSOGs250m, 4=D-fill)
# ==============================================================================

source(here::here("R", "00_setup.R"))

.pick_drainage <- function(d1, d2, d3) {
  invalid <- function(x) is.na(x) | x %in% c("", "-")
  ifelse(!invalid(d1), d1, ifelse(!invalid(d2), d2, ifelse(!invalid(d3), d3, NA_character_)))
}

drain_lookup <- readr::read_csv(
  here::here("lookup", "bc_drainage_to_hsg.csv"),
  show_col_types = FALSE
)
texture_lookup <- readr::read_csv(
  here::here("lookup", "bc_texture_to_hsg.csv"),
  show_col_types = FALSE
)

aoi_upstream_sf <- read_aoi("upstream")
template <- working_template(aoi_upstream_sf)

bc_raw <- sf::read_sf(data_path("bc_soil_survey")) |>
  sf::st_filter(aoi_upstream_sf) |>
  dplyr::mutate(
    drain = .pick_drainage(.data$Drain_1, .data$Drain_2, .data$Drain_3)
  )

# diagnostics
present_codes  <- unique(stats::na.omit(bc_raw$drain))
mapped_codes   <- drain_lookup$drain_code
unmapped_codes <- setdiff(present_codes, mapped_codes)
no_drain_n     <- sum(is.na(bc_raw$drain))

if (length(unmapped_codes) > 0L) {
  message("  ! BC SSI drainage codes in AOI not in lookup/bc_drainage_to_hsg.csv: ",
          paste(unmapped_codes, collapse = ", "),
          "  (these polygons will fall through to texture lookup)")
}
if (no_drain_n > 0L) {
  message("  ! BC SSI polygons in AOI with no Drain_1/2/3 set: ", no_drain_n,
          " (these polygons will fall through to texture lookup)")
}

soils <- bc_raw |>
  dplyr::left_join(
    drain_lookup |> dplyr::select("drain_code", "hsg"),
    by = c("drain" = "drain_code")
  ) |>
  dplyr::left_join(
    texture_lookup |> dplyr::select("texture_code", "hsg_texture"),
    by = c("TEXTURE_1" = "texture_code")
  ) |>
  dplyr::mutate(
    hsg_src_poly = dplyr::case_when(
      !is.na(hsg)         ~ 1L,  # drain mapped
      !is.na(hsg_texture) ~ 2L,  # texture fallback
      TRUE                ~ NA_integer_
    ),
    hsg = dplyr::coalesce(hsg, hsg_texture)
  ) |>
  dplyr::transmute(hsg, hsg_src_poly)

bc_r     <- terra::rasterize(terra::vect(soils), template, field = "hsg", touches = TRUE)
bc_r_src <- terra::rasterize(terra::vect(soils), template, field = "hsg_src_poly", touches = TRUE)
bc_r_src <- terra::classify(bc_r_src, cbind(-2147483648, NA))

lulc <- terra::rast(file.path(paths()$processed, "02_lulc_values.tif"))

hy <- terra::rast(data_path("hysogs250m")) |> align_to_grid(method = "near")
hy <- terra::ifel(hy %in% c(11, 12, 13, 14), 4, hy)

hsg <- terra::ifel(!is.na(bc_r), bc_r, hy)
hsg <- terra::ifel(is.na(hsg) & !is.na(lulc) & !(lulc %in% excluded_lulc), 4L, hsg)
hsg <- terra::mask(hsg, lulc %in% excluded_lulc, maskvalue = TRUE)

# hsg_src derived AFTER hsg is finalised
hsg_src <- terra::ifel(!is.na(bc_r_src), bc_r_src,
                       terra::ifel(!is.na(hy) & is.na(bc_r_src), 3L,
                                   terra::ifel(!is.na(hsg), 4L, NA)))
hsg_src <- terra::mask(hsg_src, lulc %in% excluded_lulc, maskvalue = TRUE)

safe_writeRaster(hsg, file.path(paths()$processed, "03_hsg.tif"))
safe_writeRaster(hsg_src, file.path(paths()$processed, "03_hsg_source.tif"))

src_freq <- terra::freq(hsg_src)
src_lab <- c(`1` = "BC SSI drain mapped",
             `2` = "BC SSI texture fallback",
             `3` = "HYSOGs250m fallback",
             `4` = "D-fill (no soil data)")

total_px <- sum(src_freq$count)

for (i in seq_len(nrow(src_freq))) {
  k <- as.character(src_freq$value[i])
  message(sprintf("  · %-34s %10d px  (%5.1f%%)",
                  src_lab[k], src_freq$count[i],
                  100 * src_freq$count[i] / total_px))
}

qa_png("03_hsg.png", ncol = 2, function() {
  op <- graphics::par(mfrow = c(1, 2), mar = c(2, 2, 3, 2), oma = c(0, 0, 0, 10))
  on.exit(graphics::par(op), add = TRUE)
  
  terra::plot(hsg, main = "Hydrologic Soil Group",
              type = "classes",
              levels = c("A (sand)", "B (loam)", "C (clay loam)", "D (clay)", "No data"),
              col = c("#fee08b", "#fdae61", "#d73027", "#7f0000", "lightgrey"),
              background = "lightgrey")
  
  hsg_src_labeled <- terra::categories(hsg_src,
                                       value = data.frame(
                                         id    = c(1L, 2L, 3L, 4L),
                                         label = c("BC SSI drain", "BC SSI texture", "HYSOGs250m", "D-fill")
                                       ))
  
  terra::plot(hsg_src_labeled, main = "HSG source",
              type = "classes",
              col = c("#1b9e77", "#a6d854", "#d95f02", "#7570b3"))
})

message("✓ 03_soils.R — wrote HSG raster (BC soil survey mosaic with HYSOGs250m)")