# ==============================================================================
# packages.R — One-shot installer
# ------------------------------------------------------------------------------
# Run this once after cloning. Installs every package the pipeline depends on,
# plus a couple of GitHub-only ones.
# ==============================================================================

cran_pkgs <- c(
  # core
  "tidyverse", "here", "glue",

  # spatial core
  "sf", "terra", "tidyterra", "stars",
  "exactextractr", "units",

  # remote data access (BC layers + StatsCan census + WSC HYDAT)
  "bcmaps", "bcdata", "tidyhydat", "cancensus",

  # hydrology / topology
  "whitebox", "igraph",

  # interpolation / climate
  "gstat", "weathercan",

  # I/O & viz
  "rmarkdown", "knitr", "kableExtra", "stargazer", "scales",
  "MetBrewer", "viridis", "ggspatial", "patchwork",

  # extras used in the user's prior FRE 527 work
  "rnaturalearth", "spatstat", "raster", "sp"
)

missing <- setdiff(cran_pkgs, rownames(installed.packages()))
if (length(missing) > 0) install.packages(missing, dependencies = TRUE)

# whitebox needs its WBT binary — first run only
if (requireNamespace("whitebox", quietly = TRUE)) {
  if (!whitebox::check_whitebox_binary()) whitebox::install_whitebox()
}

# rnaturalearthhires (high-res world boundaries; useful for inset maps)
if (!requireNamespace("rnaturalearthhires", quietly = TRUE)) {
  install.packages("rnaturalearthhires",
    repos = c("https://ropensci.r-universe.dev", "https://cloud.r-project.org"))
}

message("✓ all packages installed")
