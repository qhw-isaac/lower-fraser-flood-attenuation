# ==============================================================================
# packages.R
# ------------------------------------------------------------------------------
# Installs all required packages.
# ==============================================================================

cran_pkgs <- c(
  # core
  "tidyverse", "here", "glue",

  # spatial
  "sf", "terra", "exactextractr",

  # remote / boundaries
  "bcmaps",

  # topology
  "igraph",

  # reports + figs
  "rmarkdown", "knitr", "kableExtra"
)

missing <- setdiff(cran_pkgs, rownames(installed.packages()))
if (length(missing) > 0) install.packages(missing, dependencies = TRUE)

message("✓ all packages installed")