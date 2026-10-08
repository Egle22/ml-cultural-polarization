# =============================================================================
# R/00_setup.R
# Shared setup sourced at the top of every script:
#   source(here::here("R", "00_setup.R"))
# Loads packages, defines project paths (relative to the project root via
# here::here(), which finds WVS_project.Rproj / .here), and sources the
# constants and helper functions in R/.
# =============================================================================

# -- Packages -----------------------------------------------------------------
# Loading order is kept from the original setup because later packages mask
# functions of earlier ones.
library(here)
library(dplyr)
library(tidyr)
library(readr)
library(forcats)
library(FactoMineR)
library(factoextra)
library(nFactors)        # Broken Stick, Parallel Analysis
library(cluster)         # PAM, daisy (Gower), silhouette
library(mclust)          # adjustedRandIndex (ARI)
library(randomForest)
library(ranger)          # random-forest imputation (04_preprocess_impute.R)
library(ggplot2)
library(patchwork)
library(ggalluvial)
library(plotly)
library(RColorBrewer)
library(knitr)
library(kableExtra)
library(ggnewscale)
library(scales)
library(umap)
library(deldir)
library(infotheo)
library(ggrepel)
library(clue)
library(glmnet)
library(purrr)

# -- Paths --------------------------------------------------------------------
# Raw WVS data are not redistributed: see README.md for how to obtain them.
DIR_RAW     <- here("data", "raw")
DIR_DATA    <- here("data", "processed")
DIR_OBJ     <- here("results", "objects")
DIR_CLUST   <- here("results", "objects", "clusters")  # canonical cluster labels
DIR_FIG     <- here("results", "figures")
DIR_TAB     <- here("results", "tables")

# Canonical source of the cluster labels (see CLAUDE.md).
PATH_CLUSTERS <- file.path(DIR_CLUST, "all_clusters_aligned.rds")

# Create a directory if needed and return its path.
ensure_dir <- function(...) {
  path <- file.path(...)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}

# -- Constants and helpers ----------------------------------------------------
source(here("R", "constants.R"))   # CLUSTER_COLORS, missing codes, sections
source(here("R", "labels.R"))      # fct(), fct_ord(), fct_10()
source(here("R", "codebook.R"))    # build_codebook(), print_codebook_table()
source(here("R", "missing.R"))     # get_mode(), missing-value summaries
source(here("R", "sections.R"))    # build_sections()
source(here("R", "density.R"))     # density_hdr(), plot_hdr()
source(here("R", "alignment.R"))   # hungarian_align(), compute_concordance(), ...
