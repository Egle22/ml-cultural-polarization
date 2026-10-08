# =============================================================================
# scripts/06_robustness/10_umap_sensitivity.R
# Purpose : Sensitivity of KM-UMAP to the UMAP hyper-parameters: the same
#           pipeline as 02_clustering/02_umap_kmeans.R with n_neighbors = 50
#           and min_dist = 1e-6 (tighter clumps), 2D only, k = 2..6.
# Thesis  : not reported (robustness check)
# Input   : data/processed/data_imputed_selected.rds
# Output  : results/objects/UMAP_new_min_dist/{umap_kmeans_objects,
#           umap_df_clusters}.rds
#           results/figures/06_robustness/umap_sensitivity/*.pdf
# =============================================================================

source(here::here("R", "00_setup.R"))

UMAP_N_NEIGHBORS <- 50
UMAP_MIN_DIST    <- 0.000001
UMAP_EXTRA_DIMS  <- FALSE
UMAP_K_VALUES    <- c(2, 3, 4, 5, 6)
UMAP_OUT_OBJ     <- file.path(DIR_OBJ, "UMAP_new_min_dist")
UMAP_OUT_CLUST   <- UMAP_OUT_OBJ
UMAP_OUT_FIG     <- file.path(DIR_FIG, "06_robustness", "umap_sensitivity")

source(here("scripts", "02_clustering", "02_umap_kmeans.R"))

rm(UMAP_N_NEIGHBORS, UMAP_MIN_DIST, UMAP_EXTRA_DIMS, UMAP_K_VALUES,
   UMAP_OUT_OBJ, UMAP_OUT_CLUST, UMAP_OUT_FIG)
