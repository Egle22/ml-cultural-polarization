# =============================================================================
# scripts/06_robustness/11_gower_pam_with_politics.R
# Purpose : Variant of 02_clustering/03_gower_pam.R where the Gower / PAM step
#           uses the value items AND the political items (258 variables)
#           instead of the 198 value items only.
# Thesis  : not reported (robustness check)
# Input   : data/processed/data_imputed_selected_pol.rds
#           (values + political block; produced by an earlier run of the
#           imputation restricted to these variables, so it is not rebuilt by
#           04_preprocess_impute.R)
# Output  : results/objects/politics_variant/{pam_results_allk,gower_dist,
#           pam2_cluster_membership,pam_mds_df_clusters}.rds
#           results/figures/06_robustness/pam_with_politics/*.pdf
# =============================================================================

source(here::here("R", "00_setup.R"))

PAM_INPUT   <- file.path(DIR_DATA, "data_imputed_selected_pol.rds")
PAM_OUT_OBJ <- file.path(DIR_OBJ, "politics_variant")
PAM_OUT_FIG <- file.path(DIR_FIG, "06_robustness", "pam_with_politics")

source(here("scripts", "02_clustering", "03_gower_pam.R"))

rm(PAM_INPUT, PAM_OUT_OBJ, PAM_OUT_FIG)
