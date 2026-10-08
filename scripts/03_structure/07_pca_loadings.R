# =============================================================================
# scripts/03_structure/07_pca_loadings.R
# Purpose : Ten largest absolute loadings on PC1 and PC2 of the z-score PCA
#           used throughout the thesis (same computation as the min-max
#           version in 06_robustness/13_pca_minmax_desmet.R, applied to the
#           thesis PCA).
# Thesis  : Appendix A.3.1; Table A.4
# Input   : results/objects/pca_model.rds (from 02_clustering/01_pca_kmeans.R)
#           results/objects/wvs_question_labels.rds
# Output  : console; results/tables/03_structure/pca_top_loadings.csv
# =============================================================================

source(here::here("R", "00_setup.R"))

pca_model  <- readRDS(file.path(DIR_OBJ, "pca_model.rds"))
wvs_labels <- readRDS(file.path(DIR_OBJ, "wvs_question_labels.rds"))

var_exp <- pca_model$sdev^2 / sum(pca_model$sdev^2)
cat(sprintf("Var. explained PC1 = %.1f%%, PC2 = %.1f%%, PC1+PC2 = %.1f%%\n",
            var_exp[1]*100, var_exp[2]*100, sum(var_exp[1:2])*100))

loadings <- as.data.frame(pca_model$rotation[, 1:2])
loadings$variable <- rownames(loadings)
loadings$label    <- unname(wvs_labels[loadings$variable])

top_pc1 <- loadings[order(-abs(loadings$PC1)), ][1:10, c("variable", "label", "PC1")]
top_pc2 <- loadings[order(-abs(loadings$PC2)), ][1:10, c("variable", "label", "PC2")]

cat("\n=== TOP 10 loadings PC1 (z-score) ===\n")
print(top_pc1, row.names = FALSE)

cat("\n=== TOP 10 loadings PC2 (z-score) ===\n")
print(top_pc2, row.names = FALSE)

top_loadings <- rbind(
  data.frame(component = "PC1", top_pc1[, c("variable", "label")], loading = top_pc1$PC1),
  data.frame(component = "PC2", top_pc2[, c("variable", "label")], loading = top_pc2$PC2)
)
write.csv(top_loadings,
          file.path(ensure_dir(DIR_TAB, "03_structure"), "pca_top_loadings.csv"),
          row.names = FALSE)
