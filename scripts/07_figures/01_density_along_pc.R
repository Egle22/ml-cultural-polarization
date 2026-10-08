# =============================================================================
# scripts/07_figures/01_density_along_pc.R
# Purpose : Densities of the clusters along PC1 and PC2, KM-UMAP vs KM-PC5,
#           k = 2. Cluster colours are re-ordered by increasing mean PC1, so
#           the raw or aligned labels give the same figure.
# Thesis  : Figure 4.3 (b)
# Input   : results/objects/pca_scores.rds
#           results/objects/clusters/all_clusters_aligned.rds
# Output  : results/figures/07_figures/density_pc1_pc2_k2.pdf
# =============================================================================

#### SETUP #####################################################################

source(here::here("R", "00_setup.R"))

rm(list = intersect(c("df_pca", "p"), ls()))

pca_scores <- readRDS(file.path(DIR_OBJ, "pca_scores.rds"))
clusters   <- readRDS(PATH_CLUSTERS)
plot_dir   <- ensure_dir(DIR_FIG, "07_figures")

# Only the first two components now
df_pca <- data.frame(
  PC1 = pca_scores[, 1],
  PC2 = pca_scores[, 2]
)

stopifnot(nrow(df_pca) == nrow(clusters))


#### HELPER: relabel clusters by increasing mean PC1 ###########################

align_colors <- function(labels, k) {
  labels <- as.character(labels)
  means  <- data.frame(cluster = labels, PC1 = df_pca$PC1) |>
    dplyr::group_by(cluster) |>
    dplyr::summarise(mean_pc1 = mean(PC1), .groups = "drop") |>
    dplyr::arrange(mean_pc1)
  rank_map <- setNames(as.character(seq_len(k)), means$cluster)
  unname(rank_map[labels])
}


#### PLOT ######################################################################

plot_density_comparison <- function(k) {
  
  col_umap <- paste0("km_umap_", k)
  col_pca  <- paste0("km_pc5_k", k)
  stopifnot(all(c(col_umap, col_pca) %in% names(clusters)))
  
  df_alg1 <- df_pca |>
    dplyr::mutate(
      cluster   = factor(align_colors(clusters[[col_umap]], k),
                         levels = as.character(seq_len(k))),
      algorithm = "UMAP + K-means"
    )
  
  df_alg2 <- df_pca |>
    dplyr::mutate(
      cluster   = factor(align_colors(clusters[[col_pca]], k),
                         levels = as.character(seq_len(k))),
      algorithm = "PCA + K-means"
    )
  
  df_long <- dplyr::bind_rows(df_alg1, df_alg2) |>
    tidyr::pivot_longer(c(PC1, PC2),
                        names_to = "component", values_to = "score") |>
    dplyr::mutate(
      component = factor(component, levels = c("PC1", "PC2")),
      algorithm = factor(algorithm,
                         levels = c("UMAP + K-means", "PCA + K-means"))
    )
  
  ggplot(df_long, aes(x = score, fill = cluster, color = cluster)) +
    geom_density(alpha = 0.35, linewidth = 0.6) +
    facet_grid(algorithm ~ component, scales = "free_y") +
    scale_fill_manual(values = CLUSTER_COLORS[as.character(seq_len(k))],
                      name = "Cluster") +
    scale_color_manual(values = CLUSTER_COLORS[as.character(seq_len(k))],
                       name = "Cluster") +
    geom_vline(xintercept = 0, linetype = "dashed",
               color = "grey40", linewidth = 0.4) +
    labs(
      title    = NULL,
      subtitle = NULL,
      x = "Principal component score",
      y = "Density"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      strip.text       = element_text(face = "bold", size = 10),
      legend.position  = "bottom",
      panel.grid.minor = element_blank(),
      plot.title       = element_text(face = "bold", hjust = 0.5),
      plot.subtitle    = element_text(size = 8, hjust = 0.5, color = "grey40")
    ) +
    guides(fill = guide_legend(nrow = 1), color = guide_legend(nrow = 1))
}


#### RUN & SAVE ################################################################

p <- plot_density_comparison(2)
print(p)

ggsave(file.path(plot_dir, "density_pc1_pc2_k2.pdf"),
       p, width = 9, height = 5.5)