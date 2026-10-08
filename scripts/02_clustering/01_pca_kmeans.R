# =============================================================================
# scripts/02_clustering/01_pca_kmeans.R
# Purpose : z-score PCA of the 198 value items, then k-means on the first
#           2, 5, 15, 75 and 198 principal components for k = 2..8 and 16.
#           KM-PC2 and KM-PC5 are the configurations used in the thesis; the
#           other dimensions are kept for the concordance analysis.
# Thesis  : Sections 3.2.1, 3.3.1; Figure 4.7; Tables 4.2, A.3
# Input   : data/processed/data_imputed_selected.rds
# Output  : results/objects/{pca_model,pca_scores,pca_kmeans_objects}.rds
#           results/objects/clusters/pca_df_clusters.rds
#           results/figures/02_clustering/pca/*.pdf
# Seeds   : set.seed(42); kmeans(nstart = 50, iter.max = 100)
# =============================================================================

source(here::here("R", "00_setup.R"))

out_dir_pca   <- ensure_dir(DIR_OBJ)
out_dir_pca_p <- ensure_dir(DIR_FIG, "02_clustering", "pca")

# -- 1. Data ------------------------------------------------------------------
imputed_data <- readRDS(file.path(DIR_DATA, "data_imputed_selected.rds"))
w_data <- imputed_data %>%
  mutate(across(where(is.factor), ~ as.numeric(as.character(.))))

# -- 2. PCA -------------------------------------------------------------------
pca_model <- prcomp(w_data, center = TRUE, scale. = TRUE)
pca_scores <- as.data.frame(pca_model$x)

var_explained <- pca_model$sdev^2 / sum(pca_model$sdev^2)
var_pc1 <- round(var_explained[1] * 100, 1)
var_pc2 <- round(var_explained[2] * 100, 1)

# -- 3. K-means on different numbers of PCs -----------------------------------
# Thesis configuration: the grid that produced the saved objects used by
# 04_concordance.R (k = 2..8, 16; 2, 5, 15, 75, 198 PCs).
k_values  <- c(2, 3, 4, 5, 6, 7, 8, 16)
n_pc_list <- c(2, 5, 15, 75, 198)
# Reduced variant (KM-PC2 / KM-PC5 only, k <= 7):
# k_values  <- c(2, 3, 4, 5, 6, 7)
# n_pc_list <- c(2, 5)

# For each number of PCs, k-means for each k
set.seed(42)
kmeans_pca <- lapply(n_pc_list, function(n_pc) {
  X_pc <- pca_scores[, 1:n_pc]
  km_list <- lapply(k_values, function(k) {
    kmeans(X_pc, centers = k, nstart = 50, iter.max = 100)
  })
  names(km_list) <- paste0("k", k_values)
  km_list
})
names(kmeans_pca) <- paste0("pc", n_pc_list)

# -- 4. PC1-PC2 data frame with cluster labels --------------------------------
# Always displayed on PC1-PC2, whatever the number of PCs used for clustering
# (as UMAP uses its 2 coordinates). One data frame per number of PCs.
pca_dfs <- lapply(n_pc_list, function(n_pc) {
  df <- data.frame(
    PC1 = pca_scores$PC1,
    PC2 = pca_scores$PC2
  )
  for (k in k_values) {
    df[[paste0("cluster_k", k)]] <- factor(
      kmeans_pca[[paste0("pc", n_pc)]][[paste0("k", k)]]$cluster
    )
  }
  df
})
names(pca_dfs) <- paste0("pc", n_pc_list)

# -- 5. Scatter plot function ------------------------------------------------
plot_pca_clusters <- function(df, k, n_pc) {
  col_name <- paste0("cluster_k", k)
  pal   <- cluster_colors
  
  ggplot(df, aes(x = PC1, y = PC2, color = .data[[col_name]])) +
    geom_point(alpha = 0.45, size = 0.7) +
    scale_color_manual(values = pal, name = "Cluster") +
    geom_hline(yintercept = 0, color = "grey80", linewidth = 0.4) +
    geom_vline(xintercept = 0, color = "grey80", linewidth = 0.4) +
    labs(
      title    = sprintf("K-means on PC1:%d  —  K = %d", n_pc, k),
      x        = paste0("PC1 (", var_pc1, "%)"),
      y        = paste0("PC2 (", var_pc2, "%)")
    ) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position  = if (k <= 4) "right" else "bottom",
      legend.key.size  = unit(0.4, "cm"),
      panel.grid.minor = element_blank()
    ) +
    guides(color = guide_legend(
      nrow         = if (k <= 8) 2 else 4,
      override.aes = list(size = 2.5, alpha = 1)
    ))
}

# -- 6. Voronoi plot function ------------------------------------------------
plot_pca_voronoi <- function(df, k, n_pc) {
  col_name <- paste0("cluster_k", k)
  pal   <- cluster_colors
  
  centroids <- df %>%
    group_by(.data[[col_name]]) %>%
    summarise(PC1 = mean(PC1), PC2 = mean(PC2), .groups = "drop")
  
  pad  <- 0.1
  xlim <- range(df$PC1) + c(-pad, pad) * diff(range(df$PC1))
  ylim <- range(df$PC2) + c(-pad, pad) * diff(range(df$PC2))
  
  vor   <- deldir(centroids$PC1, centroids$PC2, rw = c(xlim, ylim))
  tiles <- tile.list(vor)
  
  vor_df <- do.call(rbind, lapply(seq_along(tiles), function(i) {
    data.frame(
      x       = tiles[[i]]$x,
      y       = tiles[[i]]$y,
      cluster = factor(i)
    )
  }))
  
  ggplot() +
    geom_polygon(
      data      = vor_df,
      aes(x = x, y = y, fill = cluster, group = cluster),
      alpha     = 0.25,
      color     = "grey40",
      linewidth = 0.35
    ) +
    geom_point(
      data  = df,
      aes(x = PC1, y = PC2, color = .data[[col_name]]),
      alpha = 0.45, size = 0.7
    ) +
    geom_point(
      data   = centroids,
      aes(x = PC1, y = PC2),
      shape  = 3, size = 3, color = "black", stroke = 1.2
    ) +
    geom_hline(yintercept = 0, color = "grey80", linewidth = 0.4) +
    geom_vline(xintercept = 0, color = "grey80", linewidth = 0.4) +
    scale_color_manual(values = pal, name = "Cluster") +
    scale_fill_manual( values = pal, guide = "none") +
    labs(
      title    = sprintf("K-means on PC1:%d + Voronoi — K = %d", n_pc, k),
      x        = paste0("PC1 (", var_pc1, "%)"),
      y        = paste0("PC2 (", var_pc2, "%)")
    ) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position  = if (k <= 4) "right" else "bottom",
      legend.key.size  = unit(0.4, "cm"),
      panel.grid.minor = element_blank()
    ) +
    guides(color = guide_legend(
      nrow         = if (k <= 8) 2 else 4,
      override.aes = list(size = 2.5, alpha = 1)
    ))
}

# -- 7. Grids of scatter and Voronoi plots for each number of PCs -------------
for (n_pc in n_pc_list) {
  df_cur <- pca_dfs[[paste0("pc", n_pc)]]

  # scatter
  plots_sc <- lapply(k_values, function(k) plot_pca_clusters(df_cur, k, n_pc))
  names(plots_sc) <- paste0("k", k_values)

  # voronoi
  plots_vo <- lapply(k_values, function(k) plot_pca_voronoi(df_cur, k, n_pc))
  names(plots_vo) <- paste0("k", k_values)

  caption_sc <- sprintf("K-means on first %d PCs  |  visualization: PC1–PC2", n_pc)
  caption_vo <- sprintf("K-means on first %d PCs  |  PC1–PC2  |  + = k-means centroid", n_pc)

  # scatter layout
  sc_1 <- (plots_sc$k2 | plots_sc$k3) / (plots_sc$k4 | plots_sc$k5) +
    plot_annotation(
      title   = sprintf("K-means on PCA (first %d PCs)", n_pc),
      caption = caption_sc,
      theme   = theme(plot.title = element_text(size = 13, face = "bold"))
    )
  
  sc_2 <- (plots_sc$k5 | plots_sc$k6) / (plots_sc$k7 | plots_sc$k8) +
    plot_annotation(
      title   = sprintf("K-means on PCA (first %d PCs)", n_pc),
      caption = caption_sc,
      theme   = theme(plot.title = element_text(size = 13, face = "bold"))
    )
  
  # voronoi layout
  vo_1 <- (plots_vo$k2 | plots_vo$k3) / (plots_vo$k4 | plots_vo$k5) +
    plot_annotation(
      title   = sprintf("K-means on PCA (first %d PCs) — Voronoi", n_pc),
      caption = caption_vo,
      theme   = theme(plot.title = element_text(size = 13, face = "bold"))
    )
  
  vo_2 <- (plots_vo$k5 | plots_vo$k6) / (plots_vo$k7 | plots_vo$k8) +
    plot_annotation(
      title   = sprintf("K-means on PCA (first %d PCs) — Voronoi", n_pc),
      caption = caption_vo,
      theme   = theme(plot.title = element_text(size = 13, face = "bold"))
    )
  
  # save plots
  ggsave(file.path(out_dir_pca_p, sprintf("pca_pc%d_kmeans_grid.pdf",   n_pc)), sc_1, width = 12, height = 9)
  ggsave(file.path(out_dir_pca_p, sprintf("pca_pc%d_kmeans_grid2.pdf",  n_pc)), sc_2, width = 12, height = 9)
  ggsave(file.path(out_dir_pca_p, sprintf("pca_pc%d_voronoi_grid.pdf",  n_pc)), vo_1, width = 12, height = 9)
  ggsave(file.path(out_dir_pca_p, sprintf("pca_pc%d_voronoi_grid2.pdf", n_pc)), vo_2, width = 12, height = 9)

  cat("Saved plots for PC1:", n_pc, "\n")
}

# -- 8. Save objects ----------------------------------------------------------
saveRDS(pca_model,   file = file.path(out_dir_pca, "pca_model.rds"))
saveRDS(pca_scores,  file = file.path(out_dir_pca, "pca_scores.rds"))
saveRDS(kmeans_pca,  file = file.path(out_dir_pca, "pca_kmeans_objects.rds"))
saveRDS(pca_dfs,     file = file.path(ensure_dir(DIR_CLUST), "pca_df_clusters.rds"))

# -- 9. Single-cluster Voronoi plots for PC2 and PC5 --------------------------
out_dir_single_2 <- ensure_dir(out_dir_pca_p, "PC2", "single_cluster2")
out_dir_single_5 <- ensure_dir(out_dir_pca_p, "PC5", "single_cluster5")

for (n_pc in c(2, 5)) {
  out_dir_cur <- if (n_pc == 2) out_dir_single_2 else out_dir_single_5
  df_cur      <- pca_dfs[[paste0("pc", n_pc)]]

  for (kv in 2:6) {
    p_vo <- plot_pca_voronoi(df_cur, kv, n_pc)
    ggsave(
      file.path(out_dir_cur, sprintf("k%d.pdf", kv)),
      p_vo, width = 7, height = 5.5
    )

    cat(sprintf("Saved: PC%d, k=%d\n", n_pc, kv))
  }
}

# -- 10. PCA density, filled-contour rendering --------------------------------
# z-score PCA scores (the thesis always uses the z-score PCA; the min-max PCA
# is only a comparison with Desmet et al., see 06_robustness/13_pca_minmax_desmet.R).
S <- pca_scores
pts_pca <- data.frame(X = S[, 1], Y = S[, 2])
stopifnot(!anyNA(pts_pca))

dd_pca <- density_hdr(pts_pca$X, pts_pca$Y)

n_bins <- 14L
# Interpolated palette: as many colours as bins, so no level falls through to NA.
blues <- colorRampPalette(RColorBrewer::brewer.pal(9, "Blues"))(n_bins)

p_pca_filled <- ggplot() +
  geom_contour_filled(data = dd_pca$grid, aes(X, Y, z = z), bins = n_bins) +
  geom_point(data = pts_pca, aes(X, Y), colour = "red", size = 0.45,
             alpha = 0.55) +
  scale_fill_manual(values = blues, guide = "none") +
  coord_equal(expand = FALSE) +
  labs(x = "PC 1", y = "PC 2") +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank())

print(p_pca_filled)

ggsave(file.path(out_dir_pca_p, "pca_density_filled.pdf"), p_pca_filled,
       width = 7, height = 6, device = cairo_pdf)
