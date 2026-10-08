# =============================================================================
# scripts/07_figures/07_pam_in_pca_space.R
# Purpose : Display the PAM (Gower) partitions in the z-score PCA space using
#           the ALIGNED labels: 2D scatter, Voronoi on the projected medoids,
#           3D scatter (PC1-PC3), and Voronoi on cluster means in PCA and MDS
#           space.
# Thesis  : Figure 4.7 (PAM panel)
# Input   : results/objects/{pca_model,pca_scores,pam_results_allk,
#           pam_mds_df_clusters,gower_dist}.rds
#           results/objects/clusters/all_clusters_aligned.rds
# Output  : results/objects/pam_pca_df_clusters.rds
#           results/figures/07_figures/pam/*.pdf
# =============================================================================

source(here::here("R", "00_setup.R"))
library(scatterplot3d)

out_dir_pam   <- DIR_OBJ
out_dir_pam_p <- ensure_dir(DIR_FIG, "07_figures", "pam")

pca_model            <- readRDS(file.path(DIR_OBJ, "pca_model.rds"))   # prcomp
pca_scores           <- readRDS(file.path(DIR_OBJ, "pca_scores.rds"))  # scores
results_pam          <- readRDS(file.path(DIR_OBJ, "pam_results_allk.rds"))
mds_df               <- readRDS(file.path(DIR_OBJ, "pam_mds_df_clusters.rds"))
all_clusters_aligned <- readRDS(PATH_CLUSTERS)

n_centers <- 2:6   # the values of k kept in all_clusters_aligned

scores_mat <- if (inherits(pca_scores, "prcomp")) pca_scores$x else as.matrix(pca_scores)

var_exp_pca <- round(pca_model$sdev^2 / sum(pca_model$sdev^2) * 100, 1)
cat("PCA variance explained: ", var_exp_pca[1], "% /", var_exp_pca[2], "%\n")

# MDS variance explained by the first 2 axes (same computation as 03_gower_pam.R)
gower_dist <- readRDS(file.path(DIR_CLUST, "gower_dist.rds"))
mds_fit    <- cmdscale(as.matrix(gower_dist), k = 2, eig = TRUE)
var_exp    <- round(mds_fit$eig[1:2] / sum(mds_fit$eig[mds_fit$eig > 0]) * 100, 1)

# row-by-row alignment with the PAM clusters
stopifnot(nrow(scores_mat) == nrow(all_clusters_aligned))

pca_df <- data.frame(
  PC1 = scores_mat[, 1],
  PC2 = scores_mat[, 2],
  PC3 = scores_mat[, 3]
)

# Aligned PAM labels (also used for the MDS panels)
for (k in n_centers) {
  pca_df[[paste0("cluster_k", k)]] <- factor(
    all_clusters_aligned[[paste0("pam_", k)]]
  )
  mds_df[[paste0("cluster_k", k)]] <- factor(
    all_clusters_aligned[[paste0("pam_", k)]]
  )
}

# -- 1. PAM scatter in PCA space ----------------------------------------------
plot_pam_pca <- function(df, k) {
  col_name <- paste0("cluster_k", k)
  pal <- cluster_colors

  ggplot(df, aes(x = PC1, y = PC2, color = .data[[col_name]])) +
    geom_point(
      data = df,
      aes(x = PC1, y = PC2, color = .data[[col_name]]),
      alpha = 0.45, size = 1.5
    ) +
    scale_color_manual(values = pal, name = "Cluster") +
    labs(
      title = sprintf("PAM (Gower) in PCA space — K = %d", k),
      x = sprintf("PC1 (%.1f%%)", var_exp_pca[1]),
      y = sprintf("PC2 (%.1f%%)", var_exp_pca[2])
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

plots_pca <- lapply(n_centers, function(k) plot_pam_pca(pca_df, k))
names(plots_pca) <- paste0("k", n_centers)

plot_1_pca <- (plots_pca$k2 | plots_pca$k3) / (plots_pca$k4 | plots_pca$k5) +
  plot_annotation(
    title   = "K-medoids (PAM) clustering — PCA projection",
    caption = sprintf("Dissimilarity: Gower (labels)  |  PC1 %.1f%% / PC2 %.1f%%",
                      var_exp_pca[1], var_exp_pca[2]),
    theme   = theme(plot.title = element_text(size = 13, face = "bold"))
  )

plot_2_pca <- (plots_pca$k5 | plots_pca$k6) +
  plot_annotation(
    title   = "K-medoids (PAM) clustering — PCA projection",
    caption = sprintf("Dissimilarity: Gower (labels)  |  PC1 %.1f%% / PC2 %.1f%%",
                      var_exp_pca[1], var_exp_pca[2]),
    theme   = theme(plot.title = element_text(size = 13, face = "bold"))
  )

# -- 2. Voronoi in PCA space on the PAM medoids projected onto PC1-PC2 ---------
plot_pam_pca_voronoi <- function(df, k) {
  col_name <- paste0("cluster_k", k)
  pal <- cluster_colors

  # actual PAM medoids projected on the PCA plane (not visual centroids)
  med_ids <- results_pam[[as.character(k)]]$id.med
  centroids <- data.frame(PC1 = scores_mat[med_ids, 1],
                          PC2 = scores_mat[med_ids, 2])

  pad  <- 0.1
  xlim <- range(df$PC1) + c(-pad, pad) * diff(range(df$PC1))
  ylim <- range(df$PC2) + c(-pad, pad) * diff(range(df$PC2))

  vor   <- deldir(centroids$PC1, centroids$PC2, rw = c(xlim, ylim))
  tiles <- tile.list(vor)

  vor_df <- do.call(rbind, lapply(seq_along(tiles), function(i) {
    data.frame(x = tiles[[i]]$x, y = tiles[[i]]$y, cluster = factor(i))
  }))

  ggplot() +
    geom_polygon(
      data = vor_df,
      aes(x = x, y = y, fill = cluster, group = cluster),
      alpha = 0.25, color = "grey40", linewidth = 0.35
    ) +
    geom_point(
      data = df,
      aes(x = PC1, y = PC2, color = .data[[col_name]]),
      alpha = 0.45, size = 0.7
    ) +
    geom_point(
      data = centroids,
      aes(x = PC1, y = PC2),
      shape = 3, size = 3, color = "black", stroke = 1.2
    ) +
    scale_color_manual(values = pal, name = "Cluster") +
    scale_fill_manual( values = pal, guide = "none") +
    labs(
      title = sprintf("PAM (Gower) + Voronoi in PCA space — K = %d", k),
      x = sprintf("PC1 (%.1f%%)", var_exp_pca[1]),
      y = sprintf("PC2 (%.1f%%)", var_exp_pca[2])
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

plots_pca_vor <- lapply(n_centers, function(k) plot_pam_pca_voronoi(pca_df, k))
names(plots_pca_vor) <- paste0("k", n_centers)

plot_1_pca_vor <- (plots_pca_vor$k2 | plots_pca_vor$k3) / (plots_pca_vor$k4 | plots_pca_vor$k5) +
  plot_annotation(
    title   = "K-medoids (PAM) — Voronoi tessellation in PCA space",
    caption = sprintf("Dissimilarity: Gower (labels)  |  PC1 %.1f%% / PC2 %.1f%%  |  + = medoid",
                      var_exp_pca[1], var_exp_pca[2]),
    theme   = theme(plot.title = element_text(size = 13, face = "bold"))
  )

plot_2_pca_vor <- (plots_pca_vor$k5 | plots_pca_vor$k6) +
  plot_annotation(
    title   = "K-medoids (PAM) — Voronoi tessellation in PCA space",
    caption = sprintf("Dissimilarity: Gower (labels)  |  PC1 %.1f%% / PC2 %.1f%%  |  + = medoid",
                      var_exp_pca[1], var_exp_pca[2]),
    theme   = theme(plot.title = element_text(size = 13, face = "bold"))
  )

saveRDS(pca_df, file = file.path(out_dir_pam, "pam_pca_df_clusters.rds"))

ggsave(file.path(out_dir_pam_p, "pam_pca_scatter_grid.pdf"),  plot = plot_1_pca,     width = 12, height = 9)
ggsave(file.path(out_dir_pam_p, "pam_pca_scatter_grid2.pdf"), plot = plot_2_pca,     width = 12, height = 9)
ggsave(file.path(out_dir_pam_p, "pam_pca_voronoi_grid.pdf"),  plot = plot_1_pca_vor, width = 12, height = 9)
ggsave(file.path(out_dir_pam_p, "pam_pca_voronoi_grid2.pdf"), plot = plot_2_pca_vor, width = 12, height = 9)

# -- 3. 3D scatter (PC1-PC3) ---------------------------------------------------
plot_pam_pca_3d <- function(df, k) {
  col_name <- paste0("cluster_k", k)
  cols <- CLUSTER_COLORS[as.character(df[[col_name]])]  # index by name

  scatterplot3d(
    x    = df$PC1,
    y    = df$PC2,
    z    = df$PC3,
    color       = cols,
    pch         = 16,
    cex.symbols = 0.4,
    xlab = sprintf("PC1 (%.1f%%)", var_exp_pca[1]),
    ylab = sprintf("PC2 (%.1f%%)", var_exp_pca[2]),
    zlab = sprintf("PC3 (%.1f%%)", var_exp_pca[3]),
    main = sprintf("PAM (Gower) in PCA 3D space K = %d", k),
    grid = TRUE,
    box  = FALSE
  )

  legend(
    "topright",
    legend = paste("Cluster", seq_len(k)),
    col    = CLUSTER_COLORS[as.character(seq_len(k))],
    pch    = 16,
    cex    = 0.7,
    bty    = "n"
  )
}

out_dir_single <- ensure_dir(out_dir_pam_p, "single_cluster")
# one pdf for each k
for (k in n_centers) {
  pdf(file.path(out_dir_single, sprintf("pam_pca_3d_k%d.pdf", k)),
      width = 8, height = 7)
  plot_pam_pca_3d(pca_df, k)
  dev.off()
}

# 2 x 3 grid in a single pdf
pdf(file.path(out_dir_single, "pam_pca_3d_grid.pdf"), width = 14, height = 10)
par(mfrow = c(2, 3))
for (k in n_centers) plot_pam_pca_3d(pca_df, k)
dev.off()

# -- 4. Voronoi on cluster means (aligned labels) in PCA and MDS space ---------
plot_pam_pca_voronoi_medoids <- function(df, k) {
  col_name    <- paste0("cluster_k", k)
  cluster_ids <- as.integer(as.character(df[[col_name]]))

  # after alignment the cluster labels are remapped, so each aligned
  # cluster is represented by its mean position in PCA space
  centroids <- do.call(rbind, lapply(seq_len(k), function(cl) {
    idx <- which(cluster_ids == cl)
    data.frame(
      PC1     = mean(scores_mat[idx, 1]),
      PC2     = mean(scores_mat[idx, 2]),
      cluster = factor(cl)
    )
  }))

  pad  <- 0.1
  xlim <- range(df$PC1) + c(-pad, pad) * diff(range(df$PC1))
  ylim <- range(df$PC2) + c(-pad, pad) * diff(range(df$PC2))

  vor   <- deldir(centroids$PC1, centroids$PC2, rw = c(xlim, ylim))
  tiles <- tile.list(vor)

  vor_df <- do.call(rbind, lapply(seq_along(tiles), function(i) {
    data.frame(x = tiles[[i]]$x, y = tiles[[i]]$y, cluster = factor(i))
  }))

  ggplot() +
    geom_polygon(
      data = vor_df,
      aes(x = x, y = y, fill = cluster, group = cluster),
      alpha = 0.25, color = "grey40", linewidth = 0.35
    ) +
    geom_point(
      data = df,
      aes(x = PC1, y = PC2, color = factor(cluster_ids)),
      alpha = 0.45, size = 1.5
    ) +
    geom_point(
      data = centroids,
      aes(x = PC1, y = PC2),
      shape = 3, size = 4, color = "black", stroke = 1.5
    ) +
    scale_color_manual(values = CLUSTER_COLORS[seq_len(k)], name = "Cluster") +
    scale_fill_manual( values = CLUSTER_COLORS[seq_len(k)], guide = "none") +
    labs(
      title = sprintf("PAM Voronoi (medoids) in PCA space -- K = %d", k),
      x = sprintf("PC1 (%.1f%%)", var_exp_pca[1]),
      y = sprintf("PC2 (%.1f%%)", var_exp_pca[2])
    ) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position  = if (k <= 4) "right" else "bottom",
      legend.key.size  = unit(0.4, "cm"),
      panel.grid.minor = element_blank()
    ) +
    guides(color = guide_legend(
      nrow         = 2,
      override.aes = list(size = 2.5, alpha = 1)
    ))
}

plot_pam_mds_voronoi_medoids <- function(df, k) {
  col_name    <- paste0("cluster_k", k)
  cluster_ids <- as.integer(as.character(df[[col_name]]))

  # centroids computed as means of each aligned cluster
  centroids <- do.call(rbind, lapply(seq_len(k), function(cl) {
    idx <- which(cluster_ids == cl)
    data.frame(
      MDS1    = mean(df$MDS1[idx]),
      MDS2    = mean(df$MDS2[idx]),
      cluster = factor(cl)
    )
  }))

  pad  <- 0.1
  xlim <- range(df$MDS1) + c(-pad, pad) * diff(range(df$MDS1))
  ylim <- range(df$MDS2) + c(-pad, pad) * diff(range(df$MDS2))

  vor   <- deldir(centroids$MDS1, centroids$MDS2, rw = c(xlim, ylim))
  tiles <- tile.list(vor)

  vor_df <- do.call(rbind, lapply(seq_along(tiles), function(i) {
    data.frame(x = tiles[[i]]$x, y = tiles[[i]]$y, cluster = factor(i))
  }))

  ggplot() +
    geom_polygon(
      data = vor_df,
      aes(x = x, y = y, fill = cluster, group = cluster),
      alpha = 0.25, color = "grey40", linewidth = 0.35
    ) +
    geom_point(
      data = df,
      aes(x = MDS1, y = MDS2, color = factor(cluster_ids)),
      alpha = 0.45, size = 1.5
    ) +
    geom_point(
      data = centroids,
      aes(x = MDS1, y = MDS2),
      shape = 3, size = 4, color = "black", stroke = 1.5
    ) +
    scale_color_manual(values = CLUSTER_COLORS[seq_len(k)], name = "Cluster") +
    scale_fill_manual( values = CLUSTER_COLORS[seq_len(k)], guide = "none") +
    labs(
      title = sprintf("PAM Voronoi (medoids) in MDS space -- K = %d", k),
      x = sprintf("MDS1 (%.1f%%)", var_exp[1]),
      y = sprintf("MDS2 (%.1f%%)", var_exp[2])
    ) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position  = if (k <= 4) "right" else "bottom",
      legend.key.size  = unit(0.4, "cm"),
      panel.grid.minor = element_blank()
    ) +
    guides(color = guide_legend(
      nrow         = 2,
      override.aes = list(size = 2.5, alpha = 1)
    ))
}

plots_pca_vor_med <- lapply(n_centers, function(k) plot_pam_pca_voronoi_medoids(pca_df, k))
names(plots_pca_vor_med) <- paste0("k", n_centers)

plots_mds_vor_med <- lapply(n_centers, function(k) plot_pam_mds_voronoi_medoids(mds_df, k))
names(plots_mds_vor_med) <- paste0("k", n_centers)

plot_pca_vor_med_grid <- (plots_pca_vor_med$k2 | plots_pca_vor_med$k3) /
  (plots_pca_vor_med$k4 | plots_pca_vor_med$k6) +
  plot_annotation(
    title   = "PAM Voronoi (medoids) -- PCA space",
    caption = sprintf("Dissimilarity: Gower  |  PC1 %.1f%% / PC2 %.1f%%  |  + = medoid",
                      var_exp_pca[1], var_exp_pca[2]),
    theme   = theme(plot.title = element_text(size = 13, face = "bold"))
  )

plot_mds_vor_med_grid <- (plots_mds_vor_med$k2 | plots_mds_vor_med$k3) /
  (plots_mds_vor_med$k4 | plots_mds_vor_med$k5) +
  plot_annotation(
    title   = "PAM Voronoi (medoids) -- MDS space",
    caption = sprintf("Dissimilarity: Gower  |  MDS1 %.1f%% / MDS2 %.1f%%  |  + = medoid",
                      var_exp[1], var_exp[2]),
    theme   = theme(plot.title = element_text(size = 13, face = "bold"))
  )

ggsave(file.path(out_dir_single, "pam_pca_voronoi_medoids_grid.pdf"),
       plot = plot_pca_vor_med_grid, width = 12, height = 9)
ggsave(file.path(out_dir_single, "pam_mds_voronoi_medoids_grid.pdf"),
       plot = plot_mds_vor_med_grid, width = 12, height = 9)
