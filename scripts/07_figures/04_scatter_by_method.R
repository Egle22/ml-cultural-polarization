# =============================================================================
# scripts/07_figures/04_scatter_by_method.R
# Purpose : Scatter and Voronoi plots of every configuration (aligned labels)
#           in UMAP and PCA space, for k = 2..6, one panel per method.
# Thesis  : Figures 4.6, 4.7 (supporting panels)
# Input   : results/objects/clusters/{all_clusters_aligned,umap_df_clusters}.rds
#           results/objects/pca_scores.rds
# Output  : results/figures/07_figures/{umap_in_pca,all_voronoi}/*.pdf
# =============================================================================

# UMAP PLOT FUNCTIONS----
source(here::here("R", "00_setup.R"))
all_cluster_aligned <- readRDS(PATH_CLUSTERS)
umap_df <- readRDS(file.path(DIR_CLUST, "umap_df_clusters.rds"))
k_values <- 2:6

umap_cols <- grep("^km_umap_", names(all_cluster_aligned), value = TRUE)

umap_df <-cbind(umap_df, all_cluster_aligned[, umap_cols , drop = F])
umap_df[umap_cols]<-lapply(umap_df[umap_cols], factor)

library(RColorBrewer)

plot_umap_clusters <- function(df, k) {
  col_name <- paste0("cluster_k", k)
  
  # palette: Set1 up to 8 clusters, interpolated beyond
  n_col <- max(min(k, 9), 3)
  pal   <- colorRampPalette(brewer.pal(n_col, "Set1"))(k)
  
  ggplot(df, aes(x = UMAP1, y = UMAP2, color = .data[[col_name]])) +
    geom_point(alpha = 0.45, size = 0.7) +
    scale_color_manual(values = cluster_colors, name = "Cluster") +
    labs(
      title    = sprintf("K-means on UMAP space  —  K = %d", k),
      subtitle = sprintf("%d clusters, nstart = 25", k),
      x = "UMAP 1", y = "UMAP 2"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position    = if (k <= 4) "right" else "bottom",
      legend.key.size    = unit(0.4, "cm"),
      panel.grid.minor   = element_blank()
    ) +
    guides(color = guide_legend(
      nrow          = if (k <= 8) 2 else 4,
      override.aes  = list(size = 2.5, alpha = 1)
    ))
}

plots <- lapply(k_values, function(k) plot_umap_clusters(umap_df, k))
names(plots) <- paste0("k", k_values)
# ── 7. Layout ────────────────────────────────────────────────────────

plot_1 <- (plots$k2 | plots$k3) / (plots$k4 | plots$k5) +
  plot_annotation(
    title   = "K-means clustering in UMAP 2D space",
    caption = "UMAP: n_neighbors = 15, min_dist = 0.1, metric = Euclidean",
    theme   = theme(plot.title = element_text(size = 13, face = "bold"))
  )



plot_2 <- (plots$k6 | plots$k7) / (plots$k7 | plots$k8) +
  plot_annotation(
    title   = "K-means clustering in UMAP 2D space",
    caption = "UMAP: n_neighbors = 15, min_dist = 0.1, metric = Euclidean",
    theme   = theme(plot.title = element_text(size = 13, face = "bold"))
  )



# ── 8. Voronoi in UMAP space ─────────────────────────────────────────####
# overwrite with the correct (aligned) clusters
# Replace the raw k-means labels with the aligned ones
k_values_aligned <- 2:6   # only those in all_clusters_aligned

for (k in k_values_aligned) {
  umap_df[[paste0("cluster_k", k)]] <- factor(
    all_clusters_aligned[[paste0("km_umap_", k)]]
  )
}
library(deldir)

plot_umap_voronoi <- function(df, k) {
  col_name <- paste0("cluster_k", k)
  
  n_col <- max(min(k, 9), 3)
  pal   <- cluster_colors
  
  # visual centroids: mean UMAP1/UMAP2 per cluster
  centroids <- df %>%
    group_by(.data[[col_name]]) %>%
    summarise(UMAP1 = mean(UMAP1), UMAP2 = mean(UMAP2), .groups = "drop")
  
  # Voronoi tessellation of the visual centroids
  pad  <- 0.1
  xlim <- range(df$UMAP1) + c(-pad, pad) * diff(range(df$UMAP1))
  ylim <- range(df$UMAP2) + c(-pad, pad) * diff(range(df$UMAP2))
  
  vor   <- deldir(centroids$UMAP1, centroids$UMAP2,
                  rw = c(xlim, ylim))
  tiles <- tile.list(vor)
  
  # tiles to a data frame for ggplot
  vor_df <- do.call(rbind, lapply(seq_along(tiles), function(i) {
    data.frame(
      x       = tiles[[i]]$x,
      y       = tiles[[i]]$y,
      cluster = factor(i)
    )
  }))
  
  ggplot() +
    # filled Voronoi cells
    geom_polygon(
      data      = vor_df,
      aes(x = x, y = y, fill = cluster, group = cluster),
      alpha     = 0.25,
      color     = "grey40",
      linewidth = 0.35
    ) +
    # observations
    geom_point(
      data  = df,
      aes(x = UMAP1, y = UMAP2, color = .data[[col_name]]),
      alpha = 0.45, size = 0.7
    ) +
    # visual centroids
    geom_point(
      data   = centroids,
      aes(x = UMAP1, y = UMAP2),
      shape  = 3, size = 3, color = "black", stroke = 1.2
    ) +
    scale_color_manual(values = pal, name = "Cluster") +
    scale_fill_manual( values = pal, guide = "none") +
    labs(
      title    = sprintf("K-means (UMAP) + Voronoi K = %d", k),
      x = "UMAP1", y = "UMAP2"
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

plots_vor <- lapply(k_values, function(k) plot_umap_voronoi(umap_df, k))
names(plots_vor) <- paste0("k", k_values)

# ── 9. Layout Voronoi ────────────────────────────────────────────────
plot_1_vor <- (plots_vor$k2 | plots_vor$k3) / (plots_vor$k4 | plots_vor$k5) +
  plot_annotation(
    title   = "K-means clustering in UMAP 2D space — Voronoi tessellation",
    caption = "UMAP: n_neighbors = 15, min_dist = 0.1, metric = Euclidean ",
    theme   = theme(plot.title = element_text(size = 13, face = "bold"))
  )

plot_2_vor <- (plots_vor$k5 | plots_vor$k6) / (plots_vor$k7 | plots_vor$k8) +
  plot_annotation(
    title   = "K-means clustering in UMAP 2D space — Voronoi tessellation",
    caption = "UMAP: n_neighbors = 15, min_dist = 0.1, metric = Euclidean ",
    theme   = theme(plot.title = element_text(size = 13, face = "bold"))
  )

# pca ####
# ── 4. PC1-PC2 data frame with cluster labels ──────────────────────────
# always displayed on PC1-PC2, whatever the number of PCs used
# for clustering (as UMAP uses its 2 coordinates)

pca_scores <- readRDS(file.path(DIR_OBJ, "pca_scores.rds"))
n_pc_list <- c(2,5)
# one data frame for each number of PCs
pca_dfs <- lapply(n_pc_list, function(n_pc) {
  df <- data.frame(
    PC1 = pca_scores$PC1,
    PC2 = pca_scores$PC2
  )
  for (k in k_values) {
    key <- paste0("km_pc", n_pc, "_k", k)          # e.g. "km_pc2_k3"
    df[[paste0("cluster_k", k)]] <- factor(
      all_cluster_aligned[[key]]            # aligned label vector

    )
  }
  df
})
names(pca_dfs) <- paste0("pc", n_pc_list)

# ── 5. Scatter plot function ─────────────────────────────────────────
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
      x        = "PC1",
      y        = "PC2"
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

# ── 6. Voronoi plot function ─────────────────────────────────────────
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
      x        = "PC1",
      y        = "PC2"
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

# ── 7. All plots and layouts for each number of PCs ──────────
library(deldir)

# for each n_pc: scatter + voronoi, two grids each
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
  
  sc_2 <- (plots_sc$k5 | plots_sc$k6)  +
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
  
  vo_2 <- (plots_vo$k5 | plots_vo$k6) +
    plot_annotation(
      title   = sprintf("K-means on PCA (first %d PCs) — Voronoi", n_pc),
      caption = caption_vo,
      theme   = theme(plot.title = element_text(size = 13, face = "bold"))
    )
  
  
  cat("Plots built for PC1:", n_pc, "\n")
}


# pam ####

# PAM in PCA space: scatter + Voronoi with PAM MEDOIDS as centroids.
#
# Why pull the cluster from df (not from id.med order):
#   df[[col_name]] holds the HUNGARIAN-ALIGNED PAM labels (from
#   all_clusters_aligned[["pam_k"]]). results_pam[[k]]$id.med returns the
#   medoids in PAM's ORIGINAL (unaligned) numbering. To color each Voronoi
#   cell by the correct cluster we read each medoid's aligned label straight
#   from df: df[[col_name]][med_ids]. Robust to any label permutation.
#
# Requires in scope: results_pam (keyed "2","3",...), var_exp_pca,
#   cluster_colors (named by cluster id, e.g. "1","2",...), deldir, ggplot2,
#   dplyr, patchwork. pca_df must carry PC1/PC2 and a factor column
#   cluster_k{k} aligned row-for-row with the observations PAM was run on.




stopifnot(nrow(pca_df) == nrow(all_clusters_aligned))
for (kv in n_centers) {
  pca_df[[paste0("cluster_k", kv)]] <-
    factor(all_clusters_aligned[[paste0("pam_", kv)]])
}
names(pca_df)   # expected: PC1, PC2, cluster_k2, cluster_k3, ...
# -- Shared helper: PAM medoids projected to PCA, with aligned cluster
# Returns one row per cluster: cl (aligned factor), PC1, PC2 — ordered by cl.
pam_medoids_pca <- function(df, k, col_name) {
  med_ids <- results_pam[[as.character(k)]]$id.med   # medoid ROW indices
  cen <- data.frame(
    cl  = df[[col_name]][med_ids],   # ALIGNED cluster of each medoid (from df)
    PC1 = df$PC1[med_ids],           # medoid coords = exact data-point coords
    PC2 = df$PC2[med_ids],
    stringsAsFactors = FALSE
  )
  cen[order(cen$cl), ]               # tidy order (cosmetic; mapping stays safe)
}


# -- Scatter: points by cluster + medoid markers -----------------------------
plot_pam_pca <- function(df, k) {
  col_name  <- paste0("cluster_k", k)
  pal       <- cluster_colors
  centroids <- pam_medoids_pca(df, k, col_name)
  
  ggplot() +
    geom_point(
      data = df,
      aes(x = PC1, y = PC2, color = .data[[col_name]]),
      alpha = 0.45, size = 0.7
    ) +
    # PAM medoids = centroids (black cross, unambiguous over the points)
    geom_point(
      data = centroids,
      aes(x = PC1, y = PC2),
      shape = 3, size = 3.5, stroke = 1.3, color = "black"
    ) +
    scale_color_manual(values = pal, name = "Cluster") +
    labs(
      title = sprintf("PAM (Gower) in PCA space - K = %d", k),
      x = "PC1",
      y = "PC2"
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


# -- Voronoi: cells built on the PAM medoids, colored by aligned cluster ------
plot_pam_pca_voronoi <- function(df, k) {
  col_name  <- paste0("cluster_k", k)
  pal       <- cluster_colors
  centroids <- pam_medoids_pca(df, k, col_name)   # <- real medoids, aligned cl
  
  pad  <- 0.1
  xlim <- range(df$PC1) + c(-pad, pad) * diff(range(df$PC1))
  ylim <- range(df$PC2) + c(-pad, pad) * diff(range(df$PC2))
  
  vor   <- deldir(centroids$PC1, centroids$PC2, rw = c(xlim, ylim))
  tiles <- tile.list(vor)
  
  # tile i corresponds to centroids row i -> assign its ALIGNED cluster, not i
  vor_df <- do.call(rbind, lapply(seq_along(tiles), function(i) {
    data.frame(x = tiles[[i]]$x, y = tiles[[i]]$y,
               cluster = centroids$cl[i])
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
    # medoid markers
    geom_point(
      data = centroids,
      aes(x = PC1, y = PC2),
      shape = 3, size = 3, color = "black", stroke = 1.2
    ) +
    scale_color_manual(values = pal, name = "Cluster") +
    scale_fill_manual( values = pal, guide = "none") +
    labs(
      title = sprintf("PAM (Gower) + Voronoi in PCA space - K = %d", k),
      x = "PC1",
      y = "PC2"
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


# -- Build the grids ----------------------------------------------------------
plots_pca <- lapply(n_centers, function(k) plot_pam_pca(pca_df, k))
names(plots_pca) <- paste0("k", n_centers)

plots_pca_vor <- lapply(n_centers, function(k) plot_pam_pca_voronoi(pca_df, k))
names(plots_pca_vor) <- paste0("k", n_centers)

vor_caption <-   "Dissimilarity: Gower|  + = PAM medoid"

plot_1_pca <- (plots_pca$k2 | plots_pca$k3) / (plots_pca$k4 | plots_pca$k5) +
  plot_annotation(title = "K-medoids (PAM) clustering - PCA projection",
                  caption = vor_caption,
                  theme = theme(plot.title = element_text(size = 13, face = "bold")))

plot_2_pca <- (plots_pca$k5 | plots_pca$k6) +
  plot_annotation(title = "K-medoids (PAM) clustering - PCA projection",
                  caption = vor_caption,
                  theme = theme(plot.title = element_text(size = 13, face = "bold")))

plot_1_pca_vor <- (plots_pca_vor$k2 | plots_pca_vor$k3) / (plots_pca_vor$k4 | plots_pca_vor$k5) +
  plot_annotation(title = "K-medoids (PAM) - Voronoi tessellation in PCA space",
                  caption = vor_caption,
                  theme = theme(plot.title = element_text(size = 13, face = "bold")))

plot_2_pca_vor <- (plots_pca_vor$k5 | plots_pca_vor$k6) +
  plot_annotation(title = "K-medoids (PAM) - Voronoi tessellation in PCA space",
                  caption = vor_caption,
                  theme = theme(plot.title = element_text(size = 13, face = "bold")))

# Voronoi in pca ####

# UMAP clusters shown in PCA space: scatter + Voronoi.
#
# New dataset = PC1, PC2 + the UMAP cluster label per k.
#   * Voronoi cells are built on the UMAP clusters.
#   * Centroid = MEAN of PC1/PC2 over the observations belonging to that
#     UMAP cluster (as in the dataset).
#   * Each cell is colored by its UMAP cluster, so it holds that cluster's
#     observations (approximately; see note).
#
# NOTE (document in the thesis caption): UMAP k-means clusters need not be
# convex or separable in PCA space, so a few points of a cluster can fall in a
# neighboring Voronoi cell. The centroid is a VISUAL proxy (within-cluster PC
# mean), not a model centroid.

 # 0. Inputs
# pca_scores          : PCA coordinates (>=2 cols), row-aligned with clusters
# all_clusters_aligned: aligned cluster labels incl. km_umap_*
# cluster_colors      : named by cluster id ("1","2",...)
# var_exp_pca         : variance explained per PC (for axis labels)
k_values  <- 2:6
pc_source <- pca_scores     # <- swap to scores_mat / pca_df if PC1/PC2 live there

stopifnot(nrow(pc_source) == nrow(all_clusters_aligned))


# -- 1. Build the new PC1/PC2 + UMAP-cluster dataset
umap_pca_df <- data.frame(PC1 = pc_source[, 1], PC2 = pc_source[, 2])
for (kv in k_values) {
  umap_pca_df[[paste0("cluster_k", kv)]] <-
    factor(all_clusters_aligned[[paste0("km_umap_", kv)]])
}


# -- 2. Helper: visual centroids = mean PC1/PC2 per UMAP cluster --------------
# Returns one row per cluster: .cl (aligned factor), PC1, PC2 — ordered by .cl.
umap_centroids_pca <- function(df, col_name) {
  out <- df %>%
    dplyr::group_by(.cl = .data[[col_name]]) %>%
    dplyr::summarise(PC1 = mean(PC1), PC2 = mean(PC2), .groups = "drop")
  out[order(out$.cl), ]
}


# -- 3. Scatter colored by UMAP cluster + centroid markers -------------------
plot_umap_pca <- function(df, k) {
  col_name  <- paste0("cluster_k", k)
  centroids <- umap_centroids_pca(df, col_name)
  
  ggplot() +
    geom_point(data = df, aes(x = PC1, y = PC2, color = .data[[col_name]]),
               alpha = 0.45, size = 1.2) +
    geom_point(data = centroids, aes(x = PC1, y = PC2),
               shape = 3, size = 3.5, stroke = 1.3, color = "black") +
    scale_color_manual(values = cluster_colors, name = "Cluster") +
    labs(title = sprintf("K-means (UMAP) shown in PCA space - K = %d", k),
         x = "PC1",
         y = "PC2") +
    theme_minimal(base_size = 11) +
    theme(legend.position  = if (k <= 4) "right" else "bottom",
          legend.key.size  = unit(0.4, "cm"),
          panel.grid.minor = element_blank()) +
    guides(color = guide_legend(nrow = if (k <= 8) 2 else 4,
                                override.aes = list(size = 2.5, alpha = 1)))
}


# -- 4. Voronoi on UMAP-cluster means, cells colored by cluster --------------
plot_umap_pca_voronoi <- function(df, k) {
  col_name  <- paste0("cluster_k", k)
  centroids <- umap_centroids_pca(df, col_name)
  
  pad  <- 0.1
  xlim <- range(df$PC1) + c(-pad, pad) * diff(range(df$PC1))
  ylim <- range(df$PC2) + c(-pad, pad) * diff(range(df$PC2))
  
  vor   <- deldir(centroids$PC1, centroids$PC2, rw = c(xlim, ylim))
  tiles <- tile.list(vor)
  
  # tile i corresponds to centroid i -> assign its ACTUAL UMAP cluster, not i
  vor_df <- do.call(rbind, lapply(seq_along(tiles), function(i) {
    data.frame(x = tiles[[i]]$x, y = tiles[[i]]$y,
               cluster = centroids$.cl[i])
  }))
  
  ggplot() +
    geom_polygon(data = vor_df,
                 aes(x = x, y = y, fill = cluster, group = cluster),
                 alpha = 0.25, color = "grey40", linewidth = 0.35) +
    geom_point(data = df, aes(x = PC1, y = PC2, color = .data[[col_name]]),
               alpha = 0.45, size = 0.7) +
    geom_point(data = centroids, aes(x = PC1, y = PC2),
               shape = 3, size = 3, color = "black", stroke = 1.2) +
    scale_color_manual(values = cluster_colors, name = "Cluster") +
    scale_fill_manual( values = cluster_colors, guide = "none") +
    labs(title = sprintf("K-means (UMAP) + Voronoi in PCA space - K = %d", k),
         x = "PC1",
         y = "PC2")+
    theme_minimal(base_size = 11) +
    theme(legend.position  = if (k <= 4) "right" else "bottom",
          legend.key.size  = unit(0.4, "cm"),
          panel.grid.minor = element_blank()) +
    guides(color = guide_legend(nrow = if (k <= 8) 2 else 4,
                                override.aes = list(size = 2.5, alpha = 1)))
}


# -- 5. Build the grids (k = 2..6) -------------------------------------------
plots_pca     <- lapply(k_values, function(k) plot_umap_pca(umap_pca_df, k))
plots_pca_vor <- lapply(k_values, function(k) plot_umap_pca_voronoi(umap_pca_df, k))
names(plots_pca) <- names(plots_pca_vor) <- paste0("k", k_values)

cap <- "Cluster: K-means on UMAP  |  centroid = within-cluster PC1/PC2 mean (visual proxy)"

plot_umap_pca_grid <- wrap_plots(plots_pca, ncol = 3) +
  plot_annotation(title   = "K-means (UMAP) clusters projected in PCA space",
                  caption = cap,
                  theme   = theme(plot.title = element_text(size = 13, face = "bold")))

plot_umap_pca_voronoi_grid <- wrap_plots(plots_pca_vor, ncol = 3) +
  plot_annotation(title   = "K-means (UMAP) - Voronoi tessellation in PCA space",
                  caption = cap,
                  theme   = theme(plot.title = element_text(size = 13, face = "bold")))


# -- 6. Save ------------------------------------------------------------------
out_dir_umap_p <- ensure_dir(DIR_FIG, "07_figures", "umap_in_pca")

ggsave(file.path(out_dir_umap_p, "umap_pca_scatter_grid.pdf"),
       plot_umap_pca_grid,         width = 15, height = 9)
ggsave(file.path(out_dir_umap_p, "umap_pca_voronoi_grid.pdf"),
       plot_umap_pca_voronoi_grid, width = 15, height = 9)


# -- 0. Output directory -------------------------------------------------------
out_dir <- ensure_dir(DIR_FIG, "07_figures", "all_voronoi")


# -- 1. k range + data aliases ------------------------------------------------
k_common <- 2:6
df_pc2 <- pca_dfs[["pc2"]]
df_pc5 <- pca_dfs[["pc5"]]

# Fail fast if any method's df is missing a needed cluster_k column.
for (k in k_common) {
  for (nm in c("df_pc2", "df_pc5", "umap_df", "umap_pca_df", "pca_df")) {
    if (!paste0("cluster_k", k) %in% names(get(nm)))
      stop(sprintf("'%s' is missing column cluster_k%d", nm, k))
  }
}


# -- 2. Helper: 2x2 panel with one shared legend + title ---------------------
combine4 <- function(p_tl, p_tr, p_bl, p_br, title) {
  ((p_tl | p_tr) / (p_bl | p_br)) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = title,
      theme = theme(plot.title = element_text(size = 14, face = "bold"))
    ) &
    theme(legend.position = "bottom")
}


# -- 3. Loop over k: build + save the four images -----------------------------
for (k in k_common) {
  
  # Family A — UMAP native (UMAP1/UMAP2)
  scatter_A <- combine4(
    plot_pca_clusters(df_pc2, k, 2),
    plot_pca_clusters(df_pc5, k, 5),
    plot_umap_clusters(umap_df, k),
    plot_pam_pca(pca_df, k),
    sprintf("Scatter by method  -  k = %d  (UMAP native)", k)
  )
  voronoi_A <- combine4(
    plot_pca_voronoi(df_pc2, k, 2),
    plot_pca_voronoi(df_pc5, k, 5),
    plot_umap_voronoi(umap_df, k),
    plot_pam_pca_voronoi(pca_df, k),
    sprintf("Voronoi by method  -  k = %d  (UMAP native)", k)
  )
  
  # Family B — UMAP projected in PCA (PC1/PC2)
  scatter_B <- combine4(
    plot_pca_clusters(df_pc2, k, 2),
    plot_pca_clusters(df_pc5, k, 5),
    plot_umap_pca(umap_pca_df, k),
    plot_pam_pca(pca_df, k),
    sprintf("Scatter by method  -  k = %d  (UMAP in PCA)", k)
  )
  voronoi_B <- combine4(
    plot_pca_voronoi(df_pc2, k, 2),
    plot_pca_voronoi(df_pc5, k, 5),
    plot_umap_pca_voronoi(umap_pca_df, k),
    plot_pam_pca_voronoi(pca_df, k),
    sprintf("Voronoi by method  -  k = %d  (UMAP in PCA)", k)
  )
  
  # Save (4 files per k)
  ggsave(file.path(out_dir, sprintf("scatter_methods_umapNative_k%d.pdf", k)),
         scatter_A, width = 13, height = 10)
  ggsave(file.path(out_dir, sprintf("voronoi_methods_umapNative_k%d.pdf", k)),
         voronoi_A, width = 13, height = 10)
  ggsave(file.path(out_dir, sprintf("scatter_methods_umapPCA_k%d.pdf", k)),
         scatter_B, width = 13, height = 10)
  ggsave(file.path(out_dir, sprintf("voronoi_methods_umapPCA_k%d.pdf", k)),
         voronoi_B, width = 13, height = 10)
  
  cat("Saved k =", k, "\n")
}

