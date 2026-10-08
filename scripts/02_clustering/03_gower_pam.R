# =============================================================================
# scripts/02_clustering/03_gower_pam.R
# Purpose : Gower dissimilarity on the 198 value items, k-medoids (PAM) for
#           k = 2..8 and 16, and a classical-MDS projection for display.
# Thesis  : Sections 3.2.3, 3.3.1; Figures 4.6, 4.7; Tables 4.2, A.3
# Input   : data/processed/data_imputed_selected.rds
# Output  : results/objects/clusters/gower_dist.rds
#           results/objects/{pam_results_allk,pam2_cluster_membership,
#           pam_mds_df_clusters}.rds
#           results/figures/02_clustering/pam/*.pdf
# Note    : all items enter daisy() as unordered factors, so Gower treats them
#           as nominal (see README, Known limitations). PAM is deterministic.
#           The variant on values + political items is
#           scripts/06_robustness/11_gower_pam_with_politics.R, which sets
#           PAM_INPUT / PAM_OUT_OBJ / PAM_OUT_FIG and sources this script.
# =============================================================================

source(here::here("R", "00_setup.R"))

# Defaults (thesis configuration); overridden by the robustness variant.
if (!exists("PAM_INPUT"))   PAM_INPUT   <- file.path(DIR_DATA, "data_imputed_selected.rds")
if (!exists("PAM_OUT_OBJ")) PAM_OUT_OBJ <- DIR_OBJ
if (!exists("PAM_OUT_FIG")) PAM_OUT_FIG <- file.path(DIR_FIG, "02_clustering", "pam")

out_dir_pam   <- ensure_dir(PAM_OUT_OBJ)
out_dir_pam_p <- ensure_dir(PAM_OUT_FIG)

# -- 1. Gower dissimilarity and PAM -------------------------------------------
imputed_data <- readRDS(PAM_INPUT)
gower_dist <- daisy(imputed_data, metric = "gower")

# Thesis configuration (the grid behind the saved objects): k = 2..8, 16.
n_centers <- c(2, 3, 4, 5, 6, 7, 8, 16)
# Reduced variant (k = 2..6 only):
# n_centers <- c(2, 3, 4, 5, 6)

results_pam <- list()
for (k in n_centers) {
  cat("k =", k, "\n")
  results_pam[[as.character(k)]] <- pam(x = gower_dist, k = k, diss = TRUE)
}

# -- 2. 2D coordinates for display (classical MDS on the Gower matrix) --------
mds_fit <- cmdscale(as.matrix(gower_dist), k = 2, eig = TRUE)

mds_df <- data.frame(
  MDS1 = mds_fit$points[, 1],
  MDS2 = mds_fit$points[, 2]
)

# variance explained by the first 2 axes (for the plot captions)
var_exp <- round(mds_fit$eig[1:2] / sum(mds_fit$eig[mds_fit$eig > 0]) * 100, 1)
cat("MDS variance explained: ", var_exp[1], "% /", var_exp[2], "%\n")

# Add the cluster labels for each k
for (k in n_centers) {
  mds_df[[paste0("cluster_k", k)]] <- factor(
    results_pam[[as.character(k)]]$clustering
  )
}

# -- 3. Membership and medoids (k = 2 reference) ------------------------------
pam_2 <- results_pam[["2"]]

cluster_membership <- data.frame(
  id          = seq_len(nrow(imputed_data)),
  cluster_pam = factor(pam_2$clustering)
)

table(cluster_membership$cluster_pam)
prop.table(table(cluster_membership$cluster_pam))
pam_2$id.med
imputed_data[pam_2$id.med, ]

# -- 4. Scatter plots in MDS space --------------------------------------------
plot_pam_clusters <- function(df, k) {
  col_name <- paste0("cluster_k", k)
  pal      <- cluster_colors

  ggplot(df, aes(x = MDS1, y = MDS2, color = .data[[col_name]])) +
    geom_point(alpha = 0.45, size = 0.7) +
    scale_color_manual(values = pal, name = "Cluster") +
    labs(
      title    = sprintf("PAM (Gower) — K = %d", k),
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
      nrow         = if (k <= 8) 2 else 4,
      override.aes = list(size = 2.5, alpha = 1)
    ))
}

plots <- lapply(n_centers, function(k) plot_pam_clusters(mds_df, k))
names(plots) <- paste0("k", n_centers)

# Layout (same scheme as UMAP)
plot_1 <- (plots$k2 | plots$k3) / (plots$k4 | plots$k5) +
  plot_annotation(
    title   = "K-medoids (PAM) clustering — MDS projection",
    caption = sprintf("Dissimilarity: Gower  |  MDS1 %.1f%% / MDS2 %.1f%%",
                      var_exp[1], var_exp[2]),
    theme   = theme(plot.title = element_text(size = 13, face = "bold"))
  )

plot_2 <- (plots$k5 | plots$k6) / (plots$k7 | plots$k8) +
  plot_annotation(
    title   = "K-medoids (PAM) clustering — MDS projection",
    caption = sprintf("Dissimilarity: Gower  |  MDS1 %.1f%% / MDS2 %.1f%%",
                      var_exp[1], var_exp[2]),
    theme   = theme(plot.title = element_text(size = 13, face = "bold"))
  )

# -- 5. Voronoi tessellation in MDS space -------------------------------------
plot_pam_voronoi <- function(df, k) {
  col_name <- paste0("cluster_k", k)
  pal      <- cluster_colors

  # visual centroids: mean MDS1/MDS2 per cluster
  centroids <- df %>%
    group_by(.data[[col_name]]) %>%
    summarise(MDS1 = mean(MDS1), MDS2 = mean(MDS2), .groups = "drop")

  # Voronoi tessellation of the visual centroids, window covering all points
  pad  <- 0.1
  xlim <- range(df$MDS1) + c(-pad, pad) * diff(range(df$MDS1))
  ylim <- range(df$MDS2) + c(-pad, pad) * diff(range(df$MDS2))

  vor  <- deldir(centroids$MDS1, centroids$MDS2,
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
      data    = vor_df,
      aes(x = x, y = y, fill = cluster, group = cluster),
      alpha   = 0.25,
      color   = "grey40",
      linewidth = 0.35
    ) +
    # observations
    geom_point(
      data = df,
      aes(x = MDS1, y = MDS2, color = .data[[col_name]]),
      alpha = 0.45, size = 0.7
    ) +
    # visual centroids
    geom_point(
      data  = centroids,
      aes(x = MDS1, y = MDS2),
      shape = 3, size = 3, color = "black", stroke = 1.2
    ) +
    scale_color_manual(values = pal, name = "Cluster") +
    scale_fill_manual( values = pal, guide = "none") +
    labs(
      title    = sprintf("PAM (Gower) + Voronoi — K = %d", k),
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
      nrow         = if (k <= 8) 2 else 4,
      override.aes = list(size = 2.5, alpha = 1)
    ))
}

plots_vor <- lapply(n_centers, function(k) plot_pam_voronoi(mds_df, k))
names(plots_vor) <- paste0("k", n_centers)

plot_1_vor <- (plots_vor$k2 | plots_vor$k3) / (plots_vor$k4 | plots_vor$k5) +
  plot_annotation(
    title   = "K-medoids (PAM) — Voronoi tessellation in MDS space",
    caption = sprintf("Dissimilarity: Gower  |  MDS1 %.1f%% / MDS2 %.1f%%  |  + = visual centroid",
                      var_exp[1], var_exp[2]),
    theme   = theme(plot.title = element_text(size = 13, face = "bold"))
  )

plot_2_vor <- (plots_vor$k5 | plots_vor$k6) / (plots_vor$k7 | plots_vor$k8) +
  plot_annotation(
    title   = "K-medoids (PAM) — Voronoi tessellation in MDS space",
    caption = sprintf("Dissimilarity: Gower  |  MDS1 %.1f%% / MDS2 %.1f%%  |  + = visual centroid",
                      var_exp[1], var_exp[2]),
    theme   = theme(plot.title = element_text(size = 13, face = "bold"))
  )

# -- 6. Save ------------------------------------------------------------------
saveRDS(results_pam,        file = file.path(out_dir_pam, "pam_results_allk.rds"))
saveRDS(gower_dist,         file = file.path(if (PAM_OUT_OBJ == DIR_OBJ) ensure_dir(DIR_CLUST) else out_dir_pam, "gower_dist.rds"))
saveRDS(cluster_membership, file = file.path(out_dir_pam, "pam2_cluster_membership.rds"))
saveRDS(mds_df,             file = file.path(out_dir_pam, "pam_mds_df_clusters.rds"))

ggsave(filename = file.path(out_dir_pam_p, "pam_scatter_grid.pdf"),
       plot = plot_1, width = 12, height = 9)
ggsave(filename = file.path(out_dir_pam_p, "pam_scatter_grid2.pdf"),
       plot = plot_2, width = 12, height = 9)
ggsave(filename = file.path(out_dir_pam_p, "pam_voronoi_grid.pdf"),
       plot = plot_1_vor, width = 12, height = 9)
ggsave(filename = file.path(out_dir_pam_p, "pam_voronoi_grid2.pdf"),
       plot = plot_2_vor, width = 12, height = 9)
