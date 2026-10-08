# =============================================================================
# scripts/02_clustering/02_umap_kmeans.R
# Purpose : UMAP embedding (2D, 3D, 5D) of the standardised value items, then
#           k-means on the 2D embedding (KM-UMAP) for k = 2..8 and 16.
# Thesis  : Sections 3.2.2, 3.3.1; Figures 4.3, 4.6, 4.7; Tables 4.2, A.3
# Input   : data/processed/data_imputed_selected.rds
# Output  : results/objects/umap_kmeans_objects.rds
#           results/objects/clusters/umap_df_clusters.rds (2D coordinates + raw
#             k-means labels; the aligned labels are built in 07_align_clusters.R)
#           results/objects/UMAP/umapd_5d.rds
#           results/figures/02_clustering/umap/*.pdf
# Seeds   : set.seed(42); UMAP random_state = 42; kmeans(nstart = 50)
# Note    : the HDBSCAN scripts use a separate 2D/5D embedding computed in
#           05_hdbscan.R (results/objects/UMAP/umap_2d.rds, umap_5d.rds). The two
#           2D layouts are close but not identical (see README, Known limitations).
#           The sensitivity variant (n_neighbors = 50, min_dist = 1e-6, 2D only,
#           k = 2..6) is scripts/06_robustness/10_umap_sensitivity.R, which sets
#           the UMAP_* parameters below and sources this script.
# =============================================================================

source(here::here("R", "00_setup.R"))

# Defaults (thesis configuration); overridden by the sensitivity variant.
if (!exists("UMAP_N_NEIGHBORS")) UMAP_N_NEIGHBORS <- 15
if (!exists("UMAP_MIN_DIST"))    UMAP_MIN_DIST    <- 0.1
if (!exists("UMAP_EXTRA_DIMS"))  UMAP_EXTRA_DIMS  <- TRUE   # also fit 3D and 5D
if (!exists("UMAP_K_VALUES"))    UMAP_K_VALUES    <- c(2, 3, 4, 5, 6, 7, 8, 16)
if (!exists("UMAP_OUT_OBJ"))     UMAP_OUT_OBJ     <- DIR_OBJ
if (!exists("UMAP_OUT_CLUST"))   UMAP_OUT_CLUST   <- DIR_CLUST
if (!exists("UMAP_OUT_FIG"))     UMAP_OUT_FIG     <- file.path(DIR_FIG, "02_clustering", "umap")

out_dir_umap   <- ensure_dir(UMAP_OUT_OBJ)
out_dir_umap_p <- ensure_dir(UMAP_OUT_FIG)

# -- 1. The data matrix has to be numeric -------------------------------------
X <- readRDS(file.path(DIR_DATA, "data_imputed_selected.rds")) %>%
  mutate(across(where(is.factor), ~ as.numeric(as.character(.)))) %>%
  scale()                        # standardization

str(X)

# -- 2. UMAP configuration ----------------------------------------------------
umap_config <- umap.defaults

umap_config$n_neighbors   <- UMAP_N_NEIGHBORS   # thesis: 15
umap_config$min_dist      <- UMAP_MIN_DIST      # thesis: 0.1
umap_config$metric        <- "euclidean"
umap_config$n_epochs      <- 500    # optimisation iterations
umap_config$random_state  <- 42     # reproducibility

# -- 3. Fit UMAP --------------------------------------------------------------
set.seed(42)
umap_fit   <- umap(X, config = umap_config, n_components = 2)

# -- 4. Data frames of results ------------------------------------------------
umap_df <- data.frame(
  UMAP1   = umap_fit$layout[, 1],
  UMAP2   = umap_fit$layout[, 2]
)

if (UMAP_EXTRA_DIMS) {
  umap_fit3d <- umap(X, config = umap_config, n_components = 3)
  umap_fit5d <- umap(X, config = umap_config, n_components = 5)
  str(umap_fit3d)

  umap_df3d <- data.frame(
    UMAP1   = umap_fit3d$layout[, 1],
    UMAP2   = umap_fit3d$layout[, 2],
    UMAP3   = umap_fit3d$layout[, 3]
  )

  umap_df5d <- data.frame(
    UMAP1   = umap_fit5d$layout[, 1],
    UMAP2   = umap_fit5d$layout[, 2],
    UMAP3   = umap_fit5d$layout[, 3],
    UMAP4   = umap_fit5d$layout[, 4],
    UMAP5   = umap_fit5d$layout[, 5]
  )

  saveRDS(umap_df5d, file.path(ensure_dir(UMAP_OUT_OBJ, "UMAP"), "umapd_5d.rds"))
}

# -- 5. K-means on the 2D UMAP space ------------------------------------------
# Thesis configuration (the grid behind the saved labels): k = 2..8, 16.
set.seed(42)
k_values <- UMAP_K_VALUES

kmeans_umap <- lapply(k_values, function(k) {
  kmeans(umap_df, centers = k, nstart = 50, iter.max = 100)
})
names(kmeans_umap) <- paste0("k", k_values)

# Add the labels to the data frame
for (k in k_values) {
  umap_df[[paste0("cluster_k", k)]] <- factor(kmeans_umap[[paste0("k", k)]]$cluster)
}

# -- 6. Scatter plots ---------------------------------------------------------
plot_umap_clusters <- function(df, k) {
  col_name <- paste0("cluster_k", k)

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

# -- 7. Layout ----------------------------------------------------------------
umap_caption <- sprintf("UMAP: n_neighbors = %s, min_dist = %s, metric = Euclidean",
                        UMAP_N_NEIGHBORS, UMAP_MIN_DIST)
has_grid2    <- all(c(5, 6, 7, 8) %in% k_values)   # second grid needs k = 5..8

plot_1 <- (plots$k2 | plots$k3) / (plots$k4 | plots$k5) +
  plot_annotation(
    title   = "K-means clustering in UMAP 2D space",
    caption = umap_caption,
    theme   = theme(plot.title = element_text(size = 13, face = "bold"))
  )

if (has_grid2) plot_2 <- (plots$k6 | plots$k7) / (plots$k7 | plots$k8) +
  plot_annotation(
    title   = "K-means clustering in UMAP 2D space",
    caption = umap_caption,
    theme   = theme(plot.title = element_text(size = 13, face = "bold"))
  )

# -- 8. Voronoi tessellation in UMAP space (raw k-means labels) ---------------
# The version with aligned labels and cluster descriptions is in
# scripts/07_figures/06_umap_voronoi_labelled.R.
plot_umap_voronoi <- function(df, k) {
  col_name <- paste0("cluster_k", k)
  pal      <- cluster_colors

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

# -- 9. Voronoi layout --------------------------------------------------------
plot_1_vor <- (plots_vor$k2 | plots_vor$k3) / (plots_vor$k4 | plots_vor$k5) +
  plot_annotation(
    title   = "K-means clustering in UMAP 2D space — Voronoi tessellation",
    caption = umap_caption,
    theme   = theme(plot.title = element_text(size = 13, face = "bold"))
  )

if (has_grid2) plot_2_vor <- (plots_vor$k5 | plots_vor$k6) / (plots_vor$k7 | plots_vor$k8) +
  plot_annotation(
    title   = "K-means clustering in UMAP 2D space — Voronoi tessellation",
    caption = umap_caption,
    theme   = theme(plot.title = element_text(size = 13, face = "bold"))
  )

# -- 10. Save -----------------------------------------------------------------
# full k-means objects (centroids, withinss, ...) for each k
saveRDS(kmeans_umap, file = file.path(out_dir_umap, "umap_kmeans_objects.rds"))

# UMAP coordinates + cluster_k2 ... cluster_k16 columns
saveRDS(umap_df, file = file.path(ensure_dir(UMAP_OUT_CLUST), "umap_df_clusters.rds"))

ggsave(filename = file.path(out_dir_umap_p, "umap_kmeans_grid.pdf"),
       plot = plot_1, width = 12, height = 9)
ggsave(filename = file.path(out_dir_umap_p, "umap_voronoi_grid.pdf"),
       plot = plot_1_vor, width = 12, height = 9)
if (has_grid2) {
  ggsave(filename = file.path(out_dir_umap_p, "umap_kmeans_grid2.pdf"),
         plot = plot_2, width = 12, height = 9)
  ggsave(filename = file.path(out_dir_umap_p, "umap_voronoi_grid2.pdf"),
         plot = plot_2_vor, width = 12, height = 9)
}
