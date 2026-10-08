# =============================================================================
# scripts/02_clustering/05_hdbscan.R
# Purpose : HDBSCAN on UMAP embeddings of dimension d = 2, 5, 20, 50.
#   1. UMAP dimensionality reduction to each d.
#   2. HDBSCAN on each embedding. The number of clusters is NOT imposed: it
#      emerges from the density structure. Granularity is governed only by
#      minPts (core-distance k / min cluster size). Points that belong to no
#      dense region are labelled noise (cluster 0).
#   3. Consolidated saves (compute once, reload afterwards).
#   4. Plots, each saved right after it is built.
#   A minPts grid produces a family of density-based partitions, each with
#   its own emergent number of clusters; there is no k to select by hand.
# Thesis  : Sections 3.3.2, 4.1.2; Figure 4.4; Appendix A.4 (Tables A.8, A.9)
# Input   : data/processed/data_imputed_selected.rds
# Output  : results/objects/HDBSCAN/{umap_embeddings,hdbscan_objects,
#           hdbscan_partitions,hdbscan_minpts_sweep}.rds
#           results/objects/UMAP/{umap_2d,umap_5d}.rds (used by 06_hdbscan_vs_kmeans.R)
#           results/figures/02_clustering/hdbscan/*.pdf
# Seeds   : seed = 42 (set before each UMAP fit; UMAP random_state = 42)
# =============================================================================

source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages(library(dbscan))

# -- 0. Paths -----------------------------------------------------------------
out_dir   <- ensure_dir(DIR_OBJ, "HDBSCAN")
out_dir_p <- ensure_dir(DIR_FIG, "02_clustering", "hdbscan")

# recompute switch: set FALSE to reload saved objects instead of recomputing
recompute <- TRUE

# -- 0b. Hyper-parameters (tune here only) ------------------------------------
# Thesis configuration: UMAP output dimensionalities 2, 5, 20, 50.
umap_dims   <- c(2, 5, 20, 50)
# Variant used in an intermediate version of this script:
# umap_dims <- c(2, 5, 10, 20)
minPts      <- 25                # PRIMARY run: core-dist k / min cluster size
                                 #   larger -> fewer, more robust clusters
minPts_grid <- c(5, 10, 15, 20, 25, 30, 40, 50)  # family of emergent runs
seed        <- 42

# UMAP config: kept close to the existing km_umap pipeline.
#   For HDBSCAN, min_dist ~ 0 (tighter clumps) is often preferable; left at
#   0.1 for continuity with the km_umap embeddings — tune if needed.
umap_config <- umap.defaults
umap_config$n_neighbors  <- 15
umap_config$min_dist     <- 0.1
umap_config$metric       <- "euclidean"
umap_config$n_epochs     <- 500
umap_config$random_state <- seed
# For large d (20, 50) the naive R umap can be slow; uwot::umap is a faster
# drop-in alternative if runtime becomes a problem.

# ── 0c. Palette helper (extends cluster_colors; grey = noise) ─────── #####
noise_col <- "grey75"
hdb_palette <- function(n) {
  if (exists("cluster_colors")) {
    base <- unname(cluster_colors)
    if (length(base) >= n) base[seq_len(n)]
    else grDevices::colorRampPalette(base)(n)
  } else {
    grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))(n)
  }
}

# ── 1. Data matrix (numeric + standardized) ───────────────────────── #####
X <- readRDS(file.path(DIR_DATA, "data_imputed_selected.rds")) %>%
  dplyr::mutate(dplyr::across(where(is.factor),
                              ~ as.numeric(as.character(.)))) %>%
  scale()
stopifnot(is.matrix(X), all(is.finite(X)))
n_obs <- nrow(X)

# ── 2. UMAP embeddings at each target dimensionality ──────────────── #####
emb_file <- file.path(out_dir, "umap_embeddings.rds")

if (recompute || !file.exists(emb_file)) {
  umap_embeddings <- lapply(umap_dims, function(d) {
    set.seed(seed)
    cfg <- umap_config
    cfg$n_components <- d
    umap(X, config = cfg)$layout          # n_obs x d matrix
  })
  names(umap_embeddings) <- paste0("dim", umap_dims)
  saveRDS(umap_embeddings, emb_file)
} else {
  umap_embeddings <- readRDS(emb_file)
}

# 2D canvas for ALL scatter plots (clusters may be computed in higher dims,
# but are always visualized on this shared 2D layout).
coords2d <- as.data.frame(umap_embeddings$dim2)
colnames(coords2d) <- c("UMAP1", "UMAP2")

# 2D and 5D embeddings as stand-alone matrices (read by 06_hdbscan_vs_kmeans.R
# and the HDBSCAN diagnostics in 06_robustness/).
saveRDS(umap_embeddings$dim2, file.path(ensure_dir(DIR_OBJ, "UMAP"), "umap_2d.rds"))
saveRDS(umap_embeddings$dim5, file.path(DIR_OBJ, "UMAP", "umap_5d.rds"))

# ── 3. HDBSCAN: primary objects + emergent partitions across minPts ─ #####
hdb_file   <- file.path(out_dir, "hdbscan_objects.rds")
part_file  <- file.path(out_dir, "hdbscan_partitions.rds")
sweep_file <- file.path(out_dir, "hdbscan_minpts_sweep.rds")

if (recompute || !file.exists(part_file)) {

  # (a) full native objects at the PRIMARY minPts (used by the diagnostics:
  #     condensed tree, GLOSH, membership probability, stability)
  hdbscan_objects <- lapply(umap_embeddings, function(emb) {
    dbscan::hdbscan(emb, minPts = minPts)
  })
  names(hdbscan_objects) <- names(umap_embeddings)

  # (b) family of emergent partitions over the minPts grid. Number of
  #     clusters is whatever the density yields at each minPts. Labels and
  #     a summary (n_clusters, noise fraction) are collected in one pass.
  #       hdb_umap{d}_mp{minPts} : density-based labels, 0 = noise
  hdbscan_partitions <- data.frame(row_id = seq_len(n_obs))
  summary_rows <- list()
  for (d in umap_dims) {
    emb <- umap_embeddings[[paste0("dim", d)]]
    for (mp in minPts_grid) {
      res <- dbscan::hdbscan(emb, minPts = mp)
      hdbscan_partitions[[paste0("hdb_umap", d, "_mp", mp)]] <- res$cluster
      summary_rows[[length(summary_rows) + 1L]] <- data.frame(
        dim        = d,
        minPts     = mp,
        n_clusters = length(unique(res$cluster[res$cluster != 0])),
        noise_frac = mean(res$cluster == 0)
      )
    }
  }
  minpts_sweep <- do.call(rbind, summary_rows)

  saveRDS(hdbscan_objects,    hdb_file)
  saveRDS(hdbscan_partitions, part_file)
  saveRDS(minpts_sweep,       sweep_file)
} else {
  hdbscan_objects    <- readRDS(hdb_file)
  hdbscan_partitions <- readRDS(part_file)
  minpts_sweep       <- readRDS(sweep_file)
}

# ═══════════════════════════════════════════════════════════════════════
#  PLOTS  — each figure is saved immediately after being built
# ═══════════════════════════════════════════════════════════════════════

# reusable 2D scatter coloured by a discrete cluster vector (0 -> noise grey)
plot_hdb_scatter <- function(coords, labels, title,
                             subtitle = NULL, drop_noise_legend = TRUE) {
  df         <- coords
  df$cluster <- factor(labels)
  cl_levels  <- levels(df$cluster)
  has_noise  <- "0" %in% cl_levels
  non_noise  <- setdiff(cl_levels, "0")
  pal        <- hdb_palette(length(non_noise))
  names(pal) <- non_noise
  if (has_noise) pal <- c("0" = noise_col, pal)

  ggplot(df, aes(UMAP1, UMAP2, color = cluster)) +
    geom_point(alpha = 0.45, size = 0.7) +
    scale_color_manual(
      values = pal, name = "Cluster",
      breaks = if (drop_noise_legend) non_noise else cl_levels) +
    labs(title = title, subtitle = subtitle, x = "UMAP1", y = "UMAP2") +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          legend.key.size  = unit(0.4, "cm")) +
    guides(color = guide_legend(override.aes = list(size = 2.5, alpha = 1)))
}

# ── P1. Native HDBSCAN clustering at primary minPts (one per dim) ──── #####
plots_native <- lapply(umap_dims, function(d) {
  res <- hdbscan_objects[[paste0("dim", d)]]
  ncl <- length(unique(res$cluster[res$cluster != 0]))
  plot_hdb_scatter(
    coords2d, res$cluster,
    title    = sprintf("HDBSCAN — UMAP %dd", d),
    subtitle = sprintf("minPts = %d | %d emergent clusters | %.1f%% noise",
                       minPts, ncl, 100 * mean(res$cluster == 0)))
})
names(plots_native) <- paste0("dim", umap_dims)

plot_native_grid <- (plots_native$dim2 | plots_native$dim5) /
                    (plots_native$dim20 | plots_native$dim50) +
  plot_annotation(
    title = "HDBSCAN clustering (projected on 2D UMAP)",
    theme = theme(plot.title = element_text(size = 13, face = "bold")))

ggsave(file.path(out_dir_p, "hdbscan_native_grid.pdf"),
       plot_native_grid, width = 12, height = 9)

# ── P2. Emergent clustering across the minPts grid (one grid per dim)  #####
for (d in umap_dims) {
  ps <- lapply(minPts_grid, function(mp) {
    lab <- hdbscan_partitions[[paste0("hdb_umap", d, "_mp", mp)]]
    ncl <- length(unique(lab[lab != 0]))
    plot_hdb_scatter(
      coords2d, lab,
      title = sprintf("minPts = %d  (%d clusters, %.0f%% noise)",
                      mp, ncl, 100 * mean(lab == 0)))
  })
  grid <- wrap_plots(ps, ncol = 3) +
    plot_annotation(
      title = sprintf("HDBSCAN over minPts — UMAP %dd (projected on 2D)", d),
      theme = theme(plot.title = element_text(size = 13, face = "bold")))
  ggsave(file.path(out_dir_p, sprintf("hdbscan_minpts_dim%d.pdf", d)),
         grid, width = 13, height = 10)
}

# ── P3. Condensed cluster tree (the signature HDBSCAN plot) ───────── #####
pdf(file.path(out_dir_p, "hdbscan_condensed_tree.pdf"), width = 12, height = 9)
op <- par(mfrow = c(2, 2))
for (d in umap_dims) {
  plot(hdbscan_objects[[paste0("dim", d)]], show_flat = TRUE)
  title(main = sprintf("Condensed tree — UMAP %dd (minPts = %d)", d, minPts))
}
par(op); dev.off()

# ── P4. Cluster stability / persistence scores ────────────────────── #####
stab_df <- do.call(rbind, lapply(umap_dims, function(d) {
  sc <- hdbscan_objects[[paste0("dim", d)]]$cluster_scores
  if (length(sc) == 0) return(NULL)
  data.frame(dim = paste0("UMAP ", d, "d"),
             cluster = names(sc), stability = as.numeric(sc))
}))

p_stab <- ggplot(stab_df, aes(cluster, stability, fill = cluster)) +
  geom_col(show.legend = FALSE) +
  facet_wrap(~ dim, scales = "free_x") +
  scale_fill_manual(values = hdb_palette(length(unique(stab_df$cluster)))) +
  labs(title = "HDBSCAN cluster stability (persistence)",
       x = "Cluster", y = "Stability score") +
  theme_minimal(base_size = 11)

ggsave(file.path(out_dir_p, "hdbscan_cluster_stability.pdf"),
       p_stab, width = 11, height = 8)

# ── P5. GLOSH outlier scores ──────────────────────────────────────── #####
# (a) distribution
glosh_df <- do.call(rbind, lapply(umap_dims, function(d) {
  data.frame(dim   = paste0("UMAP ", d, "d"),
             glosh = hdbscan_objects[[paste0("dim", d)]]$outlier_scores)
}))
p_glosh_hist <- ggplot(glosh_df, aes(glosh)) +
  geom_histogram(bins = 40, fill = "steelblue", color = "white") +
  facet_wrap(~ dim, scales = "free_y") +
  labs(title = "GLOSH outlier score distribution", x = "GLOSH", y = "Count") +
  theme_minimal(base_size = 11)
ggsave(file.path(out_dir_p, "hdbscan_glosh_hist.pdf"),
       p_glosh_hist, width = 11, height = 8)

# (b) spatial: 2D scatter coloured by GLOSH
plots_glosh <- lapply(umap_dims, function(d) {
  df <- coords2d
  df$glosh <- hdbscan_objects[[paste0("dim", d)]]$outlier_scores
  ggplot(df, aes(UMAP1, UMAP2, color = glosh)) +
    geom_point(size = 0.7) +
    scale_color_viridis_c(option = "inferno", name = "GLOSH") +
    labs(title = sprintf("GLOSH — UMAP %dd", d), x = "UMAP1", y = "UMAP2") +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank())
})
grid_glosh <- wrap_plots(plots_glosh, ncol = 2) +
  plot_annotation(title = "GLOSH outlier scores on 2D UMAP",
                  theme = theme(plot.title = element_text(size = 13, face = "bold")))
ggsave(file.path(out_dir_p, "hdbscan_glosh_scatter.pdf"),
       grid_glosh, width = 12, height = 9)

# ── P6. Membership probability (assignment confidence) ────────────── #####
plots_prob <- lapply(umap_dims, function(d) {
  df <- coords2d
  df$prob <- hdbscan_objects[[paste0("dim", d)]]$membership_prob
  ggplot(df, aes(UMAP1, UMAP2, color = prob)) +
    geom_point(size = 0.7) +
    scale_color_viridis_c(name = "P(member)", limits = c(0, 1)) +
    labs(title = sprintf("Membership probability — UMAP %dd", d),
         x = "UMAP1", y = "UMAP2") +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank())
})
grid_prob <- wrap_plots(plots_prob, ncol = 2) +
  plot_annotation(title = "HDBSCAN membership probabilities on 2D UMAP",
                  theme = theme(plot.title = element_text(size = 13, face = "bold")))
ggsave(file.path(out_dir_p, "hdbscan_membership_prob.pdf"),
       grid_prob, width = 12, height = 9)

# ── P7. k-NN distance plot (minPts / core-distance diagnostic) ────── #####
pdf(file.path(out_dir_p, "hdbscan_knndist.pdf"), width = 12, height = 9)
op <- par(mfrow = c(2, 2))
for (d in umap_dims) {
  dbscan::kNNdistplot(umap_embeddings[[paste0("dim", d)]], k = minPts)
  title(main = sprintf("kNN distance (k = %d) — UMAP %dd", minPts, d))
}
par(op); dev.off()

# ── P8. n_clusters & noise fraction vs minPts (emergent structure) ── #####
p_sweep <- ggplot(minpts_sweep, aes(minPts, n_clusters, color = factor(dim))) +
  geom_line() + geom_point() +
  scale_color_manual(values = hdb_palette(length(unique(minpts_sweep$dim))),
                     name = "UMAP dim") +
  labs(title = "Emergent number of HDBSCAN clusters vs minPts",
       x = "minPts", y = "n clusters (excl. noise)") +
  theme_minimal(base_size = 11)
ggsave(file.path(out_dir_p, "hdbscan_nclusters_vs_minpts.pdf"),
       p_sweep, width = 9, height = 6)

p_noise <- ggplot(minpts_sweep, aes(minPts, noise_frac, color = factor(dim))) +
  geom_line() + geom_point() +
  scale_color_manual(values = hdb_palette(length(unique(minpts_sweep$dim))),
                     name = "UMAP dim") +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Noise fraction vs minPts",
       x = "minPts", y = "% points labelled noise") +
  theme_minimal(base_size = 11)
ggsave(file.path(out_dir_p, "hdbscan_noise_vs_minpts.pdf"),
       p_noise, width = 9, height = 6)

# ═══════════════════════════════════════════════════════════════════════
#  Saved objects (results/objects/HDBSCAN):
#    umap_embeddings.rds        list of layouts, one per dim (2/5/20/50)
#    hdbscan_objects.rds        native hdbscan objects at primary minPts
#                               (hc, GLOSH, membership_prob, cluster_scores)
#    hdbscan_partitions.rds     emergent labels per dim x minPts (0 = noise)
#    hdbscan_minpts_sweep.rds   n_clusters / noise summary per dim x minPts
#  Set recompute <- FALSE to reload without recomputing. Every plot object
#  (plot_native_grid, p_stab, p_sweep, ...) can be printed in RStudio.
