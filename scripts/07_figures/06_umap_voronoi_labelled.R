# =============================================================================
# scripts/07_figures/06_umap_voronoi_labelled.R
# Purpose : Voronoi tessellations of the 2D UMAP space using the ALIGNED
#           KM-UMAP labels (consistent colours across methods), with and
#           without observations, and with a short description of each
#           cluster written inside its cell.
# Thesis  : Figures 4.6, 4.7 (UMAP panels)
# Input   : results/objects/clusters/umap_df_clusters.rds
#           results/objects/clusters/all_clusters_aligned.rds
# Output  : results/figures/07_figures/umap/umap_voronoi_aligned_grid{,2}.pdf
#           results/figures/07_figures/umap/umap_voronoi_clean_grid{,2}.pdf
# =============================================================================

source(here::here("R", "00_setup.R"))

out_dir_umap_p <- ensure_dir(DIR_FIG, "07_figures", "umap")

umap_df              <- readRDS(file.path(DIR_CLUST, "umap_df_clusters.rds"))
all_clusters_aligned <- readRDS(PATH_CLUSTERS)
stopifnot(nrow(umap_df) == nrow(all_clusters_aligned))

k_values <- 2:6   # the values of k kept in all_clusters_aligned

# -- 1. Replace the raw k-means labels with the aligned ones ------------------
for (k in k_values) {
  umap_df[[paste0("cluster_k", k)]] <- factor(
    all_clusters_aligned[[paste0("km_umap_", k)]]
  )
}

# -- 2. Voronoi with observations ---------------------------------------------
plot_umap_voronoi <- function(df, k) {
  col_name <- paste0("cluster_k", k)
  pal      <- cluster_colors

  # visual centroids: mean UMAP1/UMAP2 per cluster
  centroids <- df %>%
    group_by(.data[[col_name]]) %>%
    summarise(UMAP1 = mean(UMAP1), UMAP2 = mean(UMAP2), .groups = "drop")

  pad  <- 0.1
  xlim <- range(df$UMAP1) + c(-pad, pad) * diff(range(df$UMAP1))
  ylim <- range(df$UMAP2) + c(-pad, pad) * diff(range(df$UMAP2))

  vor   <- deldir(centroids$UMAP1, centroids$UMAP2,
                  rw = c(xlim, ylim))
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
      aes(x = UMAP1, y = UMAP2, color = .data[[col_name]]),
      alpha = 0.45, size = 0.7
    ) +
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

plot_1_vor <- (plots_vor$k2 | plots_vor$k3) / (plots_vor$k4 | plots_vor$k5) +
  plot_annotation(
    title   = "K-means clustering in UMAP 2D space — Voronoi tessellation",
    caption = "UMAP: n_neighbors = 15, min_dist = 0.1, metric = Euclidean ",
    theme   = theme(plot.title = element_text(size = 13, face = "bold"))
  )

plot_2_vor <- (plots_vor$k5 | plots_vor$k6) +
  plot_annotation(
    title   = "K-means clustering in UMAP 2D space — Voronoi tessellation",
    caption = "UMAP: n_neighbors = 15, min_dist = 0.1, metric = Euclidean ",
    theme   = theme(plot.title = element_text(size = 13, face = "bold"))
  )

ggsave(filename = file.path(out_dir_umap_p, "umap_voronoi_aligned_grid.pdf"),
       plot = plot_1_vor, width = 12, height = 9)
ggsave(filename = file.path(out_dir_umap_p, "umap_voronoi_aligned_grid2.pdf"),
       plot = plot_2_vor, width = 12, height = 9)

# -- 3. "Clean" Voronoi: cells only, with cluster descriptions ----------------
plot_umap_voronoi_clean <- function(df, k, cell_text = NULL) {
  # cell_text: named vector mapping cluster label -> string printed in the
  #            cell (e.g. c("1" = "Progressives", "2" = "...")). Use "\n"
  #            for multi-line strings. NULL prints nothing.

  col_name <- paste0("cluster_k", k)

  # visual centroids: mean UMAP1/UMAP2 per cluster (ALIGNED labels)
  centroids <- df %>%
    group_by(.cl = .data[[col_name]]) %>%
    summarise(UMAP1 = mean(UMAP1), UMAP2 = mean(UMAP2), .groups = "drop")

  pad  <- 0.1
  xlim <- range(df$UMAP1) + c(-pad, pad) * diff(range(df$UMAP1))
  ylim <- range(df$UMAP2) + c(-pad, pad) * diff(range(df$UMAP2))

  vor   <- deldir(centroids$UMAP1, centroids$UMAP2, rw = c(xlim, ylim))
  tiles <- tile.list(vor)

  # tile i -> actual cluster label (NOT the index i)
  vor_df <- do.call(rbind, lapply(seq_along(tiles), function(i) {
    data.frame(
      x       = tiles[[i]]$x,
      y       = tiles[[i]]$y,
      cluster = centroids$.cl[i]
    )
  }))
  vor_df$cluster <- factor(vor_df$cluster, levels = levels(factor(centroids$.cl)))

  p <- ggplot() +
    geom_polygon(
      data      = vor_df,
      aes(x = x, y = y, fill = cluster, group = cluster),
      alpha     = 0.45,
      color     = "grey40",
      linewidth = 0.35
    ) +
    scale_fill_manual(values = cluster_colors, name = "Cluster") +
    coord_fixed(xlim = xlim, ylim = ylim, expand = FALSE) +
    labs(
      title = sprintf("Voronoi tessellation (UMAP) — K = %d", k),
      x = "UMAP1", y = "UMAP2"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position  = if (k <= 4) "right" else "bottom",
      legend.key.size  = unit(0.4, "cm"),
      panel.grid.minor = element_blank()
    ) +
    guides(fill = guide_legend(
      nrow         = if (k <= 8) 2 else 4,
      override.aes = list(alpha = 1)
    ))

  # text inside the cells
  if (!is.null(cell_text)) {
    centroids$label <- cell_text[as.character(centroids$.cl)]
    centroids$label[is.na(centroids$label)] <- ""   # cluster without text
    p <- p + geom_text(
      data     = centroids,
      aes(x = UMAP1, y = UMAP2, label = label),
      size     = 3.2,
      fontface = "bold",
      color    = "black",
      lineheight = 0.95
    )
  }

  p
}

# -- 4. Cluster descriptions (from the one-vs-rest profiles, Section 4.2) -----
texts_k2 <- c(
  "1" = "CONSERVATIVE–RELIGIOUS\nDistrust in international organizations",
  "2" = "LIBERALS\nNon-religious\nTrust in international organizations"
)

texts_k3 <- c(
  "1" = "ECONOMIC VULNERABILITY\nCONSERVATIVES\nUnsafe neighborhood\nDistrust in international organizations\nReligion > Science",
  "2" = "SECULAR LIBERALS\nPro-immigrant\nNon-religious\nReligion < Science",
  "3" = "DEVOUT CONSERVATIVES\nActive religiosity\nSafe neighborhood\nReligion > Science"
)

texts_k4 <- c(
  "1" = "ECONOMICALLY VULNERABLE TRADITIONALISTS\nLow religious practice\nUnsafe neighborhood\nApprove of parents beating children (strongly)\nReligion > Science",
  "2" = "ACTIVE RELIGIOSITY\nTrust in international organizations\nPro-science\nEconomically stable\nReligion < Science",
  "3" = "STRONGLY NON RELIGIOUS\nPro-immigrant\nStrongly non-religious\nReligion < Science",
  "4" = "ACTIVE RELIOGIOSITY\nDistrust in international organizations\n\nReligion > Science"
)

texts_k5 <- c(
  "1" = "ECONOMICALLY VULNERABLE TRADITIONALISTS\nReligious\nUnsafe neighborhood\nApprove of parents beating children (strongly)\nReligion < Science",
  "2" = "RELIGIOUS PROGRESSIVES\nReligious\nGender equality\nTrust in international organizations (important)\nReligion < Science",
  "3" = "CONSERVATIVES AND RELIGIOUS\nConservative and religious\nEconomic stability\nDistrust in international organizations",
  "4" = "SECULAR NATIONALISTS\nNon-religious\nAnti-immigration\nAnti-NATO\nNo gender equality\nDistrust in international organizations\nReligion < Science",
  "5" = "SECULAR LIBERALS\nNon-religious\nGender equality\nPro-immigration\nReligion < Science"
)

texts_k6 <- c(
  "1" = "RELIGIOUS LIBERALS\nSafe neighborhood\nTrust in international organizations (pro-NATO)\nReligion < Science",
  "2" = "PROGRESSIVES\nNon-religious\nPro-immigration\nPro-science\nNo respect for authority\nReligion < Science",
  "3" = "RELIGIOUS INSTITUTIONALISTS\nReligious\nConfidence in national institutions\nDistrust in the press\nEconomy > Environment\n Economically stable",
  "4" = "STRONG RELIGIOSITY\nDistrust in economic system\nAnti-NATO (strongly)",
  "5" = "ECONOMICALLY VULNERABLE TRADITIONALIST\nTraditional values\nUnsafe neighborhood\nReligion > Science",
  "6" = "SECULAR TRADITIONALIST\nNon-religious\nAnti-immigration\nNot pro-science"
)

# all texts in one list, indexed by k
cell_texts <- list(
  k2 = texts_k2,
  k3 = texts_k3,
  k4 = texts_k4,
  k5 = texts_k5,
  k6 = texts_k6
)

plots_vor_clean <- lapply(k_values, function(k) {
  txt <- cell_texts[[paste0("k", k)]]   # NULL if missing -> no text
  plot_umap_voronoi_clean(umap_df, k, cell_text = txt)
})
names(plots_vor_clean) <- paste0("k", k_values)

# -- 5. Layout and save -------------------------------------------------------
plot_1_vor_clean <- (plots_vor_clean$k2 | plots_vor_clean$k3) /
  (plots_vor_clean$k4 | plots_vor_clean$k5) +
  plot_annotation(
    title   = "Voronoi tessellation in UMAP 2D space (cells only)",
    caption = "UMAP: n_neighbors = 15, min_dist = 0.1, metric = Euclidean",
    theme   = theme(plot.title = element_text(size = 13, face = "bold"))
  )

plot_2_vor_clean <- (plots_vor_clean$k5 | plots_vor_clean$k6) +
  plot_annotation(
    title   = "Voronoi tessellation in UMAP 2D space (cells only)",
    caption = "UMAP: n_neighbors = 15, min_dist = 0.1, metric = Euclidean",
    theme   = theme(plot.title = element_text(size = 13, face = "bold"))
  )

ggsave(filename = file.path(out_dir_umap_p, "umap_voronoi_clean_grid.pdf"),
       plot = plot_1_vor_clean, width = 12, height = 9)
ggsave(filename = file.path(out_dir_umap_p, "umap_voronoi_clean_grid2.pdf"),
       plot = plot_2_vor_clean, width = 12, height = 9)
