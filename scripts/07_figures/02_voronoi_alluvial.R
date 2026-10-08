# =============================================================================
# scripts/07_figures/02_voronoi_alluvial.R
# Purpose : For pairs of configurations and every k: Voronoi scatter of each
#           partition in its own space plus an alluvial plot of the flows
#           between them (Hungarian alignment); discordant respondents in
#           black. Extra blocks for PAM vs KM-UMAP / KM-PC2 / KM-PC5.
# Thesis  : Figure 4.6
# Input   : results/objects/clusters/{all_clusters,umap_df_clusters}.rds
#           results/objects/{pca_scores,pam_mds_df_clusters}.rds
# Output  : results/figures/07_figures/voronoi/<pair>/*.pdf
# =============================================================================

source(here::here("R", "00_setup.R"))

# -- 0. Load --------------------------------------------------------------------
out_dir_obj <- DIR_OBJ
mds_df <- readRDS(file.path(out_dir_obj, "pam_mds_df_clusters.rds"))
df_mds <- data.frame(MDS1 = mds_df$MDS1, MDS2 = mds_df$MDS2)
all_clusters <- readRDS(file.path(DIR_CLUST, "all_clusters.rds"))
pca_scores   <- readRDS(file.path(out_dir_obj, "pca_scores.rds"))
umap_df      <- readRDS(file.path(DIR_CLUST, "umap_df_clusters.rds"))

df_pca  <- data.frame(PC1   = pca_scores[, 1], PC2   = pca_scores[, 2])
df_umap <- data.frame(UMAP1 = umap_df$UMAP1,   UMAP2 = umap_df$UMAP2)

k_values <- c(2, 3, 4, 5, 6, 7, 8, 16)

method_pairs <- list(
  list(ref = "km_pc5",  oth = "km_umap",  ref_lab = "KMeans-5PC",
       oth_lab = "KM-UMAP",   left_space = "pca",  right_space = "umap"),
  list(ref = "km_pc5",  oth = "pam",      ref_lab = "KMeans-5PC",
       oth_lab = "PAM",       left_space = "mds",  right_space = "pca"),
  list(ref = "km_pc5",  oth = "km_pc2",   ref_lab = "KMeans-5PC",
       oth_lab = "KM-PC2",    left_space = "pca",  right_space = "umap"),
  list(ref = "km_pc5",  oth = "km_pc15",  ref_lab = "KMeans-5PC",
       oth_lab = "KM-PC15",   left_space = "pca",  right_space = "umap"),
  list(ref = "km_pc5",  oth = "km_pc75",  ref_lab = "KMeans-5PC",
       oth_lab = "KM-PC75",   left_space = "pca",  right_space = "umap"),
  list(ref = "km_pc5",  oth = "km_pc198", ref_lab = "KMeans-5PC",
       oth_lab = "KM-PC198",  left_space = "pca",  right_space = "umap"),
  list(ref = "km_pc15", oth = "pam",      ref_lab = "KMeans-15PC",
       oth_lab = "PAM",       left_space = "mds",  right_space = "pca"),
  list(ref = "km_pc15", oth = "km_umap",  ref_lab = "KMeans-15PC",
       oth_lab = "KM-UMAP",   left_space = "pca",  right_space = "umap"),
  list(ref = "km_pc75", oth = "pam",      ref_lab = "KMeans-75PC",
       oth_lab = "PAM",       left_space = "mds",  right_space = "pca"),
  list(ref = "km_pc75", oth = "km_umap",  ref_lab = "KMeans-75PC",
       oth_lab = "KM-UMAP",   left_space = "pca",  right_space = "umap"),
  # ── new: km_pc2 vs all ──────────────────────────────────────────
  list(ref = "km_pc2",  oth = "km_umap",  ref_lab = "KMeans-2PC",
       oth_lab = "KM-UMAP",   left_space = "pca",  right_space = "umap"),
  list(ref = "km_pc2",  oth = "pam",      ref_lab = "KMeans-2PC",
       oth_lab = "PAM",       left_space = "mds",  right_space = "pca"),
  list(ref = "km_pc2",  oth = "km_pc5",   ref_lab = "KMeans-2PC",
       oth_lab = "KM-5PC",    left_space = "pca",  right_space = "umap"),
  list(ref = "km_pc2",  oth = "km_pc15",  ref_lab = "KMeans-2PC",
       oth_lab = "KM-PC15",   left_space = "pca",  right_space = "umap"),
  list(ref = "km_pc2",  oth = "km_pc75",  ref_lab = "KMeans-2PC",
       oth_lab = "KM-PC75",   left_space = "pca",  right_space = "umap"),
  list(ref = "km_pc2",  oth = "km_pc198", ref_lab = "KMeans-2PC",
       oth_lab = "KM-PC198",  left_space = "pca",  right_space = "umap")
)

CLUSTER_COLORS <- cluster_colors

theme_base <- theme_minimal(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold", size = 11),
    plot.subtitle    = element_text(color = "grey50", size = 9),
    panel.grid.minor = element_blank()
  )

out_dir_comp_vor <- ensure_dir(DIR_FIG, "07_figures", "voronoi")

# ── 1. Hungarian alignment function ───────────────────────────────
align_hungarian <- function(ref, oth, k_val) {
  tab <- table(
    factor(ref, levels = as.character(seq_len(k_val))),
    factor(oth, levels = as.character(seq_len(k_val)))
  )
  asg <- solve_LSAP(tab, maximum = TRUE)
  map <- setNames(as.character(seq_len(k_val)), as.character(asg))
  as.character(map[as.character(oth)])
}

# ── 2. Voronoi tiles function ────────────────────────────────────────
voronoi_tiles <- function(x, y, labels, xlim, ylim) {
  cents <- data.frame(x = x, y = y, cl = as.character(labels)) %>%
    group_by(cl) %>%
    summarise(cx = mean(x), cy = mean(y), n = n(), .groups = "drop") %>%
    filter(n > 0)
  
  if (nrow(cents) < 2) return(NULL)
  
  # jitter if centroids coincide
  for (i in seq_len(nrow(cents))) {
    for (j in seq_len(nrow(cents))) {
      if (i >= j) next
      if (abs(cents$cx[i] - cents$cx[j]) < 1e-10 &&
          abs(cents$cy[i] - cents$cy[j]) < 1e-10) {
        cents$cx[i] <- cents$cx[i] + runif(1, 1e-8, 1e-7)
        cents$cy[i] <- cents$cy[i] + runif(1, 1e-8, 1e-7)
      }
    }
  }
  
  pad <- 0.15
  xr  <- xlim + c(-1, 1) * diff(xlim) * pad
  yr  <- ylim + c(-1, 1) * diff(ylim) * pad
  
  vor <- tryCatch(
    deldir(cents$cx, cents$cy,
           rw = c(xr[1], xr[2], yr[1], yr[2]),
           suppressMsge = TRUE),
    error = function(e) NULL
  )
  if (is.null(vor)) return(NULL)
  
  tiles <- tile.list(vor)
  rows  <- lapply(seq_along(tiles), function(i) {
    t <- tiles[[i]]
    if (length(t$x) == 0) return(NULL)
    tx    <- mean(t$x, na.rm = TRUE)
    ty    <- mean(t$y, na.rm = TRUE)
    dists <- sqrt((cents$cx - tx)^2 + (cents$cy - ty)^2)
    best  <- cents$cl[which.min(dists)]
    data.frame(x = t$x, y = t$y, cl = best, stringsAsFactors = FALSE)
  })
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0) return(NULL)
  do.call(rbind, rows)
}

# ── 3. Plot functions with Voronoi ─────────────────────────────────────

make_voronoi_scatter <- function(df_plot, x_col, y_col, cl_col,
                                 k_val, title, subtitle, xlab, ylab) {
  x      <- df_plot[[x_col]]
  y      <- df_plot[[y_col]]
  cl     <- as.character(df_plot[[cl_col]])
  disc   <- df_plot$discordant
  lev    <- as.character(seq_len(k_val))
  xlim   <- range(x, na.rm = TRUE)
  ylim   <- range(y, na.rm = TRUE)
  pal    <- CLUSTER_COLORS[seq_len(k_val)]
  
  vor_df <- voronoi_tiles(x, y, cl, xlim, ylim)
  
  df_concordant  <- df_plot[!disc, ]
  df_discordant  <- df_plot[disc,  ]
  cents_df <- data.frame(x = x, y = y, cl = cl) %>%
    group_by(cl) %>%
    summarise(cx = mean(x), cy = mean(y), .groups = "drop") %>%
    mutate(cl = factor(cl, levels = lev))
  
  p <- ggplot()
  
  if (!is.null(vor_df) && nrow(vor_df) > 0) {
    vor_df <- vor_df %>% mutate(cl = factor(cl, levels = lev))
    p <- p +
      geom_polygon(
        data = vor_df,
        aes(x = x, y = y, group = cl, fill = cl),
        alpha = 0.18, color = "grey60", linewidth = 0.4
      )
  }
  
  p +
    # concordant points coloured
    geom_point(
      data  = df_concordant,
      aes(x = .data[[x_col]], y = .data[[y_col]],
          color = factor(.data[[cl_col]], levels = lev)),
      alpha = 0.45, size = 0.9
    ) +
    # discordant points in black
    geom_point(
      data  = df_discordant,
      aes(x = .data[[x_col]], y = .data[[y_col]]),
      color = "black", alpha = 0.6, size = 1.2
    ) +
    # centroids
    geom_point(
      data  = cents_df,
      aes(x = cx, y = cy, color = cl),
      shape = 8, size = 4, stroke = 1.3
    ) +
    scale_color_manual(values = pal, name = "Cluster") +
    scale_fill_manual( values = pal, guide = "none") +
    coord_cartesian(xlim = xlim, ylim = ylim) +
    labs(title = title, subtitle = subtitle, x = xlab, y = ylab) +
    theme_base
}

make_alluvial <- function(df_plot, k_val, ref_lab, oth_lab) {
  lev <- as.character(seq_len(k_val))
  df_al <- as.data.frame(table(
    ref = factor(df_plot$cl_ref, levels = lev),
    oth = factor(df_plot$cl_oth, levels = lev)
  )) %>%
    rename(n = Freq) %>%
    filter(n > 0)
  
  ggplot(df_al, aes(axis1 = ref, axis2 = oth, y = n)) +
    geom_alluvium(aes(fill = ref), alpha = 0.75, width = 0.3) +
    geom_stratum(aes(fill = after_stat(stratum)),
                 width = 0.3, color = "white", linewidth = 0.4) +
    geom_text(stat = "stratum",
              aes(label = after_stat(stratum)),
              size = 3.5, fontface = "bold", color = "white") +
    scale_x_discrete(
      limits = c("ref", "oth"),
      labels = c(ref_lab, oth_lab),
      expand = c(0.15, 0.15)
    ) +
    scale_fill_manual(values = CLUSTER_COLORS[seq_len(k_val)]) +
    labs(
      title = paste0("Cluster flow: ", ref_lab, " <-> ", oth_lab,
                     "  (K = ", k_val, ")"),
      y = "Observations"
    ) +
    theme_base +
    theme(legend.position = "none", panel.grid = element_blank())
}

# ── 4. Main loop ───────────────────────────────────────────────
for (pair in method_pairs) {
  
  ref_prefix <- pair$ref
  oth_prefix <- pair$oth
  ref_lab    <- pair$ref_lab
  oth_lab    <- pair$oth_lab
  
  pair_dir <- file.path(out_dir_comp_vor,
                        paste0(ref_prefix, "_vs_", oth_prefix))
  if (!dir.exists(pair_dir)) dir.create(pair_dir, recursive = TRUE)
  
  for (k_val in k_values) {
    
    col_ref <- paste0(ref_prefix, "_k", k_val)
    col_oth <- paste0(oth_prefix, "_k", k_val)
    if (ref_prefix == "pam")     col_ref <- paste0("pam_", k_val)
    if (oth_prefix == "pam")     col_oth <- paste0("pam_", k_val)
    if (ref_prefix == "km_umap") col_ref <- paste0("km_umap_", k_val)
    if (oth_prefix == "km_umap") col_oth <- paste0("km_umap_", k_val)
    
    if (!col_ref %in% colnames(all_clusters) ||
        !col_oth %in% colnames(all_clusters)) next
    
    lbl_ref    <- as.character(all_clusters[[col_ref]])
    lbl_oth    <- as.character(all_clusters[[col_oth]])
    lbl_oth_al <- align_hungarian(lbl_ref, lbl_oth, k_val)
    pct_disc   <- round(mean(lbl_ref != lbl_oth_al) * 100, 1)
    
    df_plot <- data.frame(
      PC1        = df_pca$PC1,
      PC2        = df_pca$PC2,
      UMAP1      = df_umap$UMAP1,
      UMAP2      = df_umap$UMAP2,
      MDS1       = df_mds$MDS1,       # new
      MDS2       = df_mds$MDS2,       # new
      cl_ref     = lbl_ref,
      cl_oth     = lbl_oth_al,
      discordant = lbl_ref != lbl_oth_al
    )
    
    # axes and labels of the left panel
    left_space  <- pair$left_space
    right_space <- pair$right_space
    
    if (left_space == "mds") {
      left_x <- "MDS1"; left_y <- "MDS2"
      left_xlab <- "MDS 1"; left_ylab <- "MDS 2"
      left_title <- paste0("MDS-Gower space — ", oth_lab, " labels  (K = ", k_val, ")")
      left_cl    <- "cl_oth"   # PAM is oth: its labels go in its own space
    } else {
      left_x <- "PC1"; left_y <- "PC2"
      left_xlab <- "PC 1"; left_ylab <- "PC 2"
      left_title <- paste0("PCA space — ", ref_lab, " labels  (K = ", k_val, ")")
      left_cl    <- "cl_ref"
    }
    
    if (right_space == "pca") {
      right_x <- "PC1"; right_y <- "PC2"
      right_xlab <- "PC 1"; right_ylab <- "PC 2"
      right_title <- paste0("PCA space — ", ref_lab, " labels  (K = ", k_val, ")")
      right_cl    <- "cl_ref"
    } else {
      right_x <- "UMAP1"; right_y <- "UMAP2"
      right_xlab <- "UMAP 1"; right_ylab <- "UMAP 2"
      right_title <- paste0("UMAP space — ", oth_lab, " labels  (K = ", k_val, ")")
      right_cl    <- "cl_oth"
    }
    
    p_left <- make_voronoi_scatter(
      df_plot, left_x, left_y, left_cl, k_val,
      title    = left_title,
      subtitle = paste0("Black = discordant  [", pct_disc, "%]"),
      xlab = left_xlab, ylab = left_ylab
    )
    
    p_right <- make_voronoi_scatter(
      df_plot, right_x, right_y, right_cl, k_val,
      title    = right_title,
      subtitle = paste0("Black = discordant  [", pct_disc, "%]"),
      xlab = right_xlab, ylab = right_ylab
    )
    # left panel: PCA space, colour = ref label
    p_left <- make_voronoi_scatter(
      df_plot, "PC1", "PC2", "cl_ref", k_val,
      title    = paste0("PCA space — ", ref_lab, " labels  (K = ", k_val, ")"),
      subtitle = paste0("Black = discordant with ", oth_lab,
                        "  [", pct_disc, "%]"),
      xlab = "PC 1", ylab = "PC 2"
    )
    
    # central panel: alluvial
    p_mid <- make_alluvial(df_plot, k_val, ref_lab, oth_lab)
    
    # right panel: UMAP space, colour = oth label
    p_right <- make_voronoi_scatter(
      df_plot, "UMAP1", "UMAP2", "cl_oth", k_val,
      title    = paste0("UMAP space — ", oth_lab, " labels  (K = ", k_val, ")"),
      subtitle = paste0("Black = discordant with ", ref_lab,
                        "  [", pct_disc, "%]"),
      xlab = "UMAP 1", ylab = "UMAP 2"
    )
    
    panel <- (p_left | p_mid | p_right) +
      plot_annotation(
        title   = paste0(ref_lab, "  vs  ", oth_lab,
                         "  —  K = ", k_val, "  (Voronoi)"),
        caption = paste0("Reference: ", ref_lab,
                         "  |  Hungarian alignment  |  ",
                         pct_disc, "% discordant  |  * = centroid"),
        theme   = theme(
          plot.title   = element_text(face = "bold", size = 13, hjust = 0.5),
          plot.caption = element_text(color = "grey50", size = 8)
        )
      )
    
    fname <- sprintf("voronoi_%s_vs_%s_k%d.pdf",
                     ref_prefix, oth_prefix, k_val)
    ggsave(file.path(pair_dir, fname),
           panel, width = 16, height = 5.5, bg = "white")
    
    cat(sprintf("Voronoi saved: %s vs %s  k=%d  disc=%.1f%%\n",
                ref_lab, oth_lab, k_val, pct_disc))
  }
}

cat("\nDone — Script 2 (Voronoi)\n")

# pam check ----
# ── EXTRA BLOCK: PAM vs KM-UMAP, PAM vs KM-PC2, PAM vs KM-PC5 ──────────

extra_pairs <- list(
  
  # ── PAM (ref) vs KM-UMAP (oth) — Voronoi in UMAP space ─────────────────────
  # left panel: UMAP Voronoi, colour = PAM, black = discordant with KM-UMAP
  # right panel: alluvial PAM <-> KM-UMAP
  list(
    ref        = "pam",
    oth        = "km_umap",
    ref_lab    = "PAM",
    oth_lab    = "KM-UMAP",
    vor_space  = "umap",   # space in which the Voronoi is drawn
    vor_cl     = "ref"     # colour = ref label (PAM)
  ),
  
  # ── PAM (ref) vs KM-PC2 (oth) — Voronoi in PC space ───────────────────────
  # left panel: PCA Voronoi (PC1/PC2), colour = PAM, black = discordant with KM-PC2
  # right panel: alluvial PAM <-> KM-PC2
  list(
    ref        = "pam",
    oth        = "km_pc2",
    ref_lab    = "PAM",
    oth_lab    = "KM-PC2",
    vor_space  = "pca",
    vor_cl     = "ref"
  ),
  
  # ── PAM (ref) vs KM-PC5 (oth) — Voronoi in PC space, colour = KM-PC5 ──────
  # left panel: PCA Voronoi (PC1/PC2), colour = KM-PC5, black = discordant with PAM
  # right panel: alluvial PAM <-> KM-PC5
  list(
    ref        = "pam",
    oth        = "km_pc5",
    ref_lab    = "PAM",
    oth_lab    = "KM-PC5",
    vor_space  = "pca",
    vor_cl     = "oth"    # colour = oth label (KM-PC5)
  )
)

out_dir_extra <- file.path(out_dir_comp_vor, "extra_pam")
if (!dir.exists(out_dir_extra)) dir.create(out_dir_extra, recursive = TRUE)

for (pair in extra_pairs) {
  
  ref_prefix <- pair$ref
  oth_prefix <- pair$oth
  ref_lab    <- pair$ref_lab
  oth_lab    <- pair$oth_lab
  vor_space  <- pair$vor_space   # "umap" or "pca"
  vor_cl     <- pair$vor_cl      # "ref" or "oth"
  
  pair_dir <- file.path(out_dir_extra,
                        paste0(ref_prefix, "_vs_", oth_prefix))
  if (!dir.exists(pair_dir)) dir.create(pair_dir, recursive = TRUE)
  
  for (k_val in k_values) {
    
    # ── column names in all_clusters ────────────────────────────────
    col_ref <- paste0("pam_", k_val)                        # always PAM = ref
    col_oth <- if (oth_prefix == "km_umap") {
      paste0("km_umap_", k_val)
    } else {
      paste0(oth_prefix, "_k", k_val)                       # km_pc2, km_pc5 …
    }
    
    if (!col_ref %in% colnames(all_clusters) ||
        !col_oth %in% colnames(all_clusters)) next
    
    # ── labels and alignment ────────────────────────────────────
    lbl_ref    <- as.character(all_clusters[[col_ref]])
    lbl_oth    <- as.character(all_clusters[[col_oth]])
    lbl_oth_al <- align_hungarian(lbl_ref, lbl_oth, k_val)
    pct_disc   <- round(mean(lbl_ref != lbl_oth_al) * 100, 1)
    
    df_plot <- data.frame(
      PC1        = df_pca$PC1,
      PC2        = df_pca$PC2,
      UMAP1      = df_umap$UMAP1,
      UMAP2      = df_umap$UMAP2,
      MDS1       = df_mds$MDS1,
      MDS2       = df_mds$MDS2,
      cl_ref     = lbl_ref,
      cl_oth     = lbl_oth_al,
      discordant = lbl_ref != lbl_oth_al
    )
    
    # ── axes of the Voronoi panel ───────────────────────────────────
    if (vor_space == "umap") {
      vx <- "UMAP1"; vy <- "UMAP2"
      vxlab <- "UMAP 1"; vylab <- "UMAP 2"
      space_label <- "UMAP space"
    } else {
      vx <- "PC1"; vy <- "PC2"
      vxlab <- "PC 1"; vylab <- "PC 2"
      space_label <- "PCA space"
    }
    
    # colour of the concordant points
    cl_col     <- if (vor_cl == "ref") "cl_ref" else "cl_oth"
    color_lab  <- if (vor_cl == "ref") ref_lab  else oth_lab
    discord_with <- if (vor_cl == "ref") oth_lab else ref_lab
    
    vor_title <- paste0(space_label, " — ", color_lab,
                        " labels  (K = ", k_val, ")")
    vor_sub   <- paste0("Black = discordant with ", discord_with,
                        "  [", pct_disc, "%]")
    
    p_vor <- make_voronoi_scatter(
      df_plot, vx, vy, cl_col, k_val,
      title    = vor_title,
      subtitle = vor_sub,
      xlab     = vxlab,
      ylab     = vylab
    )
    
    # ── alluvial ────────────────────────────────────────────────────
    p_al <- make_alluvial(df_plot, k_val, ref_lab, oth_lab)
    
    # ── 2-panel layout ─────────────────────────────────────
    panel <- (p_vor | p_al) +
      plot_annotation(
        title   = paste0(ref_lab, "  vs  ", oth_lab,
                         "  —  K = ", k_val, "  (Voronoi)"),
        caption = paste0("Colour: ", color_lab,
                         "  |  Hungarian alignment  |  ",
                         pct_disc, "% discordant  |  * = centroid"),
        theme = theme(
          plot.title   = element_text(face = "bold", size = 13, hjust = 0.5),
          plot.caption = element_text(color = "grey50", size = 8)
        )
      )
    
    fname <- sprintf("voronoi_%s_vs_%s_k%d.pdf",
                     ref_prefix, oth_prefix, k_val)
    ggsave(file.path(pair_dir, fname),
           panel, width = 12, height = 5.5, bg = "white")
    
    cat(sprintf("Extra Voronoi saved: %s vs %s  k=%d  disc=%.1f%%\n",
                ref_lab, oth_lab, k_val, pct_disc))
  }
}

cat("\nDone — Blocco extra PAM comparisons\n")

# pam vs all ----
# ── BLOCK PAM vs KM-UMAP and PAM vs KM-PC5: fixed Voronoi cells ───────────────
# cells and centroids = always from the KM method (KM-UMAP or KM-PC5)
# point colour      = always PAM
# black             = discordant PAM vs KM
# ── Voronoi function with cells fixed on KM, colour = PAM ─────────────────────
make_voronoi_fixed_cells <- function(df_plot, x_col, y_col,
                                     cell_cl,          # cells and centroids: KM
                                     color_cl,         # point colour: PAM
                                     black_discordant, # TRUE = panel with black points
                                     k_val, title, subtitle, xlab, ylab) {
  x        <- df_plot[[x_col]]
  y        <- df_plot[[y_col]]
  cl_cell  <- as.character(df_plot[[cell_cl]])
  lev      <- as.character(seq_len(k_val))
  xlim     <- range(x, na.rm = TRUE)
  ylim     <- range(y, na.rm = TRUE)
  pal      <- CLUSTER_COLORS[seq_len(k_val)]
  disc     <- df_plot$discordant
  
  # Voronoi cells and centroids built on cell_cl (KM)
  vor_df <- voronoi_tiles(x, y, cl_cell, xlim, ylim)
  
  cents_df <- data.frame(x = x, y = y, cl = cl_cell) %>%
    group_by(cl) %>%
    summarise(cx = mean(x), cy = mean(y), .groups = "drop") %>%
    mutate(cl = factor(cl, levels = lev))
  
  p <- ggplot()
  
  if (!is.null(vor_df) && nrow(vor_df) > 0) {
    vor_df <- vor_df %>% mutate(cl = factor(cl, levels = lev))
    p <- p + geom_polygon(
      data  = vor_df,
      aes(x = x, y = y, group = cl, fill = cl),
      alpha = 0.18, color = "grey60", linewidth = 0.4
    )
  }
  
  if (black_discordant) {
    df_conc <- df_plot[!disc, ]
    df_disc <- df_plot[ disc, ]
    p <- p +
      geom_point(
        data  = df_conc,
        aes(x = .data[[x_col]], y = .data[[y_col]],
            color = factor(.data[[color_cl]], levels = lev)),
        alpha = 0.45, size = 0.9
      ) +
      geom_point(
        data  = df_disc,
        aes(x = .data[[x_col]], y = .data[[y_col]]),
        color = "black", alpha = 0.6, size = 1.2
      )
  } else {
    p <- p +
      geom_point(
        data = df_plot,
        aes(x = .data[[x_col]], y = .data[[y_col]],
            color = factor(.data[[color_cl]], levels = lev)),
        alpha = 0.45, size = 0.9
      )
  }
  
  p +
    geom_point(
      data  = cents_df,
      aes(x = cx, y = cy, color = cl),
      shape = 8, size = 4, stroke = 1.3
    ) +
    scale_color_manual(values = pal, name = "Cluster") +
    scale_fill_manual(values  = pal, guide = "none") +
    coord_cartesian(xlim = xlim, ylim = ylim) +
    labs(title = title, subtitle = subtitle, x = xlab, y = ylab) +
    theme_base
}
fixed_pairs <- list(
  
  # PAM vs KM-UMAP — UMAP space
  list(
    ref       = "pam",
    oth       = "km_umap",
    ref_lab   = "PAM",
    oth_lab   = "KM-UMAP",
    x_col     = "UMAP1",  y_col = "UMAP2",
    xlab      = "UMAP 1", ylab  = "UMAP 2",
    space_lab = "UMAP space",
    out_sub   = "pam_vs_kmumap_fixed"
  ),
  
  # PAM vs KM-PC5 — PCA space (PC1/PC2)
  list(
    ref       = "pam",
    oth       = "km_pc5",
    ref_lab   = "PAM",
    oth_lab   = "KM-PC5",
    x_col     = "PC1",  y_col = "PC2",
    xlab      = "PC 1", ylab  = "PC 2",
    space_lab = "PCA space",
    out_sub   = "pam_vs_kmpc5_fixed"
  ),
  # PAM vs KM-PC2 — PCA space (PC1/PC2)
  list(
    ref       = "pam",
    oth       = "km_pc2",
    ref_lab   = "PAM",
    oth_lab   = "KM-PC2",
    x_col     = "PC1",  y_col = "PC2",
    xlab      = "PC 1", ylab  = "PC 2",
    space_lab = "PCA space",
    out_sub   = "pam_vs_kmpc2_fixed")
)
for (fp in fixed_pairs) {
  
  ref_prefix <- fp$ref
  oth_prefix <- fp$oth
  ref_lab    <- fp$ref_lab
  oth_lab    <- fp$oth_lab
  x_col      <- fp$x_col
  y_col      <- fp$y_col
  xlab       <- fp$xlab
  ylab       <- fp$ylab
  space_lab  <- fp$space_lab
  
  out_dir_fp <- file.path(out_dir_comp_vor, fp$out_sub)
  if (!dir.exists(out_dir_fp)) dir.create(out_dir_fp, recursive = TRUE)
  
  # legend adapted to the pair
  legend_text_fp <- paste0(
    "Voronoi cells & centroids (*): computed from ", oth_lab,
    " cluster assignments in ", space_lab, ".  |  ",
    "Point colour: ", ref_lab, " cluster label.  |  ",
    "Black points (right panel): ", ref_lab, " and ", oth_lab, " disagree.  |  ",
    "Discordance after Hungarian alignment."
  )
  
  make_legend_fp <- function(txt) {
    ggplot() +
      annotate("text", x = 0.5, y = 0.5,
               label     = txt,
               hjust     = 0.5, vjust = 0.5,
               size      = 3.0, lineheight = 1.4) +
      theme_void() +
      theme(plot.margin = unit(c(2, 2, 2, 2), "pt"))
  }
  
  for (k_val in k_values) {
    
    # ── column names ──────────────────────────────────────────────
    col_ref <- paste0("pam_", k_val)
    col_oth <- if (oth_prefix == "km_umap") {
      paste0("km_umap_", k_val)
    } else {
      paste0(oth_prefix, "_k", k_val)
    }
    
    if (!col_ref %in% colnames(all_clusters) ||
        !col_oth %in% colnames(all_clusters)) next
    
    # ── labels and alignment ──────────────────────────────────
    lbl_ref    <- as.character(all_clusters[[col_ref]])
    lbl_oth    <- as.character(all_clusters[[col_oth]])
    lbl_oth_al <- align_hungarian(lbl_ref, lbl_oth, k_val)
    pct_disc   <- round(mean(lbl_ref != lbl_oth_al) * 100, 1)
    
    df_plot <- data.frame(
      PC1        = df_pca$PC1,
      PC2        = df_pca$PC2,
      UMAP1      = df_umap$UMAP1,
      UMAP2      = df_umap$UMAP2,
      cl_ref     = lbl_ref,
      cl_oth     = lbl_oth_al,
      discordant = lbl_ref != lbl_oth_al
    )
    
    # ── panel 1: all coloured by PAM, cells from KM ───────────
    p_allcolor <- make_voronoi_fixed_cells(
      df_plot,
      x_col            = x_col,    y_col = y_col,
      cell_cl          = "cl_oth", # cells from KM
      color_cl         = "cl_ref", # colour from PAM
      black_discordant = FALSE,
      k_val    = k_val,
      title    = paste0(space_lab, " — ", ref_lab, " labels  (K = ", k_val, ")"),
      subtitle = paste0("All points coloured by ", ref_lab,
                        "  |  Voronoi cells from ", oth_lab),
      xlab = xlab, ylab = ylab
    )
    
    # ── panel 2: alluvial ──────────────────────────────────────
    p_al <- make_alluvial(df_plot, k_val, ref_lab, oth_lab)
    
    # ── panel 3: discordant in black ──────────────────────────────
    p_black <- make_voronoi_fixed_cells(
      df_plot,
      x_col            = x_col,    y_col = y_col,
      cell_cl          = "cl_oth", # cells from KM
      color_cl         = "cl_ref", # colour from PAM
      black_discordant = TRUE,
      k_val    = k_val,
      title    = paste0(space_lab, " — ", ref_lab, " labels  (K = ", k_val, ")"),
      subtitle = paste0("Black = discordant with ", oth_lab,
                        "  [", pct_disc, "%]  |  Voronoi cells from ", oth_lab),
      xlab = xlab, ylab = ylab
    )
    
    # ── horizontal legend ───────────────────────────────────────
    p_leg <- make_legend_fp(legend_text_fp)
    
    # ── layout ──────────────────────────────────────────────
    panel <- (p_allcolor | p_al | p_black) /
      p_leg +
      plot_layout(heights = c(6, 1)) +
      plot_annotation(
        title   = paste0(ref_lab, "  vs  ", oth_lab,
                         "  —  K = ", k_val,
                         "  (Voronoi cells from ", oth_lab, ")"),
        caption = paste0("Hungarian alignment  |  ",
                         pct_disc, "% discordant  |  * = ",
                         oth_lab, " centroid"),
        theme = theme(
          plot.title   = element_text(face = "bold", size = 13, hjust = 0.5),
          plot.caption = element_text(color = "grey50", size = 8)
        )
      )
    
    fname <- sprintf("voronoi_%s_vs_%s_fixed_k%d.pdf",
                     ref_prefix, oth_prefix, k_val)
    ggsave(file.path(out_dir_fp, fname),
           panel, width = 17, height = 7, bg = "white")
    
    cat(sprintf("Fixed-cell Voronoi saved: %s vs %s  k=%d  disc=%.1f%%\n",
                ref_lab, oth_lab, k_val, pct_disc))
  }
}

cat("\nDone — Fixed-cell Voronoi (PAM vs KM-UMAP + PAM vs KM-PC5)\n")
