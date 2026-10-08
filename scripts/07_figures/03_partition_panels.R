# =============================================================================
# scripts/07_figures/03_partition_panels.R
# Purpose : Voronoi panels of the aligned partitions (KM-UMAP, PAM in MDS and
#           PCA space, KM-PC5, KM-PC2) with a short description of each
#           cluster written in its cell, and the four-method comparison in
#           PCA space (PAM cells built on its medoids).
# Thesis  : Figure 4.7
# Input   : results/objects/clusters/{all_clusters_aligned,umap_df_clusters}.rds
#           results/objects/{pca_scores,pam_mds_df_clusters,pam_results_allk}.rds
# Output  : results/figures/07_figures/partition_panels/**/*.pdf
# =============================================================================

source(here::here("R", "00_setup.R"))


# helper: geometric centroid of a polygon (shoelace formula)
poly_centroid <- function(x, y) {
  n <- length(x)
  x2 <- c(x[-1], x[1]); y2 <- c(y[-1], y[1])
  cross <- x * y2 - x2 * y
  A <- sum(cross) / 2
  if (abs(A) < 1e-9) return(c(mean(x), mean(y)))  # degenerate case
  cx <- sum((x + x2) * cross) / (6 * A)
  cy <- sum((y + y2) * cross) / (6 * A)
  c(cx, cy)
}

plot_voronoi_aligned <- function(coord_df, aligned_df, k,
                                 cl_prefix, xcol, ycol, title,
                                 cell_text = NULL, fixed = TRUE,medoid_ids = NULL) {
  
  cl_col <- paste0(cl_prefix, k)
  stopifnot(cl_col %in% names(aligned_df))
  stopifnot(nrow(coord_df) == nrow(aligned_df))
  
  df <- data.frame(
    X   = coord_df[[xcol]],
    Y   = coord_df[[ycol]],
    .cl = factor(aligned_df[[cl_col]])
  )
  
  if (is.null(medoid_ids)) {
    centroids <- df %>%
      group_by(.cl) %>%
      summarise(X = mean(X), Y = mean(Y), .groups = "drop")
  } else {
    # actual medoids: use the medoid row, with its ALIGNED label
    centroids <- data.frame(
      .cl = df$.cl[medoid_ids],
      X   = df$X[medoid_ids],
      Y   = df$Y[medoid_ids]
    )
    centroids <- centroids[order(centroids$.cl), ]
  }
  print(centroids) 
  pad  <- 0.1
  xlim <- range(df$X) + c(-pad, pad) * diff(range(df$X))
  ylim <- range(df$Y) + c(-pad, pad) * diff(range(df$Y))
  
  vor   <- deldir(centroids$X, centroids$Y, rw = c(xlim, ylim))
  tiles <- tile.list(vor)
  
  vor_df <- do.call(rbind, lapply(seq_along(tiles), function(i) {
    data.frame(x = tiles[[i]]$x, y = tiles[[i]]$y,
               cluster = centroids$.cl[i])
  }))
  vor_df$cluster <- factor(vor_df$cluster, levels = levels(centroids$.cl))
  
  p <- ggplot() +
    geom_polygon(
      data = vor_df,
      aes(x = x, y = y, fill = cluster, group = cluster),
      alpha = 0.45, color = "grey40", linewidth = 0.35
    ) +
    scale_fill_manual(values = cluster_colors, name = "Cluster") +
    labs(title = sprintf("%s - K = %d", title, k), x = xcol, y = ycol) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position  = if (k <= 4) "right" else "bottom",
      legend.key.size  = unit(0.4, "cm"),
      panel.grid.minor = element_blank()
    ) +
    guides(fill = guide_legend(
      nrow = if (k <= 8) 2 else 4,
      override.aes = list(alpha = 1)
    ))
  
  if (fixed) p <- p + coord_fixed(xlim = xlim, ylim = ylim, expand = FALSE)
  
  # ── text: anchored to the centroid of the CELL, with repel ──
  if (!is.null(cell_text)) {
    lab_pos <- do.call(rbind, lapply(split(vor_df, vor_df$cluster), function(d) {
      cc <- poly_centroid(d$x, d$y)
      data.frame(.cl = d$cluster[1], X = cc[1], Y = cc[2])
    }))
    lab_pos$label <- cell_text[as.character(lab_pos$.cl)]
    lab_pos$label[is.na(lab_pos$label)] <- ""
    
    p <- p + ggrepel::geom_text_repel(
      data        = lab_pos,
      aes(x = X, y = Y, label = label),
      size        = if (k <= 4) 3.0 else 2.5,
      fontface    = "bold",
      lineheight  = 0.95,
      color       = "black",
      box.padding = 0.4,
      point.size  = NA,
      min.segment.length = Inf,
      seed        = 1
    )
  }
  
  p
}

# needed objects 



all_cluster_aligned <- readRDS(PATH_CLUSTERS)
pca_scores  <- readRDS(file.path(DIR_OBJ, "pca_scores.rds"))
mds_df      <- readRDS(file.path(DIR_OBJ, "pam_mds_df_clusters.rds"))
results_pam <- readRDS(file.path(DIR_OBJ, "pam_results_allk.rds"))
umap_df     <- readRDS(file.path(DIR_CLUST, "umap_df_clusters.rds"))

scores_mat <- if (inherits(pca_scores, "prcomp")) pca_scores$x else as.matrix(pca_scores)

n_center<- k_values <- c(2,3, 4,5,6)
pca_df <- data.frame( PC1 = pca_scores[,1],
                      PC2 = pca_scores[,2])

# TEXTS ####

# PCA 2 
texts_pc2_k2 <- c(
  "1" = "RELIGIOUS CONSERVATIVES\nAnti-immigration\nDistrust intl. institutions\nReligion > Science",
  "2" = "SECULAR PROGRESSIVES\nPro-immigration, pro-science\nTrust intl. institutions"
)
texts_pc2_k3 <- c(
  "1" = "ECONOMICALLY VULNERABLE\nInstitutional distrust\nUnsafe, low trust",
  "2" = "SECULAR PROGRESSIVES\nPro-immigration\nEconomically stable",
  "3" = "DEVOUT CONSERVATIVES\nReligious, pro-authority\nRelatively safe"
)
texts_pc2_k4 <- c(
  "1" = "ECONOMICALLY VULNERABLE\nInstitutional distrust\nLow social trust",
  "2" = "RELIGIOUS CONSERVATIVES\nAnti-immigration\nTrust in government",
  "3" = "SECULAR PROGRESSIVES\nPro-environment\nNon-religious",
  "4" = "RELIGIOUS MODERATES\nBelieve God/Hell\nModerate, pro-authority"
)
texts_pc2_k5 <- c(
  "1" = "CENTRAL / MODERATE\nWeakly differentiated",
  "2" = "INSTITUTIONAL TRUSTERS\nSafe, faith not fanatic",
  "3" = "RELIGIOUS CONSERVATIVES\nAnti-immigration",
  "4" = "ECONOMICALLY VULNERABLE\nDisorder, high distrust",
  "5" = "SECULAR PROGRESSIVES\nStrongly pro-immigration"
)
texts_pc2_k6 <- c(
  "1" = "PROGRESSIVE MODERATES\nPro-science, egalitarian",
  "2" = "SECULAR PROGRESSIVES\nPro-immigration, wealthy",
  "3" = "RELIGIOUS CONSERVATIVES\nAnti-immigration",
  "4" = "VULNERABLE / DISORDER\nDistrust",
  "5" = "SECULAR PROGRESSIVES\nLow wellbeing, anti-intl.",
  "6" = "CENTRAL / MODERATE\nWeakly differentiated"
)
cell_texts_pc2 <- list(k2=texts_pc2_k2, k3=texts_pc2_k3, k4=texts_pc2_k4,
                       k5=texts_pc2_k5, k6=texts_pc2_k6)

#PCA5
texts_pc5_k2 <- c(
  "1" = "RELIGIOUS CONSERVATIVES\nAnti-immigration\nDistrust intl. institutions\nReligion > Science",
  "2" = "SECULAR PROGRESSIVES\nPro-immigration, pro-science\nTrust intl. institutions"
)
texts_pc5_k3 <- c(
  "1" = "ECONOMICALLY VULNERABLE\nInstitutional distrust\nUnsafe, low trust",
  "2" = "SECULAR PROGRESSIVES\nPro-immigration\nEconomically stable",
  "3" = "DEVOUT CONSERVATIVES\nReligious, pro-authority\nInstitutional trust, safe"
)
texts_pc5_k4 <- c(
  "1" = "ECONOMICALLY VULNERABLE\nInstitutional distrust\nLow social trust",
  "2" = "TOLERANT / CIVIC\nEconomically secure\nEgalitarian, pro-immigrant",
  "3" = "SECULAR PROGRESSIVES\nPro-immigration, wealthy\nNon-religious",
  "4" = "RELIGIOUS MODERATES\nConservative, pro-authority\nInstitutional trust, safe"
)
texts_pc5_k5 <- c(
  "1" = "CENTRAL / UNSTABLE\nNon-egalitarian gender\nReligion > Science",
  "2" = "MODERATE RELIGIOSITY\nInstitutional trust\nPro-authority",
  "3" = "ACTIVE RELIGIOSITY\nConservative, anti-immigration\nSelective anti-intl.",
  "4" = "ECONOMICALLY VULNERABLE\nDisorder, high distrust",
  "5" = "SECULAR / ANTI-ESTABLISHMENT\nPro-environment\nAnti-authority"
)
texts_pc5_k6 <- c(
  "1" = "MODERATE RELIGIOSITY\nInstitutional trust\nMixed on immigration",
  "2" = "SECULAR PROGRESSIVES\nPro-immigration\nNon-religious",
  "3" = "ACTIVE RELIGIOSITY\nSelective anti-intl. (IMF)\nTrust national gov.",
  "4" = "ECONOMICALLY VULNERABLE\nDisorder, high distrust",
  "5" = "VULNERABLE / UNSTABLE\nNorm-violation contrast",
  "6" = "SECULAR NATIONALISTS\nNon-religious\nAnti-immigration"
)
cell_texts_pc5 <- list(k2=texts_pc5_k2, k3=texts_pc5_k3, k4=texts_pc5_k4,
                       k5=texts_pc5_k5, k6=texts_pc5_k6)


# UMAP
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
  "1" = "Economic vulnerability\nlow religious practice\nUnsafe neighbors\npro parents beating their childre(a lot)\nReligion > Science",
  "2" = "Strong religiosity\nInternational organization trust\nPro science attitude\n(wealthy)\nReligion < Science",
  "3" = "Liberals\nPro-immigrant\nStrongly non religious\nReligion < Science",
  "4" = "No trust in international organizations\nActive Religiosity\nConservatives\n(wealthy)\nReligion > Science"
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
  "1" = "Religious\nTrust in international organization (pro NATO)\nLiberal\nsafe neighbors\n/nReligion > Science  ",
  "2" = "Non Religious\nNo respect for authority\nPro immigrations\nProscience\nLiberals/nReligion < Science",
  "3" = "Religious/nConfidence in National institutions/nNo trust in the press/nEconomic > Enviroment\nWealthy",
  "4" = "Strong Religiosity/nNo trust in Economical system/nAginst NATO (strongly)/nEconomic Instability/nUNsafe nieghbors/nReligion < Science",
  "5" = "Unsafe neighbors/nEconomic instability/nTraditional values/nReligion > Science"  ,
  "6" = "Non religious\nAginst immigrants\n No pro science"
)

texts_k6 <- c(
  "1" = "RELIGIOUS LIBERALS\nSafe neighborhood\nTrust in international organizations (pro-NATO)\nReligion < Science",
  "2" = "PROGRESSIVES\nNon-religious\nPro-immigration\nPro-science\nNo respect for authority\nReligion < Science",
  "3" = "RELIGIOUS INSTITUTIONALISTS\nReligious\nConfidence in national institutions\nDistrust in the press\nEconomy > Environment\n Economically stable",
  "4" = "STRONG RELIGIOSITY\nDistrust in economic system\nAnti-NATO (strongly)",
  "5" = "ECONOMICALLY VULNERABLE TRADITIONALIST\nTraditional values\nUnsafe neighborhood\nReligion > Science",
  "6" = "SECULAR TRADITIONALIST\nNon-religious\nAnti-immigration\nNot pro-science"
)

cell_texts_umap <- list(
  k2 = texts_k2,
  k3 = texts_k3,
  k4 = texts_k4,
  k5 = texts_k5,
  k6 = texts_k6
)

# PAM 
texts_pam_k2 <- c(
  "1" = "ACTIVE RELIGIOSITY\nConservative, anti-immigration\nDistrust intl. institutions",
  "2" = "SECULAR / TRUSTING\nPro-immigration, safe\nTrust intl. institutions"
)
texts_pam_k3 <- c(
  "1" = "LOW SOCIAL TRUST\nSecurity-focused\nDistrust government",
  "2" = "SECULAR / PRO-SCIENCE\nEconomically stable, safe",
  "3" = "ACTIVE RELIGIOSITY\nConservative, nationalist"
)
texts_pam_k4 <- c(
  "1" = "INSTITUTIONAL DISTRUST\nSecurity concerns\nWary of others",
  "2" = "SECULAR MODERATES\nCivically engaged\nPro-immigration",
  "3" = "SECULAR PROGRESSIVES\nPro-immigration\nDigital privacy",
  "4" = "ACTIVE RELIGIOSITY\nConservative, nationalist"
)
texts_pam_k5 <- c(
  "1" = "SECULAR / INSECURE\nLiberal, unsafe area",
  "2" = "SECULAR / TOLERANT\nCivically engaged, happy",
  "3" = "ACTIVE RELIGIOSITY\nConservative, nationalist",
  "4" = "MODERATE RELIGIOSITY\nDistrust, security concerns",
  "5" = "SECULAR / ANTI-INTL\nNeighborhood disorder"
)
texts_pam_k6 <- c(
  "1" = "SECULAR MODERATES\nSocially trusting",
  "2" = "SECULAR PROGRESSIVES\nPro-immigration",
  "3" = "ACTIVE RELIGIOSITY\nConservative, nationalist",
  "4" = "MODERATE RELIGIOSITY\nDistrust, insecurity",
  "5" = "ANTI-INTERNATIONALISM\nAnti-immigration, insecure",
  "6" = "SECULAR-LIBERAL\nPersonally insecure"
)
cell_texts_pam <- list(k2=texts_pam_k2, k3=texts_pam_k3, k4=texts_pam_k4,
                       k5=texts_pam_k5, k6=texts_pam_k6)

# ── UMAP (via the generic function) ──####


plots_umap <- lapply(k_values, function(k)
  plot_voronoi_aligned(umap_df, all_cluster_aligned, k,
                       cl_prefix = "km_umap_",
                       xcol = "UMAP1", ycol = "UMAP2",
                       title = "Voronoi (UMAP)",
                      fixed = FALSE,
                      cell_text = cell_texts_umap[[paste0("k", k)]]))
names(plots_umap) <- paste0("k", k_values)
# ── PAM - #####

plots_pam_mds <- lapply(n_centers, function(k)
  plot_voronoi_aligned(mds_df, all_cluster_aligned, k,
                       cl_prefix = "pam_",
                       xcol = "MDS1", ycol = "MDS2",
                       title = "PAM (Gower) - MDS",
                       fixed = FALSE,
                       cell_text = cell_texts_pam[[paste0("k", k)]],
                       medoid_ids = results_pam[[as.character(k)]]$id.med))   # MDS: axes on different scales
names(plots_pam_mds) <- paste0("k", n_centers)

# ── PAM - pca ####

plots_pam_pca <- lapply(n_centers, function(k)
  plot_voronoi_aligned(pca_df, all_cluster_aligned, k,
                       cl_prefix = "pam_",
                       xcol = "PC1", ycol = "PC2",
                       title = "PAM (Gower) - PCA",
                       fixed = FALSE,
                       cell_text = cell_texts_pam[[paste0("k", k)]],
                       medoid_ids = results_pam[[as.character(k)]]$id.med))
names(plots_pam_pca) <- paste0("k", n_centers)

# ── K-means on 5 PCs ──####
plots_kpca5 <- lapply(k_values, function(k)
  plot_voronoi_aligned(pca_scores, all_cluster_aligned, k,
                       cl_prefix = "km_pc5_k",   
                       xcol = "PC1", ycol = "PC2",
                       title = "K-means PCA (5 PC)",
                       fixed = FALSE,
                       cell_texts_pc5[[paste0("k", k)]]))
names(plots_kpca5) <- paste0("k", k_values)

# ── K-means on 2 PCs ──####
 
  plots_kpca2 <- lapply(k_values, function(k)
    plot_voronoi_aligned(pca_scores, all_cluster_aligned, k,
                         cl_prefix = "km_pc2_k",   # <-- check against the column names
                         xcol = "PC1", ycol = "PC2",
                         title = "K-means PCA (2 PC)",
                         fixed = FALSE,
                         cell_text = cell_texts_pc2[[paste0("k", k)]]))
  names(plots_kpca2) <- paste0("k", k_values)
  
  
  
# SAVINGS ####
  # ── Output folders by method ────────────────────────────────────
  base_dir <- ensure_dir(DIR_FIG, "07_figures", "partition_panels")
  
  out_dirs <- list(
    umap = file.path(base_dir, "UMAP"),
    pc2  = file.path(base_dir, "PCA", "2PC"),
    pc5  = file.path(base_dir, "PCA", "5PC"),
    pam.mds  = file.path(base_dir, "PAM", "mds"),
    pam.pca  = file.path(base_dir, "PAM", "pca")
  )
  
  for (d in out_dirs) if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  
  # ── Helper: save a list of plots (named k2, k3, ...) ─────────────
  save_vor_plots <- function(plot_list, out_dir, width = 8, height = 6.5) {
    for (nm in names(plot_list)) {          # nm = "k2", "k3", ...
      ggsave(
        filename = file.path(out_dir, paste0(nm, ".pdf")),
        plot     = plot_list[[nm]],
        width    = width, height = height
      )
    }
    cat("Saved", length(plot_list), "grafici in", out_dir, "\n")
  }
  
# ── Save ──────────────────────────────────────────────────────
  save_vor_plots(plots_umap,  out_dirs$umap)
  save_vor_plots(plots_kpca2, out_dirs$pc2)
  save_vor_plots(plots_kpca5, out_dirs$pc5)   # texts for PC5
  save_vor_plots(plots_pam_mds,   out_dirs$pam.mds)
  save_vor_plots(plots_pam_pca,   out_dirs$pam.pca)

# UMAP in PCA ####
  # ── Seed in UMAP space (centroid or medoid), one per cluster ───────
  umap_seed_ids <- function(umap_df, aligned_df, k,
                            cl_prefix = "km_umap_",
                            method = c("nearest_mean", "medoid")) {
    method <- match.arg(method)
    cl  <- factor(aligned_df[[paste0(cl_prefix, k)]])
    U1  <- umap_df$UMAP1; U2 <- umap_df$UMAP2
    ids <- integer(0)
    for (lv in levels(cl)) {
      idx <- which(cl == lv)
      if (method == "nearest_mean") {
        cx <- mean(U1[idx]); cy <- mean(U2[idx])
        d2 <- (U1[idx] - cx)^2 + (U2[idx] - cy)^2
        ids <- c(ids, idx[which.min(d2)])
      } else {                                  # UMAP medoid (robust)
        D <- as.matrix(dist(cbind(U1[idx], U2[idx])))
        ids <- c(ids, idx[which.min(rowSums(D))])
      }
    }
    ids
  }
  
  # ── How many observations fall OUTSIDE their own cell (in PCA) ────
  voronoi_outside_count <- function(coord_df, aligned_df, k, seed_ids,
                                    cl_prefix, xcol, ycol) {
    cl <- factor(aligned_df[[paste0(cl_prefix, k)]])
    X  <- coord_df[[xcol]]; Y <- coord_df[[ycol]]
    
    seeds <- data.frame(.cl = cl[seed_ids], sx = X[seed_ids], sy = Y[seed_ids])
    seeds <- seeds[order(seeds$.cl), ]
    
    D <- outer(X, seeds$sx, `-`)^2 + outer(Y, seeds$sy, `-`)^2   # squared distance to the seeds
    nearest <- seeds$.cl[max.col(-D, ties.method = "first")]     # assigned cell
    
    outside <- tapply(nearest != cl, cl, sum)
    ntot    <- tapply(rep(1L, length(cl)), cl, sum)
    data.frame(
      cluster = names(outside),
      outside = as.integer(outside),
      n       = as.integer(ntot),
      pct_out = round(100 * as.integer(outside) / as.integer(ntot), 1),
      row.names = NULL
    )
  }
  
  # ── UMAP clusters in PCA space, seed = central UMAP observation ───
  SEED_METHOD <- "nearest_mean"   # or "medoid"
  
  plots_umap_in_pca <- lapply(k_values, function(k) {
    sid <- umap_seed_ids(umap_df, all_cluster_aligned, k, method = SEED_METHOD)
    
    oc <- voronoi_outside_count(pca_df, all_cluster_aligned, k, sid,
                                cl_prefix = "km_umap_", xcol = "PC1", ycol = "PC2")
    cat(sprintf("\n--- UMAP clusters in PCA space, K = %d (seed: %s) ---\n",
                k, SEED_METHOD))
    print(oc)
    cat(sprintf("Total outside their cell: %d / %d (%.1f%%)\n",
                sum(oc$outside), sum(oc$n), 100*sum(oc$outside)/sum(oc$n)))
    
    plot_voronoi_aligned(pca_df, all_cluster_aligned, k,
                         cl_prefix = "km_umap_",
                         xcol = "PC1", ycol = "PC2",
                         title = "UMAP clusters in PCA space",
                         fixed = FALSE,
                         cell_text  = cell_texts_umap[[paste0("k", k)]],
                         medoid_ids = sid)
  })
  names(plots_umap_in_pca) <- paste0("k", k_values)
  
  # save
  out_dirs$umap.in.pca <- file.path(base_dir, "UMAP", "in_pca")
  if (!dir.exists(out_dirs$umap.in.pca)) dir.create(out_dirs$umap.in.pca, recursive = TRUE)
  save_vor_plots(plots_umap_in_pca, out_dirs$umap.in.pca)
  

  
# comparison of the 4 methods in PCA space ####
  
  library(patchwork)
  
  # reuses umap_seed_ids(): despite the name, it computes the seed in ANY
  # space given the right data frame. Here it is given the PCA one.
  # (it could be renamed central_obs_ids)
  
  # ── 4 lists, ALL in PCA space, one per method ─────────────────────
  build_in_pca <- function(cl_prefix, use_pam_medoids = FALSE) {
    lapply(k_values, function(k) {
      if (use_pam_medoids) {
        sid <- results_pam[[as.character(k)]]$id.med      # actual medoids (PAM)
      } else {
        # observation closest to the cluster mean, computed in PCA
        cl  <- factor(all_cluster_aligned[[paste0(cl_prefix, k)]])
        sid <- integer(0)
        for (lv in levels(cl)) {
          idx <- which(cl == lv)
          cx <- mean(pca_df$PC1[idx]); cy <- mean(pca_df$PC2[idx])
          d2 <- (pca_df$PC1[idx] - cx)^2 + (pca_df$PC2[idx] - cy)^2
          sid <- c(sid, idx[which.min(d2)])
        }
      }
      plot_voronoi_aligned(pca_df, all_cluster_aligned, k,
                           cl_prefix = cl_prefix,
                           xcol = "PC1", ycol = "PC2",
                           title = "", fixed = FALSE,
                           cell_text  = NULL,         # no text: small panel
                           medoid_ids = sid)
    })
  }
  
  plots_umap_pca <- build_in_pca("km_umap_")
  plots_pc2_pca  <- build_in_pca("km_pc2_k")
  plots_pc5_pca  <- build_in_pca("km_pc5_k")
  plots_pam_pca4 <- build_in_pca("pam_", use_pam_medoids = TRUE)
  
  names(plots_umap_pca) <- names(plots_pc2_pca) <-
    names(plots_pc5_pca) <- names(plots_pam_pca4) <- paste0("k", k_values)
  
  # ── 2x2 panel for one k ────────────────────────────────────
  make_4panel_pca <- function(k) {
    nm <- paste0("k", k)
    p1 <- plots_umap_pca[[nm]] + labs(subtitle = "KM-UMAP")
    p2 <- plots_pc2_pca[[nm]]  + labs(subtitle = "KM-PC2")
    p3 <- plots_pc5_pca[[nm]]  + labs(subtitle = "KM-PC5")
    p4 <- plots_pam_pca4[[nm]] + labs(subtitle = "PAM (Gower, medoids)")
    
    ((p1 | p2) / (p3 | p4)) +
      plot_layout(guides = "collect") +
      plot_annotation(
        title = sprintf("Cluster solutions in PCA space — K = %d", k),
        theme = theme(plot.title = element_text(face = "bold", size = 13))
      )
  }
  
  plots_4panel_pca <- lapply(k_values, make_4panel_pca)
  names(plots_4panel_pca) <- paste0("k", k_values)
  
  # ── Save ──────────────────────────────────────────────────────
  out_4pca <- file.path(base_dir, "COMPARISON_4METHODS", "in_pca")
  if (!dir.exists(out_4pca)) dir.create(out_4pca, recursive = TRUE)
  
  for (nm in names(plots_4panel_pca)) {
    ggsave(file.path(out_4pca, paste0(nm, ".pdf")),
           plots_4panel_pca[[nm]], width = 14, height = 11)
  }
  cat("Saved", length(plots_4panel_pca), "pannelli in", out_4pca, "\n")
  
