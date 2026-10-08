# =====================================================================
# WITHIN-CLUSTER SHAPE & GOWER-DISTANCE DIAGNOSTICS
# ---------------------------------------------------------------------
# scripts/06_robustness/14_within_cluster_diagnostics.R
# Thesis  : not reported (within-cluster shape diagnostics)
# Input   : results/objects/clusters/{all_clusters_aligned,gower_dist,
#           umap_df_clusters}.rds, results/objects/pca_scores.rds
# Output  : results/objects/diagnostics/*.rds,
#           results/figures/06_robustness/cluster_diagnostics/*.pdf
#
# Two questions:
#   (A) Is each cluster concentrated / unimodal / normal-like around its
#       own mean?  -> Mahalanobis D2 vs chi-square, Mardia, Shapiro,
#                     Hartigan dip test, marginal densities.
#   (B) How are Gower dissimilarities distributed INSIDE each cluster?
#       -> full within-pair distribution, within vs between, distance to
#          the cluster medoid.  Reported for k = 2, 4, 5.
#
# Labels are read exclusively from all_clusters_aligned (authoritative).
# All Gower statistics use the FULL set of pairs (no subsampling);
# subsampling is applied ONLY to the vectors kept for plotting.
# =====================================================================

set.seed(42)
source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})
has_dip  <- requireNamespace("diptest", quietly = TRUE)
has_MASS <- requireNamespace("MASS",    quietly = TRUE)
if (!has_dip) message("NOTE: package 'diptest' not installed -> dip test skipped.")

# ── 0. CONFIG ────────────────────────────────────────────────────────
# ── PATHS ────────────────────────────────────────────────────────────
OBJ_DIR <- DIR_OBJ
CLU_DIR <- DIR_CLUST

# Objects go under results/objects, plots under results/figures.
OUT_OBJ   <- ensure_dir(DIR_OBJ, "diagnostics")
OUT_PLOTS <- ensure_dir(DIR_FIG, "06_robustness", "cluster_diagnostics")
cat("OUT_OBJ  :", OUT_OBJ,   "\nOUT_PLOTS:", OUT_PLOTS, "\n")

# PATH_CLUSTERS comes from R/00_setup.R (canonical cluster labels).
PATH_GOWER  <- file.path(DIR_CLUST, "gower_dist.rds")
PATH_UMAP   <- file.path(DIR_CLUST, "umap_df_clusters.rds")
PATH_PCA    <- file.path(DIR_OBJ,   "pca_scores.rds")
PATH_PCA_DF <- file.path(DIR_CLUST, "pca_df_clusters.rds")     # fallback (PC1..PC3 only)

cat("Resolved paths:\n")
print(c(clusters = PATH_CLUSTERS, gower = PATH_GOWER,
        umap = PATH_UMAP, pca_scores = PATH_PCA, pca_df = PATH_PCA_DF))

stopifnot(!is.na(PATH_CLUSTERS), !is.na(PATH_GOWER), !is.na(PATH_UMAP))

K_TARGET <- c(2, 4, 5)          # k values analysed for the Gower section

# Statistics run on ALL pairs. This cap only thins the plotting data.
MAX_PLOT_POINTS <- 100000       # set to Inf to keep every pair (heavy RAM)


# ── 1. OBJECT RESOLUTION: memory first, disk as fallback ─────────────
# Interactive RStudio sessions keep objects alive; reuse them instead of
# re-reading, but log the provenance so stale objects are visible.
.resolve <- function(obj_name, path, validator = NULL) {
  if (exists(obj_name, envir = .GlobalEnv, inherits = FALSE)) {
    x <- get(obj_name, envir = .GlobalEnv)
    ok <- is.null(validator) || isTRUE(validator(x))
    if (ok) {
      cat(sprintf("  [memory] %-22s\n", obj_name))
      return(x)
    }
    cat(sprintf("  [memory->REJECTED] %-22s (failed validation, reloading)\n",
                obj_name))
  }
  if (!file.exists(path)) stop("Object '", obj_name, "' not in memory and file not found: ", path)
  cat(sprintf("  [disk]   %-22s <- %s\n", obj_name, basename(path)))
  readRDS(path)
}

cat("Resolving objects:\n")

all_clusters_aligned <- .resolve(
  "all_clusters_aligned", PATH_CLUSTERS,
  validator = function(x) is.data.frame(x) && any(grepl("_\\d+$", names(x)))
)

# Gower: accept an already-materialised matrix (G) or a dist/daisy object.
if (exists("G", envir = .GlobalEnv, inherits = FALSE) &&
    is.matrix(get("G", envir = .GlobalEnv))) {
  cat("  [memory] G (Gower matrix)\n")
  G <- get("G", envir = .GlobalEnv)
} else {
  gower_dist <- .resolve(
    "gower_dist", PATH_GOWER,
    validator = function(x) inherits(x, "dist") || inherits(x, "dissimilarity") || is.matrix(x)
  )
  G <- if (is.matrix(gower_dist)) gower_dist else as.matrix(gower_dist)
}
storage.mode(G) <- "double"

pca_scores_raw <- .resolve(
  "pca_scores", PATH_PCA,
  validator = function(x) inherits(x, "prcomp") || is.data.frame(x) || is.matrix(x)
)

pca_scores <- if (inherits(pca_scores_raw, "prcomp")) {
  as.data.frame(pca_scores_raw$x)
} else {
  as.data.frame(pca_scores_raw)
}

stopifnot(is.data.frame(pca_scores), ncol(pca_scores) >= 5)
umap_raw <- .resolve(
  "umap_df", PATH_UMAP,
  validator = function(x) is.data.frame(x) && all(c("UMAP1", "UMAP2") %in% names(x))
)
umap_coords <- as.matrix(umap_raw[, c("UMAP1", "UMAP2"), drop = FALSE])

N <- nrow(all_clusters_aligned)
stopifnot(nrow(pca_scores)  == N,
          nrow(umap_coords) == N,
          nrow(G)           == N,
          ncol(G)           == N)
print(cat(sprintf("\nN = %d observations, PCA dims available = %d\n", N, ncol(pca_scores))))

# ── 2. DETECT METHOD x K COLUMNS ─────────────────────────────────────
# Two naming conventions coexist in all_clusters_aligned:
#   "km_umap_5", "pam_5"   -> suffix "_<k>"
#   "km_pc5_k5"            -> suffix "_k<k>"
# The pattern below accepts both.
PAT_K <- "_k?\\d+$"

cl_cols <- grep(PAT_K, names(all_clusters_aligned), value = TRUE)

cl_index <- data.frame(
  col    = cl_cols,
  method = sub(PAT_K, "", cl_cols),
  k      = as.integer(sub(paste0("^.*", PAT_K), "\\1",
                          sub("^.*_k?(\\d+)$", "\\1", cl_cols))),
  stringsAsFactors = FALSE
)
# simpler and safer: extract k directly
cl_index$k <- as.integer(sub("^.*_k?(\\d+)$", "\\1", cl_cols))

stopifnot(!anyNA(cl_index$k))

cat("\nAll label columns detected:\n")
print(table(cl_index$method))

# Keep only the four methods used in the thesis, and the target k values
METHODS_KEEP <- c("km_pc2", "km_pc5", "km_umap", "pam")

cl_index <- cl_index[cl_index$method %in% METHODS_KEEP &
                       cl_index$k      %in% K_TARGET, ]
cl_index <- cl_index[order(match(cl_index$method, METHODS_KEEP), cl_index$k), ]

cat("\nMethod x k cells analysed:\n")
print(cl_index, row.names = FALSE)
stopifnot(nrow(cl_index) == length(METHODS_KEEP) * length(K_TARGET))
cat("\nLabel columns detected in all_clusters_aligned:\n")
print(cl_index, row.names = FALSE)

cl_index <- cl_index[cl_index$k %in% K_TARGET, ]
stopifnot(nrow(cl_index) > 0)

# ── 3. MDS EMBEDDING OF THE GOWER SPACE (for PAM) ────────────────────
# PAM has no native Euclidean space; a 5-axis classical MDS gives an
# approximate one in which "distance from the mean" is meaningful.
if (exists("mds_coords_diag", envir = .GlobalEnv, inherits = FALSE) &&
    is.matrix(get("mds_coords_diag", envir = .GlobalEnv)) &&
    nrow(get("mds_coords_diag", envir = .GlobalEnv)) == N) {
  cat("  [memory] mds_coords_diag\n")
  mds_coords <- get("mds_coords_diag", envir = .GlobalEnv)
} else {
  cat("  [compute] classical MDS on Gower (5 axes)...\n")
  mds_fit_diag <- stats::cmdscale(G, k = 5, eig = TRUE)
  mds_coords   <- mds_fit_diag$points
  colnames(mds_coords) <- paste0("MDS", seq_len(ncol(mds_coords)))
  mds_coords_diag <- mds_coords          # cached for re-runs
}

# ── 4. METHOD -> EMBEDDING MAP (edit here if prefixes differ) ────────
get_embedding <- function(method) {
  switch(
    method,
    "km_pc2"  = as.matrix(pca_scores[, 1:2, drop = FALSE]),
    "km_pc5"  = as.matrix(pca_scores[, 1:5, drop = FALSE]),
    "km_umap" = umap_coords,
    "pam"     = mds_coords,        # approximate: MDS(Gower)
    NULL
  )
}

# ── 5. SHAPE HELPERS ─────────────────────────────────────────────────
mardia <- function(Xc) {
  n <- nrow(Xc); p <- ncol(Xc)
  S <- stats::cov(Xc) * (n - 1) / n
  Sinv <- tryCatch(solve(S), error = function(e) {
    if (has_MASS) MASS::ginv(S) else stop("Singular covariance; install MASS.")
  })
  Z <- scale(Xc, center = TRUE, scale = FALSE)
  D <- Z %*% Sinv %*% t(Z)
  
  b1p <- sum(D^3) / n^2
  b2p <- sum(diag(D)^2) / n
  
  stat_skew <- n * b1p / 6
  df_skew   <- p * (p + 1) * (p + 2) / 6
  p_skew    <- stats::pchisq(stat_skew, df_skew, lower.tail = FALSE)
  
  z_kurt <- (b2p - p * (p + 2)) / sqrt(8 * p * (p + 2) / n)
  p_kurt <- 2 * stats::pnorm(abs(z_kurt), lower.tail = FALSE)
  
  c(mardia_skew_p = p_skew, mardia_kurt = b2p,
    mardia_kurt_z = z_kurt, mardia_kurt_p = p_kurt)
}

shape_one_cluster <- function(X, idx) {
  Xc <- X[idx, , drop = FALSE]
  n  <- nrow(Xc); p <- ncol(Xc)
  
  d2 <- stats::mahalanobis(Xc, center = colMeans(Xc),
                           cov = stats::cov(Xc), tol = 1e-20)
  
  md <- if (n > p + 5) mardia(Xc) else c(mardia_skew_p = NA_real_,
                                         mardia_kurt   = NA_real_,
                                         mardia_kurt_z = NA_real_,
                                         mardia_kurt_p = NA_real_)
  
  # Shapiro-Wilk is capped at n = 5000 by stats::shapiro.test
  ii <- if (n > 5000) sample(seq_len(n), 5000) else seq_len(n)
  shap <- suppressWarnings(
    vapply(seq_len(p), function(j) stats::shapiro.test(Xc[ii, j])$p.value, numeric(1))
  )
  dipp <- if (has_dip) {
    vapply(seq_len(p), function(j) diptest::dip.test(Xc[, j])$p.value, numeric(1))
  } else rep(NA_real_, p)
  
  list(
    d2     = d2,
    coords = Xc,
    stats  = data.frame(
      n             = n,
      p             = p,
      # share of points inside the 95% Gaussian ellipsoid:
      #   ~95  = Gaussian-like
      #   >95  = more concentrated than Gaussian (Voronoi truncation)
      #   <95  = heavy tails
      pct_chisq95   = 100 * mean(d2 <= stats::qchisq(0.95, p)),
      mardia_skew_p = unname(md["mardia_skew_p"]),
      mardia_kurt_z = unname(md["mardia_kurt_z"]),
      mardia_kurt_p = unname(md["mardia_kurt_p"]),
      shapiro_min_p = min(shap),
      dip_min_p     = if (all(is.na(dipp))) NA_real_ else min(dipp)
    )
  )
}

# ── 6. SHAPE DIAGNOSTICS LOOP (sequential: shared server) ────────────
shape_tab  <- list()
qq_plots   <- list()
marg_plots <- list()

for (i in seq_len(nrow(cl_index))) {
  col    <- cl_index$col[i]
  method <- cl_index$method[i]
  k      <- cl_index$k[i]
  
  X <- get_embedding(method)
  if (is.null(X)) {
    message("No embedding defined for method '", method, "' -> skipped.")
    next
  }
  
  lab <- as.integer(as.character(all_clusters_aligned[[col]]))
  stopifnot(length(lab) == N, !anyNA(lab))
  
  qq_df <- list(); mg_df <- list()
  
  for (cl in sort(unique(lab))) {
    idx <- which(lab == cl)
    if (length(idx) < ncol(X) + 6) next
    
    res <- shape_one_cluster(X, idx)
    
    shape_tab[[length(shape_tab) + 1]] <- cbind(
      data.frame(method = method, k = k, cluster = cl,
                 share_pct = round(100 * length(idx) / N, 1)),
      res$stats
    )
    
    p <- ncol(X); n <- length(idx)
    qq_df[[as.character(cl)]] <- data.frame(
      cluster = factor(cl),
      theo    = stats::qchisq(stats::ppoints(n), df = p),
      obs     = sort(res$d2)
    )
    
    mg_df[[as.character(cl)]] <- as.data.frame(res$coords) %>%
      tidyr::pivot_longer(dplyr::everything(),
                          names_to = "dim", values_to = "value") %>%
      dplyr::mutate(cluster = factor(cl))
  }
  
  qq_all <- dplyr::bind_rows(qq_df)
  mg_all <- dplyr::bind_rows(mg_df)
  if (nrow(qq_all) == 0) next
  
  # (a) chi-square QQ-plot of Mahalanobis D2
  qq_plots[[paste0(method, "_k", k)]] <-
    ggplot(qq_all, aes(theo, obs, color = cluster)) +
    geom_abline(slope = 1, intercept = 0, color = "grey40", linetype = 2) +
    geom_point(size = 0.6, alpha = 0.6) +
    scale_color_manual(values = cluster_colors, guide = "none") +
    facet_wrap(~ cluster, scales = "free", nrow = 2) +
    labs(
      title    = sprintf("Chi-square QQ-plot of Mahalanobis D2 -- %s, k = %d",
                         method, k),
      subtitle = "On the dashed line = the cluster is multivariate normal around its own mean",
      x = expression(paste("theoretical ", chi^2, " quantiles")),
      y = expression(paste("observed ", D^2))
    ) +
    theme_minimal(base_size = 10)
  
  # (b) marginal densities with a fitted Normal overlay
  norm_curve <- mg_all %>%
    dplyr::group_by(cluster, dim) %>%
    dplyr::summarise(m  = mean(value), s = stats::sd(value),
                     lo = min(value),  hi = max(value), .groups = "drop") %>%
    dplyr::rowwise() %>%
    dplyr::do({
      xs <- seq(.$lo, .$hi, length.out = 200)
      data.frame(cluster = .$cluster, dim = .$dim,
                 value = xs, dens = stats::dnorm(xs, .$m, .$s))
    }) %>%
    dplyr::ungroup()
  
  marg_plots[[paste0(method, "_k", k)]] <-
    ggplot(mg_all, aes(value)) +
    geom_density(aes(color = cluster), linewidth = 0.6) +
    geom_line(data = norm_curve, aes(y = dens), color = "grey30",
              linetype = 2, linewidth = 0.4) +
    scale_color_manual(values = cluster_colors, guide = "none") +
    facet_grid(cluster ~ dim, scales = "free") +
    labs(
      title    = sprintf("Within-cluster marginal densities -- %s, k = %d", method, k),
      subtitle = "Dashed grey = Normal(mean, sd) fitted on the same cluster",
      x = NULL, y = "density"
    ) +
    theme_minimal(base_size = 9)
  
  rm(X, lab, qq_df, mg_df, qq_all, mg_all, norm_curve)
}

shape_table <- dplyr::bind_rows(shape_tab) %>%
  dplyr::arrange(method, k, cluster)
cat("\n--- SHAPE DIAGNOSTICS ---\n")
print(shape_table, row.names = FALSE, digits = 3)

# ── 7. GOWER DISTANCES: FULL within / between / medoid ───────────────
within_pairs <- function(idx) {
  sub <- G[idx, idx]
  sub[upper.tri(sub)]              # all n(n-1)/2 within-cluster pairs
}
between_pairs <- function(idx) {
  as.vector(G[idx, -idx])          # all n * (N - n) cluster-vs-rest pairs
}
thin <- function(v, m = MAX_PLOT_POINTS) {
  if (is.finite(m) && length(v) > m) sample(v, m) else v
}
skew1 <- function(x) mean((x - mean(x))^3) / stats::sd(x)^3

gower_long  <- list()
medoid_long <- list()
gower_tab   <- list()

for (i in seq_len(nrow(cl_index))) {
  col    <- cl_index$col[i]
  method <- cl_index$method[i]
  k      <- cl_index$k[i]
  lab    <- as.integer(as.character(all_clusters_aligned[[col]]))
  
  for (cl in sort(unique(lab))) {
    idx <- which(lab == cl)
    if (length(idx) < 5) next
    
    d_within <- within_pairs(idx)    # FULL - no subsampling
    d_betw   <- between_pairs(idx)   # FULL - no subsampling
    
    # medoid = observation minimising total within-cluster Gower distance
    # (defined for k-means clusters too, not only PAM)
    sub      <- G[idx, idx]
    med_pos  <- which.min(colSums(sub))
    med_id   <- idx[med_pos]
    d_medoid <- sub[, med_pos]       # n values: kept in full
    
    dw_plot <- thin(d_within)
    db_plot <- thin(d_betw)
    
    gower_long[[length(gower_long) + 1]] <- data.frame(
      method  = method, k = k, cluster = factor(cl),
      type    = c(rep("within",  length(dw_plot)),
                  rep("between", length(db_plot))),
      dist    = c(dw_plot, db_plot)
    )
    medoid_long[[length(medoid_long) + 1]] <- data.frame(
      method = method, k = k, cluster = factor(cl), dist = d_medoid
    )
    
    gower_tab[[length(gower_tab) + 1]] <- data.frame(
      method = method, k = k, cluster = cl, n = length(idx),
      share_pct        = round(100 * length(idx) / N, 1),
      n_within_pairs   = length(d_within),   # audit: full set, no sampling
      n_between_pairs  = length(d_betw),
      within_mean      = mean(d_within),
      within_median    = stats::median(d_within),
      within_sd        = stats::sd(d_within),
      within_iqr       = stats::IQR(d_within),
      within_skew      = skew1(d_within),
      between_mean     = mean(d_betw),
      # < 1 = cluster more homogeneous than its complement
      # ~ 1 = cluster is no more homogeneous than the sample at large
      sep_ratio        = mean(d_within) / mean(d_betw),
      medoid_id        = med_id,
      medoid_dist_mean = mean(d_medoid),
      medoid_dist_sd   = stats::sd(d_medoid),
      medoid_dist_skew = skew1(d_medoid),
      medoid_dip_p     = if (has_dip) diptest::dip.test(d_medoid)$p.value else NA_real_
    )
    
    rm(d_within, d_betw, dw_plot, db_plot, sub, d_medoid)  # keep RAM flat
  }
  rm(lab)
}

gower_df    <- dplyr::bind_rows(gower_long)
medoid_df   <- dplyr::bind_rows(medoid_long)
gower_table <- dplyr::bind_rows(gower_tab) %>% dplyr::arrange(method, k, cluster)

cat("\n--- GOWER WITHIN-CLUSTER DIAGNOSTICS (all pairs) ---\n")
print(gower_table, row.names = FALSE, digits = 3)

global_mean_gower <- mean(G[upper.tri(G)])
cat(sprintf("\nGlobal mean Gower distance = %.4f\n", global_mean_gower))

# ── 8. GOWER PLOTS (one figure per k, faceted by method) ─────────────
gower_plots  <- list()
medoid_plots <- list()
box_plots    <- list()

plot_caption <- if (is.finite(MAX_PLOT_POINTS)) {
  sprintf(paste("Densities estimated on a random subsample of up to %s pairs;",
                "all reported statistics use the full set of pairs."),
          format(MAX_PLOT_POINTS, big.mark = ","))
} else "All pairs used."

for (kk in K_TARGET) {
  d_k <- dplyr::filter(gower_df,  k == kk)
  m_k <- dplyr::filter(medoid_df, k == kk)
  if (nrow(d_k) == 0) next
  
  gower_plots[[paste0("k", kk)]] <-
    ggplot(d_k, aes(dist, color = cluster, linetype = type)) +
    geom_density(linewidth = 0.6) +
    geom_vline(xintercept = global_mean_gower, color = "grey40",
               linetype = 3, linewidth = 0.4) +
    scale_color_manual(values = cluster_colors, name = "Cluster") +
    scale_linetype_manual(values = c(within = "solid", between = "22"),
                          name = "Pairs") +
    facet_wrap(~ method, scales = "free_y") +
    labs(
      title    = sprintf("Gower distance distributions -- k = %d", kk),
      subtitle = "Solid = within-cluster pairs; dashed = cluster vs rest; dotted grey = global mean",
      caption  = plot_caption,
      x = "Gower distance", y = "density"
    ) +
    theme_minimal(base_size = 10)
  
  medoid_plots[[paste0("k", kk)]] <-
    ggplot(m_k, aes(dist, color = cluster, fill = cluster)) +
    geom_density(alpha = 0.15, linewidth = 0.6) +
    scale_color_manual(values = cluster_colors, name = "Cluster") +
    scale_fill_manual(values = cluster_colors, guide = "none") +
    facet_wrap(~ method, scales = "free_y") +
    labs(
      title    = sprintf("Gower distance to the cluster medoid -- k = %d", kk),
      subtitle = "Unimodal & right-skewed = single dense core; bimodal = the cluster cuts through a continuum",
      x = "Gower distance to medoid", y = "density"
    ) +
    theme_minimal(base_size = 10)
  
  box_plots[[paste0("k", kk)]] <-
    ggplot(dplyr::filter(d_k, type == "within"),
           aes(cluster, dist, fill = cluster)) +
    geom_boxplot(outlier.size = 0.3, alpha = 0.6) +
    geom_hline(yintercept = global_mean_gower, color = "grey40", linetype = 3) +
    scale_fill_manual(values = cluster_colors, guide = "none") +
    facet_wrap(~ method) +
    labs(
      title    = sprintf("Within-cluster Gower distances -- k = %d", kk),
      subtitle = "Dotted line = global mean Gower distance",
      caption  = plot_caption,
      x = "Cluster", y = "Gower distance"
    ) +
    theme_minimal(base_size = 10)
}

# ── 9. SAVINGS (all writes consolidated here) ────────────────────────
saveRDS(shape_table, file.path(OUT_OBJ, "shape_diagnostics_table.rds"))
saveRDS(gower_table, file.path(OUT_OBJ, "gower_within_table.rds"))
saveRDS(gower_df,    file.path(OUT_OBJ, "gower_distance_long.rds"))
saveRDS(medoid_df,   file.path(OUT_OBJ, "gower_medoid_distance_long.rds"))

write.csv(shape_table, file.path(OUT_OBJ, "shape_diagnostics_table.csv"),
          row.names = FALSE)
write.csv(gower_table, file.path(OUT_OBJ, "gower_within_table.csv"),
          row.names = FALSE)

for (nm in names(qq_plots))
  ggsave(file.path(OUT_PLOTS, paste0("qq_chisq_", nm, ".pdf")),
         qq_plots[[nm]], width = 10, height = 6)

for (nm in names(marg_plots))
  ggsave(file.path(OUT_PLOTS, paste0("marginal_density_", nm, ".pdf")),
         marg_plots[[nm]], width = 11, height = 7)

for (nm in names(gower_plots))
  ggsave(file.path(OUT_PLOTS, paste0("gower_within_between_", nm, ".pdf")),
         gower_plots[[nm]], width = 11, height = 7)

for (nm in names(medoid_plots))
  ggsave(file.path(OUT_PLOTS, paste0("gower_medoid_density_", nm, ".pdf")),
         medoid_plots[[nm]], width = 11, height = 7)

for (nm in names(box_plots))
  ggsave(file.path(OUT_PLOTS, paste0("gower_within_boxplot_", nm, ".pdf")),
         box_plots[[nm]], width = 11, height = 7)

cat("\nDone.\n  Plots   ->", OUT_PLOTS,
    "\n  Objects ->", OUT_OBJ, "\n")
