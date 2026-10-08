# =============================================================================
# scripts/03_structure/06_hdbscan_report.R
# Purpose : HDBSCAN results, console report and interpretation. Loads the
#           objects saved by 02_clustering/05_hdbscan.R and prints the minPts
#           sweep, the primary run (weak membership < 0.5), the ARI between
#           partitions across embedding dimensions and cluster profiles.
#           Nothing here recomputes clustering; it only reads and summarises.
# Thesis  : Appendix A.4; Tables A.8, A.9
# Input   : results/objects/HDBSCAN/{umap_embeddings,hdbscan_objects,
#           hdbscan_partitions,hdbscan_minpts_sweep}.rds
#           data/processed/data_imputed_selected.rds
# Output  : console only
# =============================================================================

source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages(library(dbscan))

# -- 0. Paths (same out_dir used in the pipeline) ------------------------------
out_dir <- file.path(DIR_OBJ, "HDBSCAN")
umap_embeddings    <- readRDS(file.path(out_dir, "umap_embeddings.rds"))
hdbscan_objects    <- readRDS(file.path(out_dir, "hdbscan_objects.rds"))
hdbscan_partitions <- readRDS(file.path(out_dir, "hdbscan_partitions.rds"))
minpts_sweep       <- readRDS(file.path(out_dir, "hdbscan_minpts_sweep.rds"))

umap_dims <- as.integer(sub("dim", "", names(umap_embeddings)))
minPts    <- hdbscan_objects[[1]]$minPts        # primary minPts of the objects

# original standardized data for substantive interpretation
X <- readRDS(file.path(DIR_DATA, "data_imputed_selected.rds")) %>%
  dplyr::mutate(dplyr::across(where(is.factor),
                              ~ as.numeric(as.character(.)))) %>%
  scale()
stopifnot(nrow(X) == nrow(hdbscan_partitions))

# ── printing helpers ──────────────────────────────────────────────── #####
hr  <- function(ch = "-", n = 72) cat(strrep(ch, n), "\n", sep = "")
sec <- function(txt) { cat("\n"); hr("="); cat(txt, "\n"); hr("=") }
sub <- function(txt) { cat("\n"); cat(txt, "\n"); hr("-") }

# base-R Adjusted Rand Index (no extra dependency)
adj_rand <- function(a, b) {
  tab   <- table(a, b)
  comb2 <- function(x) sum(choose(x, 2))
  ai <- rowSums(tab); bj <- colSums(tab); n <- sum(tab)
  idx  <- comb2(tab)
  exp_ <- comb2(ai) * comb2(bj) / choose(n, 2)
  max_ <- (comb2(ai) + comb2(bj)) / 2
  (idx - exp_) / (max_ - exp_)
}

# per-cluster profile on the original (standardized) variables
cluster_profile <- function(labels, Xmat, top = 12, drop_noise = TRUE) {
  top <- min(top, ncol(Xmat))
  cls <- sort(unique(labels[if (drop_noise) labels != 0 else TRUE]))
  for (cc in cls) {
    idx <- which(labels == cc)
    mu  <- colMeans(Xmat[idx, , drop = FALSE])          # mean z-score per var
    ord <- order(abs(mu), decreasing = TRUE)[seq_len(top)]
    cat(sprintf("\n  Cluster %s  (n = %d)  — top %d distinctive variables\n",
                cc, length(idx), top))
    prof <- data.frame(variable = colnames(Xmat)[ord],
                       mean_z   = round(mu[ord], 2),
                       dir      = ifelse(mu[ord] >= 0, "high", "low"))
    print(prof, row.names = FALSE)
  }
}

# ═══════════════════════════════════════════════════════════════════════
# 1. minPts SWEEP — how the emergent structure moves with minPts
# ═══════════════════════════════════════════════════════════════════════
sec("1. minPts SWEEP  (emergent n_clusters and noise per UMAP dimension)")

cat("\nn_clusters  (rows = minPts, cols = UMAP dim):\n")
print(xtabs(n_clusters ~ minPts + dim, data = minpts_sweep))

cat("\nnoise fraction:\n")
print(round(xtabs(noise_frac ~ minPts + dim, data = minpts_sweep), 3))

cat("\nFull long table:\n")
print(transform(minpts_sweep, noise_frac = round(noise_frac, 3)),
      row.names = FALSE)

# heuristic hint (NOT an automatic selection): fewest-noise runs with >=2 clusters
sub("Low-noise configurations with >= 2 clusters (heuristic reference only)")
hint <- minpts_sweep %>%
  dplyr::filter(n_clusters >= 2) %>%
  dplyr::group_by(dim) %>%
  dplyr::slice_min(noise_frac, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(noise_frac = round(noise_frac, 3))
print(as.data.frame(hint), row.names = FALSE)

# ═══════════════════════════════════════════════════════════════════════
# 2. PRIMARY RUN per UMAP dimension  (minPts = primary)
# ═══════════════════════════════════════════════════════════════════════
sec(sprintf("2. PRIMARY HDBSCAN RUN  (minPts = %d)", minPts))

for (d in umap_dims) {
  res    <- hdbscan_objects[[paste0("dim", d)]]
  lab    <- res$cluster
  ncl    <- length(unique(lab[lab != 0]))
  nnoise <- sum(lab == 0)

  sub(sprintf("UMAP %dd  |  %d clusters  |  noise = %d (%.1f%%)",
              d, ncl, nnoise, 100 * nnoise / length(lab)))

  cat("Cluster sizes (0 = noise):\n")
  print(table(lab))

  cat("\nStability / persistence score per cluster:\n")
  sc <- res$cluster_scores
  if (length(sc)) print(round(sort(sc, decreasing = TRUE), 3)) else cat("  (none)\n")

  cat("\nMembership probability (assignment confidence):\n")
  mp_by <- tapply(res$membership_prob, lab, function(v) round(mean(v), 3))
  print(mp_by)
  cat(sprintf("  points with prob < 0.5 (weak assignment): %d (%.1f%%)\n",
              sum(res$membership_prob < 0.5),
              100 * mean(res$membership_prob < 0.5)))

  cat("\nGLOSH outlier score summary:\n")
  print(round(summary(res$outlier_scores), 3))
  topo <- order(res$outlier_scores, decreasing = TRUE)[1:5]
  cat("  top-5 outliers (row -> GLOSH, cluster):\n")
  print(data.frame(row     = topo,
                   glosh   = round(res$outlier_scores[topo], 3),
                   cluster = lab[topo]),
        row.names = FALSE)
}

# ═══════════════════════════════════════════════════════════════════════
# 3. CONSISTENCY ACROSS UMAP DIMENSIONS  (Adjusted Rand Index)
# ═══════════════════════════════════════════════════════════════════════
sec("3. AGREEMENT BETWEEN DIMENSIONS  (ARI, primary minPts, noise incl.)")

labs_list <- lapply(umap_dims, function(d) hdbscan_objects[[paste0("dim", d)]]$cluster)
ari <- matrix(NA_real_, length(umap_dims), length(umap_dims),
              dimnames = list(paste0("d", umap_dims), paste0("d", umap_dims)))
for (i in seq_along(umap_dims))
  for (j in seq_along(umap_dims))
    ari[i, j] <- adj_rand(labs_list[[i]], labs_list[[j]])
cat("\nARI matrix (1 = identical partitions):\n")
print(round(ari, 3))
cat("\nReading: high off-diagonal ARI => the emergent partition is robust to\n",
    "the UMAP dimensionality; low values => strongly dimension-dependent.\n", sep = "")

# ═══════════════════════════════════════════════════════════════════════
# 4. SUBSTANTIVE INTERPRETATION  (what each cluster means)
#    per-cluster mean z-score on the ORIGINAL variables; the variables with
#    the largest |mean z-score| are the ones that most define the cluster.
# ═══════════════════════════════════════════════════════════════════════
sec("4. CLUSTER PROFILES ON ORIGINAL VARIABLES  (primary minPts)")

for (d in umap_dims) {
  sub(sprintf("UMAP %dd", d))
  cluster_profile(hdbscan_objects[[paste0("dim", d)]]$cluster, X,
                  top = 12, drop_noise = TRUE)
}

cat("\n")
hr("=")
cat("Notes:\n")
cat(" - mean_z is the cluster mean on standardized variables; 0 = overall\n")
cat("   average, so |mean_z| ranks how distinctive a variable is per cluster.\n")
cat(" - Profiles here are descriptive; cross-check with your RF / LASSO /\n")
cat("   OVR characterization for a supervised view of the same clusters.\n")
cat(" - Cluster 0 is HDBSCAN noise (points in no dense region), not a group.\n")
hr("=")

