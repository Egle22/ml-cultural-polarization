# =============================================================================
# scripts/03_structure/05_silhouette.R
#   Cluster model selection via a hedonic-game preference analogy
#   GLOBAL variant: joint choice of (method, k) across all candidate
#   partitions, aggregated by majority vote over respondents.
#
#   Preference of respondent i for a candidate partition = silhouette width
#   s_i in that partition, where
#      a_i = mean Gower dissimilarity to co-members  (cost of staying)
#      b_i = mean Gower dissimilarity to nearest foreign cluster
#            (cost of the best individual deviation = Nash-blocking target)
#      s_i = (b_i - a_i) / max(a_i, b_i)   in [-1, 1]
#   Higher s_i  <=>  respondent i is better placed  <=>  partition preferred.
#   Rationale: raw cohesion a_i alone is monotonically biased toward high k;
#   the deviation baseline b_i removes that bias (see raw-cohesion diagnostic
#   reported at the end, kept ONLY to expose the bias, not to select).
# Thesis  : Section 4.1.4; Table A.6 (average silhouette by method and k)
# Input   : results/objects/clusters/{gower_dist,all_clusters_aligned}.rds
# Output  : results/objects/preference_voting/{silhouette_matrix,
#           raw_cohesion_matrix,person_top_choice_global,
#           person_pref_ranking_global,vote_silhouette_global,
#           mean_silhouette_per_partition}.rds
# Note    : deterministic.
# =============================================================================

source(here::here("R", "00_setup.R"))

# -- Paths ---------------------------------------------------------------------
clusters_dir <- DIR_CLUST
out_dir      <- ensure_dir(DIR_OBJ, "preference_voting")

# ── Load data (rm() guard against stale interactive objects) ───────────
rm(list = intersect(c("gower_dist", "all_clusters_aligned", "gower_mat"), ls()))

gower_dist           <- readRDS(file.path(clusters_dir, "gower_dist.rds"))
all_clusters_aligned <- readRDS(file.path(clusters_dir, "all_clusters_aligned.rds"))
stopifnot(
  "gower_dist must be a 'dissimilarity'/'dist' object" =
    inherits(gower_dist, c("dissimilarity", "dist")),
  "gower_dist size must match nrow(all_clusters_aligned)" =
    attr(gower_dist, "Size") == nrow(all_clusters_aligned)
)


## 1. keep only the wanted families
nm <- names(all_clusters_aligned)
keep_prefix <- grepl("^(km_umap_|pam_|km_pc2_k|km_pc5_k)", nm) | nm == "hdb_umap2_mp25"

## 2. extract trailing k (NA for hdb, which has no k)
k_num <- suppressWarnings(as.integer(sub("^.*?([0-9]+)$", "\\1", nm)))
k_num[nm == "hdb_umap2_mp25"] <- NA_integer_

## 3. keep k <= 6, plus hdbscan (no k)
keep_idx <- keep_prefix & (is.na(k_num) | k_num <= 6)

all_clusters_aligned <- all_clusters_aligned[keep_idx]



k_values <- 2:6

cols_keep <- c(
  paste0("km_pc5_k", k_values),
  paste0("km_umap_", k_values),
  paste0("km_pc2_k", k_values),
  paste0("pam_", k_values)
)

all_clusters_aligned <- all_clusters_aligned[, cols_keep]

gower_mat    <- as.matrix(gower_dist)          # 2579 x 2579, ~53 MB
n            <- nrow(all_clusters_aligned)
cluster_cols <- names(all_clusters_aligned)

# ── Raw within-cluster mean dissimilarity a_i (DIAGNOSTIC ONLY) ─────────
# Divides by (cluster size - 1) = mean dissimilarity to OTHER members,
# matching the silhouette convention. Singletons -> NA (undefined).
avg_within_cluster_dist <- function(labels, dist_mat) {
  out <- rep(NA_real_, length(labels))
  for (cl in unique(labels)) {
    idx <- which(labels == cl)
    if (length(idx) == 1L) next
    sub <- dist_mat[idx, idx, drop = FALSE]
    out[idx] <- rowSums(sub) / (length(idx) - 1L)
  }
  out
}

# ── Per-partition silhouette (rows returned in original respondent order) ──
silhouette_profile <- function(labels, dist_obj, dist_mat) {
  sil <- cluster::silhouette(labels, dist_obj)
  stopifnot("silhouette() did not return per-observation widths" =
              is.matrix(sil) && nrow(sil) == length(labels))
  data.frame(
    sil_width = as.numeric(sil[, "sil_width"]),
    neighbour = as.integer(sil[, "neighbor"]),
    a_raw     = avg_within_cluster_dist(labels, dist_mat)
  )
}

profiles <- lapply(cluster_cols, function(col)
  silhouette_profile(all_clusters_aligned[[col]], gower_dist, gower_mat))
names(profiles) <- cluster_cols

# ── Assemble matrices: respondents x candidate partitions ──────────────
sil_mat <- vapply(profiles, function(d) d$sil_width, numeric(n))
a_mat   <- vapply(profiles, function(d) d$a_raw,     numeric(n))
rownames(sil_mat) <- rownames(a_mat) <- rownames(all_clusters_aligned)

stopifnot(
  "silhouette matrix must be respondents x n_partitions" =
    all(dim(sil_mat) == c(n, length(cluster_cols)))
)

# ── Per-respondent most-preferred partition (argmax silhouette) ────────
top_choice <- function(row, decreasing = TRUE) {
  if (all(is.na(row))) return(NA_character_)
  names(row)[if (decreasing) which.max(row) else which.min(row)]
}

person_top_choice <- apply(sil_mat, 1, top_choice, decreasing = TRUE)

# Full preference profile per respondent (most -> least preferred):
# the empirical analogue of the ordinal preference lists in Example 1.
person_pref_ranking <- t(apply(sil_mat, 1, function(x)
  names(x)[order(x, decreasing = TRUE, na.last = TRUE)]))

# ── Majority vote (plurality) across respondents ───────────────────────
vote_silhouette   <- sort(table(person_top_choice), decreasing = TRUE)
winner_silhouette <- names(vote_silhouette)[1]

# ── Cross-checks / diagnostics ─────────────────────────────────────────
mean_sil_per_partition <- sort(colMeans(sil_mat, na.rm = TRUE), decreasing = TRUE)

# Raw-cohesion vote: kept ONLY to demonstrate the k-bias (lower a_i "wins").
person_top_choice_rawA <- apply(a_mat, 1, top_choice, decreasing = FALSE)
vote_rawA <- sort(table(person_top_choice_rawA), decreasing = TRUE)

na_counts <- colSums(is.na(a_mat))

# ── Report ─────────────────────────────────────────────────────────────
cat("=== GLOBAL model selection (method x k) ===\n\n")
cat("Winner by SILHOUETTE majority vote:", winner_silhouette, "\n\n")
cat("Silhouette vote distribution (top 15):\n"); print(head(vote_silhouette, 15))
cat("\nMean silhouette per partition (top 10, cross-check):\n")
print(head(mean_sil_per_partition, 10))
cat("\n[Diagnostic] Raw-cohesion vote winner (EXPECTED high-k, biased):",
    names(vote_rawA)[1], "\n")
cat("[Diagnostic] Raw-cohesion vote (top 10):\n"); print(head(vote_rawA, 10))
if (any(na_counts > 0)) {
  cat("\n[Diagnostic] Singleton-cluster NA counts per partition:\n")
  print(na_counts[na_counts > 0])
}

# ── Save (consolidated at end) ─────────────────────────────────────────
saveRDS(sil_mat,                file.path(out_dir, "silhouette_matrix.rds"))
saveRDS(a_mat,                  file.path(out_dir, "raw_cohesion_matrix.rds"))
saveRDS(person_top_choice,      file.path(out_dir, "person_top_choice_global.rds"))
saveRDS(person_pref_ranking,    file.path(out_dir, "person_pref_ranking_global.rds"))
saveRDS(vote_silhouette,        file.path(out_dir, "vote_silhouette_global.rds"))
saveRDS(mean_sil_per_partition, file.path(out_dir, "mean_silhouette_per_partition.rds"))
