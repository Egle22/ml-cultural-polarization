# =============================================================================
# scripts/06_robustness/01_fhg_equilibrium.R
# Thesis  : Appendix A.5; Table A.10 (deviators under the fractional hedonic
#           game), Table 4.3
# Input   : results/objects/clusters/{gower_dist,all_clusters_aligned}.rds
# Output  : results/objects/preference_voting/{fhg_profiles,fhg_summary,
#           fhg_winner_global,fhg_per_method,fhg_per_k}.rds
# =============================================================================
#  Cluster model selection via FRACTIONAL HEDONIC GAME equilibrium
#  (Aziz, Brandt & Harrenstein). Unlike scripts 01-03, which rank
#  partitions by respondent "comfort" (silhouette), this script asks a
#  different question: which candidate partition is (closest to) a Nash
#  equilibrium of the induced fractional hedonic game?
#
#  Fractional hedonic utility of respondent i in coalition C (i in C):
#     u_i(C) = ( sum_{j in C} v_i(j) ) / |C| ,   v_i(i) = 0
#  with similarity values  v_i(j) = 1 - d_G(i, j)  in [0, 1]
#  (standard NON-NEGATIVE FHG: players prefer groups whose members are, on
#   average, similar to them; the singleton payoff is u_i({i}) = 0).
#
#  Respondent i "deviates" (Nash-blocks) partition pi if she can strictly
#  raise her utility by moving, on her own, either
#     - into another existing cluster C:  u_i(C u {i}) > u_i(pi(i)) , or
#     - out to a singleton:               0            > u_i(pi(i)) .
#  A partition is Nash-stable iff NO respondent deviates (n_blockers = 0).
#
#  Selection rule (all three granularities): choose the partition with the
#  FEWEST deviators; ties broken by the smallest aggregate deviation
#  incentive (total_gain). "Fewest deviators" = "closest to equilibrium".
#
#  Column names are PARSED (not hardcoded): any method set / k grid works.
#
#  NOTE (honesty): k-means / PAM partitions optimise a different objective,
#  so an EXACT equilibrium (0 deviators) will usually not exist. The script
#  reports whether one does; if not, it selects the closest and prints the
#  residual deviator fraction so the distance from equilibrium is explicit.
# ═══════════════════════════════════════════════════════════════════════

# ── Setup ──────────────────────────────────────────────────────────────
# (base R only; no clustering package needed here)
source(here::here("R", "00_setup.R"))

# ── Paths ──────────────────────────────────────────────────────────────
clusters_dir <- DIR_CLUST
out_dir      <- ensure_dir(DIR_OBJ, "preference_voting")

# ── Load data (rm() guard against stale interactive objects) ───────────
rm(list = intersect(c("gower_dist", "all_clusters_aligned", "gower_mat"), ls()))

gower_dist           <- readRDS(file.path(clusters_dir, "gower_dist.rds"))
all_clusters_aligned <- readRDS(file.path(clusters_dir, "all_clusters_aligned.rds"))
nm <- names(all_clusters_aligned)
keep_prefix <- grepl("^(km_umap_|pam_|km_pc2_k|km_pc5_k)", nm) | nm == "hdb_umap2_mp25"

## 2. extract trailing k (NA for hdb, which has no k)
k_num <- suppressWarnings(as.integer(sub("^.*?([0-9]+)$", "\\1", nm)))
k_num[nm == "hdb_umap2_mp25"] <- NA_integer_

## 3. keep k <= 6, plus hdbscan (no k)
keep_idx <- keep_prefix & (is.na(k_num) | k_num <= 6)

all_clusters_aligned <- all_clusters_aligned[keep_idx]
stopifnot(
  "gower_dist must be a 'dissimilarity'/'dist' object" =
    inherits(gower_dist, c("dissimilarity", "dist")),
  "gower_dist size must match nrow(all_clusters_aligned)" =
    attr(gower_dist, "Size") == nrow(all_clusters_aligned)
)

gower_mat    <- as.matrix(gower_dist)          # 2579 x 2579, ~53 MB
n            <- nrow(all_clusters_aligned)
cluster_cols <- names(all_clusters_aligned)

# ── Parse (method, k) from each column name ─────────────────────────────
k_of      <- function(x) as.integer(regmatches(x, regexpr("[0-9]+$", x)))
method_of <- function(x) sub("_k?[0-9]+$", "", x)

col_k <- vapply(cluster_cols, k_of,      integer(1))
col_m <- vapply(cluster_cols, method_of, character(1))
stopifnot(
  "failed to parse a k from some column name"       = !any(is.na(col_k)),
  "parsing left an empty method for some column"    = all(nzchar(col_m))
)

# ── FHG deviation profile for one partition ─────────────────────────────
# Returns per-partition equilibrium diagnostics plus per-respondent vectors.
fhg_deviation_profile <- function(labels, gm, eps = 1e-9) {
  cl    <- sort(unique(labels)); K <- length(cl); nU <- length(labels)
  sizes <- as.integer(table(factor(labels, levels = cl)))
  names(sizes) <- as.character(cl)

  # dist_sum[i, c] = sum_{j in cluster c} d_G(i, j)   (self-distance = 0)
  dist_sum <- vapply(cl, function(c) rowSums(gm[, labels == c, drop = FALSE]),
                     numeric(nU))
  colnames(dist_sum) <- as.character(cl)

  own_idx   <- match(as.character(labels), colnames(dist_sum))
  own_sizes <- sizes[own_idx]
  own_dsum  <- dist_sum[cbind(seq_len(nU), own_idx)]

  # current utility: sum of similarities to co-members / |own|
  #   sum_sim_to_others = (|own| - 1) - own_dsum   [since sim = 1 - d]
  u_curr <- ((own_sizes - 1) - own_dsum) / own_sizes           # singleton -> 0

  # utility of JOINING each cluster c: (|c| - dist_sum) / (|c| + 1)
  U_join <- sweep(-dist_sum, 2, sizes, "+")                    # |c| - dist_sum
  U_join <- sweep(U_join,   2, sizes + 1, "/")
  U_join[cbind(seq_len(nU), own_idx)] <- -Inf                 # cannot join own cluster

  best_alt <- pmax(apply(U_join, 1, max), 0)                  # 0 = singleton fallback
  gain     <- pmax(best_alt - u_curr, 0)                      # improvement available to i
  blocker  <- gain > eps

  list(
    n_blockers   = sum(blocker),
    frac_blockers= mean(blocker),
    total_gain   = sum(gain),
    mean_gain    = mean(gain),
    nash_stable  = sum(blocker) == 0L,
    blocker      = blocker,
    gain         = gain
  )
}

# ── Run FHG check on every candidate partition ─────────────────────────
profiles <- lapply(cluster_cols, function(col)
  fhg_deviation_profile(all_clusters_aligned[[col]], gower_mat))
names(profiles) <- cluster_cols

fhg_summary <- data.frame(
  partition     = cluster_cols,
  method        = col_m,
  k             = col_k,
  n_blockers    = vapply(profiles, function(x) x$n_blockers,   integer(1)),
  frac_blockers = vapply(profiles, function(x) x$frac_blockers, numeric(1)),
  total_gain    = vapply(profiles, function(x) x$total_gain,   numeric(1)),
  mean_gain     = vapply(profiles, function(x) x$mean_gain,    numeric(1)),
  nash_stable   = vapply(profiles, function(x) x$nash_stable,  logical(1)),
  row.names     = NULL,
  stringsAsFactors = FALSE
)

# ── Selection helper: fewest deviators, tie-break on total_gain ────────
pick_closest <- function(df) {
  o <- order(df$frac_blockers, df$total_gain)   # ascending on both
  df$partition[o[1]]
}

# (a) GLOBAL: choose method AND k jointly
winner_global <- pick_closest(fhg_summary)

# (b) PER-METHOD: fix method, choose k
per_method <- do.call(rbind, lapply(split(fhg_summary, fhg_summary$method), function(d) {
  best <- d[order(d$frac_blockers, d$total_gain), ][1, ]
  data.frame(method = best$method, chosen_partition = best$partition,
             chosen_k = best$k, frac_blockers = best$frac_blockers,
             nash_stable = best$nash_stable, row.names = NULL)
}))

# (c) PER-k: fix k, choose method
per_k <- do.call(rbind, lapply(split(fhg_summary, fhg_summary$k), function(d) {
  best <- d[order(d$frac_blockers, d$total_gain), ][1, ]
  data.frame(k = best$k, chosen_partition = best$partition,
             chosen_method = best$method, frac_blockers = best$frac_blockers,
             nash_stable = best$nash_stable, row.names = NULL)
}))

# ── Report ─────────────────────────────────────────────────────────────
cat("=== FRACTIONAL HEDONIC GAME: distance-from-equilibrium selection ===\n\n")

any_stable <- any(fhg_summary$nash_stable)
if (any_stable) {
  cat("Exact Nash-stable partition(s) found (0 deviators):\n")
  print(fhg_summary[fhg_summary$nash_stable, c("partition","frac_blockers","total_gain")],
        row.names = FALSE)
} else {
  cat("No EXACT equilibrium exists among the candidates (expected for\n",
      "k-means/PAM). Selecting the partition CLOSEST to equilibrium.\n", sep = "")
}

cat("\n(a) GLOBAL winner (method & k):", winner_global,
    sprintf("  [deviators: %.1f%%]\n",
            100 * fhg_summary$frac_blockers[fhg_summary$partition == winner_global]))

cat("\n(b) PER-METHOD (fix method, choose k):\n")
print(per_method, row.names = FALSE)

cat("\n(c) PER-k (fix k, choose method):\n")
print(per_k, row.names = FALSE)

cat("\nFull FHG summary (sorted by deviator fraction):\n")
print(fhg_summary[order(fhg_summary$frac_blockers, fhg_summary$total_gain),
                  c("partition","n_blockers","frac_blockers","total_gain","nash_stable")],
      row.names = FALSE)

best_frac <- min(fhg_summary$frac_blockers)
if (!any_stable && best_frac > 0.10) {
  cat(sprintf("\nWARNING: even the best candidate has %.1f%% deviators, i.e. it is\n",
              100 * best_frac),
      "far from a true equilibrium. Report the selection as 'closest to\n",
      "equilibrium', NOT as an equilibrium, and cross-read with the\n",
      "silhouette results (scripts 01-03) and Phi_ST.\n", sep = "")
}

# ── Save (consolidated at end) ─────────────────────────────────────────
saveRDS(profiles,       file.path(out_dir, "fhg_profiles.rds"))
saveRDS(fhg_summary,    file.path(out_dir, "fhg_summary.rds"))
saveRDS(winner_global,  file.path(out_dir, "fhg_winner_global.rds"))
saveRDS(per_method,     file.path(out_dir, "fhg_per_method.rds"))
saveRDS(per_k,          file.path(out_dir, "fhg_per_k.rds"))
