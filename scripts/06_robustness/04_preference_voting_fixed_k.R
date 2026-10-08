# =============================================================================
# scripts/06_robustness/04_preference_voting_fixed_k.R
# Thesis  : Section 4.1.4 (robustness); not a thesis table
# Input   : results/objects/preference_voting/silhouette_matrix.rds
#           (from 03_structure/05_silhouette.R)
# Output  : results/objects/preference_voting/{fixed_k_results,
#           fixed_k_winners_summary}.rds
# =============================================================================
#  Cluster model selection via a hedonic-game preference analogy
#  FIXED-k variant: for EACH number of clusters k, put the competing
#  methods head-to-head and choose the preferred METHOD.
#
#  This is the fairest method comparison available: at a fixed k the
#  silhouette scores of all methods live on the same [-1,1] scale, are
#  computed from the same Gower dissimilarity, and carry no residual
#  k-bias (k is held constant), so "which method places respondent i
#  best?" is well posed.
#
#  Preference = silhouette width s_i (same definition as scripts 01-02).
#  Reuses silhouette_matrix.rds from 01_preference_voting_global.R.
#  Column names are PARSED (not hardcoded), so any set of methods / any
#  k grid present in the matrix is handled automatically.
# ═══════════════════════════════════════════════════════════════════════

# ── Setup ──────────────────────────────────────────────────────────────
source(here::here("R", "00_setup.R"))

out_dir <- file.path(DIR_OBJ, "preference_voting")

rm(list = intersect(c("sil_mat"), ls()))
sil_mat <- readRDS(file.path(out_dir, "silhouette_matrix.rds"))
cols    <- colnames(sil_mat)

stopifnot("silhouette_matrix.rds must load with rows = respondents" = nrow(sil_mat) > 0)

# ── Parse (method, k) from each column name ─────────────────────────────
# k       = trailing run of digits (e.g. km_pc5_k4 -> 4, pam_16 -> 16)
# method  = everything before the trailing "_k?<digits>" (km_pc5_k4 -> km_pc5,
#           pam_16 -> pam, km_umap_4 -> km_umap)
k_of      <- function(x) as.integer(regmatches(x, regexpr("[0-9]+$", x)))
method_of <- function(x) sub("_k?[0-9]+$", "", x)

col_k <- vapply(cols, k_of,      integer(1))
col_m <- vapply(cols, method_of, character(1))

stopifnot(
  "failed to parse a k from some column name" = !any(is.na(col_k)),
  "parsing left an empty method for some column" = all(nzchar(col_m))
)

ks <- sort(unique(col_k))

top_choice <- function(row) {
  if (all(is.na(row))) return(NA_character_)
  names(row)[which.max(row)]   # higher silhouette = preferred
}

# ── For each fixed k: methods head-to-head ─────────────────────────────
per_k_results <- lapply(ks, function(k) {
  cols_k <- cols[col_k == k]
  if (length(cols_k) < 2L) return(NULL)  # need >=2 methods to compare

  sub <- sil_mat[, cols_k, drop = FALSE]
  colnames(sub) <- col_m[col_k == k]     # label columns by METHOD

  # per-respondent preferred method at this k
  person_pref_method <- apply(sub, 1, top_choice)
  vote      <- sort(table(person_pref_method), decreasing = TRUE)
  mean_sil  <- sort(colMeans(sub, na.rm = TRUE), decreasing = TRUE)

  list(
    k           = k,
    n_methods   = ncol(sub),
    vote        = vote,
    winner_vote = names(vote)[1],
    mean_sil    = mean_sil,
    winner_mean = names(mean_sil)[1],
    agree       = identical(names(vote)[1], names(mean_sil)[1])
  )
})
per_k_results <- Filter(Negate(is.null), per_k_results)
names(per_k_results) <- vapply(per_k_results, function(x) paste0("k", x$k), character(1))

# ── Summary table: winning method per k, by vote and by mean silhouette ─
fixed_k_winners <- data.frame(
  k                 = vapply(per_k_results, function(x) x$k, integer(1)),
  n_methods         = vapply(per_k_results, function(x) x$n_methods, integer(1)),
  winner_by_vote    = vapply(per_k_results, function(x) x$winner_vote, character(1)),
  votes             = vapply(per_k_results, function(x) as.integer(x$vote[1]), integer(1)),
  winner_by_meansil = vapply(per_k_results, function(x) x$winner_mean, character(1)),
  top_mean_sil      = vapply(per_k_results, function(x) unname(x$mean_sil[1]), numeric(1)),
  vote_meansil_agree= vapply(per_k_results, function(x) x$agree, logical(1)),
  row.names         = NULL,
  stringsAsFactors  = FALSE
)

# ── Report ─────────────────────────────────────────────────────────────
cat("=== FIXED-k model selection (choose the METHOD at each k) ===\n\n")
print(fixed_k_winners, row.names = FALSE)

if (any(!fixed_k_winners$vote_meansil_agree)) {
  cat("\nNOTE: at some k the plurality winner differs from the mean-silhouette",
      "winner.\nThese k values are where the methods are near-equivalent; do not",
      "rely on\nthe vote alone there. Cross-check with Phi_ST / Brier OOF.\n")
}

cat("\nFull breakdown per k (votes = # respondents best placed by that method;\n",
    "mean_sil = average silhouette of that method's partition at this k):\n", sep = "")
for (nm in names(per_k_results)) {
  r <- per_k_results[[nm]]
  cat(sprintf("\n-- k = %d --\n", r$k))
  cat("  votes:\n"); print(r$vote)
  cat("  mean silhouette:\n"); print(round(r$mean_sil, 5))
}

# ── Save (consolidated at end) ─────────────────────────────────────────
saveRDS(per_k_results,    file.path(out_dir, "fixed_k_results.rds"))
saveRDS(fixed_k_winners,  file.path(out_dir, "fixed_k_winners_summary.rds"))

