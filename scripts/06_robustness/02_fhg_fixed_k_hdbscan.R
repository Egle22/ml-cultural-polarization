# =============================================================================
# scripts/06_robustness/02_fhg_fixed_k_hdbscan.R
# Thesis  : Appendix A.5; Table A.10 (HDBSCAN row)
# Input   : results/objects/clusters/{gower_dist,all_clusters_aligned}.rds
#           results/objects/pca_scores.rds, results/objects/UMAP/umap_2d.rds
# Output  : results/objects/preference_voting/fhg_*.rds; deviator plots
# =============================================================================
#  FRACTIONAL HEDONIC GAME — FIXED-k METHOD COMPARISON
#  with explicit handling of the HDBSCAN noise group.
#
#  Three questions answered in one script:
#   (Q1) At each fixed k, which METHOD is closest to a Nash equilibrium
#        of the induced fractional hedonic game?
#   (Q2) HDBSCAN on the 2D UMAP embedding produces 3 substantive clusters
#        plus a residual noise group coded 0. Treating the noise as a
#        coalition inflates it to k = 4 and is game-theoretically wrong:
#        noise points are not a coalition, they are unassigned players.
#        We therefore drop them and evaluate HDBSCAN as a k = 3 partition
#        defined only on the ASSIGNED respondents.
#   (Q3) Where do the deviators sit in PCA score space and in UMAP space?
#
#  IMPORTANT — comparability warning:
#  Restricting HDBSCAN to its assigned points removes precisely the
#  respondents that sit in low-density regions, i.e. the ones most likely
#  to deviate. Comparing that against k-means/PAM computed on the full
#  n = 2579 is biased in HDBSCAN's favour. The script therefore reports
#  BOTH:
#    (A) full-sample deviator fractions for the complete-partition methods;
#    (B) the SAME comparison with every method restricted to the HDBSCAN
#        assigned subset — this is the fair like-for-like comparison.
#  Only (B) supports a sentence of the form "HDBSCAN has fewer deviators
#  than PAM at k = 3".
#
#  Utility (non-negative FHG, Aziz-Brandt-Harrenstein):
#     u_i(C) = ( sum_{j in C, j != i} v_i(j) ) / |C| ,  v_i(j) = 1 - d_G(i,j)
#     singleton payoff u_i({i}) = 0
#  i Nash-blocks the partition iff it can strictly gain by unilaterally
#  joining another existing coalition or leaving to a singleton.
# ═══════════════════════════════════════════════════════════════════════

source(here::here("R", "00_setup.R"))

# ── Paths ──────────────────────────────────────────────────────────────
obj_dir      <- DIR_OBJ
clusters_dir <- DIR_CLUST
out_dir      <- ensure_dir(obj_dir, "preference_voting")

# ADAPT HERE if the file names differ. Leave as NA to let the script
# search for a single plausible .rds under the corresponding folder.
pca_scores_path <- file.path(obj_dir, "pca_scores.rds")
umap2_path      <- file.path(obj_dir, "UMAP", "umap_2d.rds")

# ── Housekeeping: guard against stale interactive objects ──────────────
rm(list = intersect(
  c("gower_dist", "gower_mat", "clusters", "pca_scores", "umap2",
    "profiles_full", "profiles_rest", "fhg_full", "fhg_rest",
    "plot_df", "p_pca", "p_umap"),
  ls()))

#### 1. LOAD ############################################################

gower_dist <- readRDS(file.path(clusters_dir, "gower_dist.rds"))
clusters   <- readRDS(file.path(clusters_dir, "all_clusters_aligned.rds"))

stopifnot(
  "gower_dist must be a 'dissimilarity'/'dist' object" =
    inherits(gower_dist, c("dissimilarity", "dist")),
  "gower_dist size must match nrow(clusters)" =
    attr(gower_dist, "Size") == nrow(clusters)
)

gower_mat <- as.matrix(gower_dist)      # 2579 x 2579, ~53 MB
n         <- nrow(clusters)

# Resolve an embedding path: use the explicit one if given, otherwise
# search for exactly one candidate and fail loudly if ambiguous.
resolve_path <- function(explicit, search_dir, pattern) {
  if (!is.na(explicit) && file.exists(explicit)) return(explicit)
  cand <- list.files(search_dir, pattern = "\\.rds$", recursive = TRUE,
                     full.names = TRUE, ignore.case = TRUE)
  cand <- cand[grepl(pattern, basename(cand), ignore.case = TRUE)]
  if (length(cand) == 1L) return(cand)
  stop("Could not resolve the path for pattern '", pattern,
       "'. Candidates found:\n  ",
       paste(cand, collapse = "\n  "),
       "\nSet the path explicitly at the top of the script.")
}

# Coerce prcomp / FactoMineR::PCA / umap / uwot / matrix / data.frame
# objects to a plain 2-column numeric matrix of coordinates.
as_coords <- function(x, want = 2L) {
  m <- if (inherits(x, "prcomp"))          x$x
  else if (inherits(x, "PCA"))             x$ind$coord
  else if (is.list(x) && !is.null(x$layout)) x$layout   # umap::umap()
  else if (is.data.frame(x))               as.matrix(x)
  else                                     as.matrix(x)
  stopifnot("embedding has fewer columns than requested" = ncol(m) >= want)
  m[, seq_len(want), drop = FALSE]
}

pca_scores <- as_coords(readRDS(resolve_path(pca_scores_path,
                                             obj_dir, "pca|score")), 2L)
umap2      <- as_coords(readRDS(resolve_path(umap2_path,
                                             file.path(obj_dir, "UMAP"),
                                             "umap.*2|2.*umap")), 2L)

stopifnot(
  "PCA scores must have one row per respondent"  = nrow(pca_scores) == n,
  "UMAP embedding must have one row per respondent" = nrow(umap2)   == n
)

#### 2. PARSE (method, k) ###############################################

# Labels can arrive as factor/character/integer depending on how the
# Hungarian alignment stored them; normalise to integer once.
as_labels <- function(x) {
  if (is.factor(x)) return(as.integer(as.character(x)))
  as.integer(x)
}

cluster_cols <- names(clusters)
hdb_cols     <- grep("^hdb", cluster_cols, value = TRUE)
std_cols     <- setdiff(cluster_cols, hdb_cols)

# Standard columns: k = trailing digits, method = the rest.
# The optional "k" before the digits covers km_pc5_k4 vs pam_4.
k_of      <- function(x) as.integer(regmatches(x, regexpr("[0-9]+$", x)))
method_of <- function(x) sub("_k?[0-9]+$", "", x)

std_k <- vapply(std_cols, k_of,      integer(1))
std_m <- vapply(std_cols, method_of, character(1))
stopifnot("failed to parse k from a standard column" = !any(is.na(std_k)),
          "parsing left an empty method name"        = all(nzchar(std_m)))

# HDBSCAN: k is NOT in the name (mp25 = minPts, not k). Recover it from
# the labels themselves, counting only non-noise clusters.
hdb_effective_k <- vapply(hdb_cols, function(cc) {
  lab <- as_labels(clusters[[cc]])
  length(setdiff(unique(lab[!is.na(lab)]), 0L))
}, integer(1))

cat("HDBSCAN columns and their effective k (noise excluded):\n")
print(data.frame(column = hdb_cols, effective_k = hdb_effective_k,
                 row.names = NULL))

# Build the label list used everywhere below. NA = unassigned player.
label_list <- c(
  setNames(lapply(std_cols, function(cc) as_labels(clusters[[cc]])), std_cols),
  setNames(lapply(hdb_cols, function(cc) {
    lab <- as_labels(clusters[[cc]]); lab[lab == 0L] <- NA_integer_; lab
  }), hdb_cols)
)
meta <- data.frame(
  partition = c(std_cols, hdb_cols),
  method    = c(unname(std_m), sub("_mp[0-9]+$", "", hdb_cols)),
  k         = c(unname(std_k), unname(hdb_effective_k)),
  row.names = NULL, stringsAsFactors = FALSE
)

#### 3. FHG DEVIATION PROFILE ###########################################

# labels: integer vector of length n; NA marks players excluded from the
#         game (HDBSCAN noise, or respondents outside a restriction set).
# Utilities are computed on the SUBMATRIX of assigned players only, so the
# coalitions a deviator can join are exactly the coalitions that exist.
fhg_deviation_profile <- function(labels, gm, eps = 1e-9) {
  keep <- !is.na(labels)
  idx  <- which(keep)
  lab  <- labels[idx]
  gsub_mat <- gm[idx, idx, drop = FALSE]
  nU <- length(idx)

  cl    <- sort(unique(lab)); sizes <- as.integer(table(factor(lab, levels = cl)))
  names(sizes) <- as.character(cl)

  # dist_sum[i, c] = sum_{j in cluster c} d_G(i, j)   (self-distance = 0)
  dist_sum <- vapply(cl, function(c) rowSums(gsub_mat[, lab == c, drop = FALSE]),
                     numeric(nU))
  colnames(dist_sum) <- as.character(cl)

  own_idx   <- match(as.character(lab), colnames(dist_sum))
  own_sizes <- sizes[own_idx]
  own_dsum  <- dist_sum[cbind(seq_len(nU), own_idx)]

  # sum of similarities to co-members = (|own| - 1) - own_dsum, since sim = 1 - d
  u_curr <- ((own_sizes - 1) - own_dsum) / own_sizes        # singleton -> 0

  # utility of JOINING cluster c: (|c| - dist_sum) / (|c| + 1)
  U_join <- sweep(-dist_sum, 2, sizes, "+")
  U_join <- sweep(U_join,    2, sizes + 1, "/")
  U_join[cbind(seq_len(nU), own_idx)] <- -Inf               # cannot rejoin own

  best_alt <- pmax(apply(U_join, 1, max), 0)                # 0 = leave to singleton
  gain     <- pmax(best_alt - u_curr, 0)
  blocker  <- gain > eps

  # Re-expand to full length so results align with the original rows.
  blocker_full <- rep(NA, length(labels)); blocker_full[idx] <- blocker
  gain_full    <- rep(NA_real_, length(labels)); gain_full[idx] <- gain

  list(
    n_assigned    = nU,
    n_blockers    = sum(blocker),
    frac_blockers = mean(blocker),          # denominator = ASSIGNED players
    pct_blockers  = 100 * mean(blocker),
    total_gain    = sum(gain),
    mean_gain     = mean(gain),
    nash_stable   = sum(blocker) == 0L,
    blocker       = blocker_full,
    gain          = gain_full
  )
}

summarise_profiles <- function(prof_list, meta_df) {
  data.frame(
    partition     = meta_df$partition,
    method        = meta_df$method,
    k             = meta_df$k,
    n_assigned    = vapply(prof_list, function(x) x$n_assigned,    integer(1)),
    n_blockers    = vapply(prof_list, function(x) x$n_blockers,    integer(1)),
    pct_blockers  = vapply(prof_list, function(x) x$pct_blockers,  numeric(1)),
    total_gain    = vapply(prof_list, function(x) x$total_gain,    numeric(1)),
    mean_gain     = vapply(prof_list, function(x) x$mean_gain,     numeric(1)),
    nash_stable   = vapply(prof_list, function(x) x$nash_stable,   logical(1)),
    row.names = NULL, stringsAsFactors = FALSE
  )
}

#### 4. (A) FULL-SAMPLE ANALYSIS ########################################

profiles_full <- lapply(label_list, fhg_deviation_profile, gm = gower_mat)
fhg_full      <- summarise_profiles(profiles_full, meta)

cat("\n=== (A) FULL SAMPLE — deviator percentages ===\n")
cat("NOTE: HDBSCAN rows use only their assigned respondents as denominator,\n",
    "so they are NOT directly comparable with the other rows here.\n", sep = "")
print(fhg_full[order(fhg_full$k, fhg_full$pct_blockers),
               c("partition", "k", "n_assigned", "n_blockers",
                 "pct_blockers", "total_gain", "nash_stable")],
      row.names = FALSE, digits = 4)

# Winner at each fixed k: fewest deviators, ties broken on total_gain.
winner_by_k <- do.call(rbind, lapply(split(fhg_full, fhg_full$k), function(d) {
  d <- d[order(d$pct_blockers, d$total_gain), ]
  data.frame(k = d$k[1], n_methods = nrow(d),
             winner = d$partition[1], winner_pct = d$pct_blockers[1],
             runner_up = if (nrow(d) > 1) d$partition[2] else NA_character_,
             runner_up_pct = if (nrow(d) > 1) d$pct_blockers[2] else NA_real_,
             row.names = NULL)
}))

cat("\n--- (A) Preferred METHOD at each fixed k ---\n")
print(winner_by_k, row.names = FALSE, digits = 4)

#### 5. (B) FAIR COMPARISON ON THE HDBSCAN-ASSIGNED SUBSET ##############

hdb_col <- hdb_cols[hdb_effective_k == 3L]
stopifnot("expected exactly one HDBSCAN column with 3 non-noise clusters" =
            length(hdb_col) == 1L)
hdb_col <- hdb_col[1]

assigned_idx <- which(!is.na(label_list[[hdb_col]]))
cat(sprintf("\nHDBSCAN (%s): %d assigned / %d total (%.1f%% noise dropped)\n",
            hdb_col, length(assigned_idx), n,
            100 * (1 - length(assigned_idx) / n)))

# Restrict EVERY method to the same respondents, keeping its own labels.
restrict_to <- function(labels, idx) {
  out <- rep(NA_integer_, length(labels)); out[idx] <- labels[idx]; out
}
label_list_rest <- lapply(label_list, restrict_to, idx = assigned_idx)

profiles_rest <- lapply(label_list_rest, fhg_deviation_profile, gm = gower_mat)
fhg_rest      <- summarise_profiles(profiles_rest, meta)

k3_rest <- fhg_rest[fhg_rest$k == 3L, ]
k3_rest <- k3_rest[order(k3_rest$pct_blockers, k3_rest$total_gain), ]

cat("\n=== (B) FAIR COMPARISON AT k = 3 ===\n")
cat("All methods restricted to the ", length(assigned_idx),
    " respondents assigned by HDBSCAN.\n", sep = "")
cat("Percentages are computed over these assigned respondents only.\n\n")
print(k3_rest[, c("partition", "method", "n_assigned", "n_blockers",
                  "pct_blockers", "total_gain", "nash_stable")],
      row.names = FALSE, digits = 4)

cat(sprintf("\nHDBSCAN deviators at k = 3: %d / %d = %.2f%%\n",
            fhg_rest$n_blockers[fhg_rest$partition == hdb_col],
            fhg_rest$n_assigned[fhg_rest$partition == hdb_col],
            fhg_rest$pct_blockers[fhg_rest$partition == hdb_col]))
cat(sprintf("Best method at k = 3 on this subset: %s (%.2f%% deviators)\n",
            k3_rest$partition[1], k3_rest$pct_blockers[1]))

if (!any(k3_rest$nash_stable)) {
  cat("\nNo exact Nash-stable partition at k = 3: report the winner as\n",
      "'closest to equilibrium', never as 'an equilibrium'.\n", sep = "")
}

# Side-by-side: same method, full sample vs restricted subset. The gap is
# the size of the selection effect introduced by dropping the noise.
cmp <- merge(
  fhg_full[fhg_full$k == 3L, c("partition", "pct_blockers")],
  fhg_rest[fhg_rest$k == 3L, c("partition", "pct_blockers")],
  by = "partition", suffixes = c("_full", "_restricted")
)
cmp$delta <- cmp$pct_blockers_restricted - cmp$pct_blockers_full
cat("\n--- Selection effect of dropping HDBSCAN noise (k = 3) ---\n")
print(cmp[order(cmp$pct_blockers_restricted), ], row.names = FALSE, digits = 4)

#### 6. PLOTS — WHERE THE DEVIATORS SIT #################################

# "restricted" keeps every method on the same respondents (fair, matches
# table B); switch to "full" to see the deviators of the complete partitions.
plot_basis <- "restricted"          # "restricted" or "full"
prof_plot  <- if (plot_basis == "restricted") profiles_rest else profiles_full

k3_partitions <- meta$partition[meta$k == 3L]

plot_df <- do.call(rbind, lapply(k3_partitions, function(p) {
  b <- prof_plot[[p]]$blocker
  data.frame(
    partition = p,
    PC1 = pca_scores[, 1], PC2 = pca_scores[, 2],
    U1  = umap2[, 1],      U2  = umap2[, 2],
    status = factor(ifelse(is.na(b), "unassigned",
                           ifelse(b, "deviator", "stable")),
                    levels = c("stable", "deviator", "unassigned")),
    row.names = NULL, stringsAsFactors = FALSE
  )
}))

# Facet labels carry the deviator percentage, so the panels are readable
# without going back to the table.
pct_lookup <- setNames(
  if (plot_basis == "restricted") fhg_rest$pct_blockers else fhg_full$pct_blockers,
  meta$partition)
plot_df$facet <- sprintf("%s  (%.1f%% dev.)",
                         plot_df$partition, pct_lookup[plot_df$partition])

pal <- c(stable = "grey75", deviator = "#C0392B", unassigned = "grey92")

make_scatter <- function(df, xv, yv, xlab, ylab, title) {
  ggplot(df, aes(.data[[xv]], .data[[yv]], colour = status)) +
    # draw stable/unassigned first so deviators are not hidden underneath
    geom_point(data = subset(df, status != "deviator"), size = 0.5, alpha = 0.5) +
    geom_point(data = subset(df, status == "deviator"), size = 0.7, alpha = 0.9) +
    scale_colour_manual(values = pal, drop = FALSE) +
    facet_wrap(~ facet, ncol = 3) +
    labs(x = xlab, y = ylab, title = title, colour = NULL) +
    guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
    theme_minimal(base_size = 9) +
    theme(legend.position = "bottom",
          panel.grid.minor = element_blank(),
          strip.text = element_text(face = "bold", size = 8))
}

p_pca <- make_scatter(plot_df, "PC1", "PC2",
                      "PC1", "PC2",
                      sprintf("FHG deviators at k = 3 — PCA score space (%s basis)",
                              plot_basis))
p_umap <- make_scatter(plot_df, "U1", "U2",
                       "UMAP 1", "UMAP 2",
                       sprintf("FHG deviators at k = 3 — UMAP space (%s basis)",
                               plot_basis))

print(p_pca)
print(p_umap)

#### 7. SAVE (consolidated at end; plots deliberately NOT saved) ########

saveRDS(profiles_full, file.path(out_dir, "fhg_profiles_full.rds"))
saveRDS(profiles_rest, file.path(out_dir, "fhg_profiles_hdb_subset.rds"))
saveRDS(fhg_full,      file.path(out_dir, "fhg_summary_full.rds"))
saveRDS(fhg_rest,      file.path(out_dir, "fhg_summary_hdb_subset.rds"))
saveRDS(winner_by_k,   file.path(out_dir, "fhg_winner_by_k.rds"))
saveRDS(k3_rest,       file.path(out_dir, "fhg_k3_fair_comparison.rds"))
saveRDS(cmp,           file.path(out_dir, "fhg_k3_selection_effect.rds"))
