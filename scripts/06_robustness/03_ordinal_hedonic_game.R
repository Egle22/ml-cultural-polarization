# =============================================================================
# scripts/06_robustness/03_ordinal_hedonic_game.R
# Thesis  : Appendix A.5 (extension with ordinal preferences); not a thesis table
# Input   : results/objects/clusters/{gower_dist,all_clusters_aligned}.rds
# Output  : results/objects/hedonic_ordinal/*.rds
# =============================================================================
#  Hedonic coalition formation: ORDINAL preference profiles
#  Karakaya (2011), Math. Soc. Sci. 61, 157-165.
#
#  The candidate partitions are read as hedonic games G = (N, >=), where
#  every respondent holds a RANKED LIST over the coalitions she could
#  belong to. The list is CONSTRUCTED from the Gower distances:
#
#     d_G(i,j)  ->  v_i(j) = tau - d_G(i,j)  ->  u_i(C) = mean_j v_i(j)
#                                             ->  sort  ->  the list
#
#  so that e.g.   C2 >_i C1 >_i {i}   whenever  u_i(C2) > u_i(C1) > 0.
#  Karakaya's definitions use only the ORDER; the numbers are the device
#  that produces it.
#
#  WHY tau. With the script-05 choice v_i(j) = 1 - d_G(i,j) >= 0 we get
#  u_i(C) >= 0 = u_i({i}) for every C: nobody can ever strictly prefer to
#  stand alone, individual rationality (Def. 3) holds vacuously, and the
#  singleton is structurally last in every list. Allowing v < 0 for pairs
#  further apart than tau lets a badly-matched coalition fall below the
#  singleton, which is what puts {i} in mid-list as in Example 1.
#  tau = 1 recovers script 05 exactly.
#
#  Utilities are FRACTIONAL (mean, not sum), as in Aziz et al.:
#     stay:    u_i(C_i)     = sum_{j in C_i} v_i(j) / |C_i|
#     join C:  u_i(C u {i}) = sum_{j in C}   v_i(j) / (|C| + 1)
#     alone:   u_i({i})     = 0
#
#  CAVEAT for the write-up: the division by |C| makes the game NOT
#  additively separable, so Bogomolnaia & Jackson's (2002) existence
#  result for additive + symmetric games does not apply. Non-existence of
#  a Nash stable partition is not ruled out a priori.
#
#  SELECTION RULE: per-k. For each k the methods are ranked by the share
#  of blocking agents and the least-contested method is retained. No
#  global (method, k) winner is declared: deviator shares are only
#  commensurable at a fixed k.
# ═══════════════════════════════════════════════════════════════════════

#### 0. Paths and printing helpers ######################################

source(here::here("R", "00_setup.R"))

clusters_dir <- DIR_CLUST
out_dir      <- ensure_dir(DIR_OBJ, "hedonic_ordinal")

hr <- function(title = NULL) {
  cat("\n", strrep("=", 72), "\n", sep = "")
  if (!is.null(title)) cat(title, "\n", strrep("=", 72), "\n", sep = "")
}
sub_hr <- function(title) cat("\n--- ", title, " ",
                              strrep("-", max(0, 66 - nchar(title))), "\n", sep = "")

#### 1. Load, with stale-object guard ###################################

rm(list = intersect(ls(), c("gower_dist", "all_clusters_aligned", "gower_mat",
                            "V", "tau", "pref_profiles", "hedonic_summary",
                            "exchange_res", "per_k", "winners", "tau_sens")))

hr("1. DATA")

gower_dist           <- readRDS(file.path(clusters_dir, "gower_dist.rds"))
all_clusters_aligned <- readRDS(file.path(clusters_dir, "all_clusters_aligned.rds"))

cat("Partitions available :", length(all_clusters_aligned), "\n")

## keep the retained families, k <= 6, plus HDBSCAN
nm          <- names(all_clusters_aligned)
keep_prefix <- grepl("^(km_umap_|pam_|km_pc2_k|km_pc5_k)", nm) |
  nm == "hdb_umap2_mp25"
k_num       <- suppressWarnings(as.integer(sub("^.*?([0-9]+)$", "\\1", nm)))
k_num[nm == "hdb_umap2_mp25"] <- NA_integer_
all_clusters_aligned <- all_clusters_aligned[keep_prefix &
                                               (is.na(k_num) | k_num <= 6)]

cat("Partitions retained  :", length(all_clusters_aligned), "\n")
cat("Dropped              :",
    paste(setdiff(nm, names(all_clusters_aligned)), collapse = ", "), "\n")

stopifnot(
  "gower_dist must be a 'dissimilarity'/'dist' object" =
    inherits(gower_dist, c("dissimilarity", "dist")),
  "at least two candidate partitions required" =
    length(all_clusters_aligned) >= 2L
)

gower_mat <- as.matrix(gower_dist)
n         <- nrow(gower_mat)

stopifnot(
  "partition length must match gower_dist size" =
    all(vapply(all_clusters_aligned, length, integer(1)) == n)
)

cluster_cols <- names(all_clusters_aligned)
is_hdb       <- grepl("^hdb", cluster_cols)

cat("Respondents          :", n, "\n")
cat("Pairwise distances   :", format(n * (n - 1) / 2, big.mark = ","), "\n\n")
cat("Gower distance distribution:\n")
print(round(summary(as.numeric(gower_dist)), 4))

sub_hr("Cluster sizes per partition")
for (col in cluster_cols) {
  tb <- table(all_clusters_aligned[[col]])
  cat(sprintf("  %-16s  %s\n", col,
              paste(sprintf("%s:%d", names(tb), as.integer(tb)), collapse = "  ")))
}

#### 2. Value function v_i(j) = tau - d_G(i,j) ##########################

hr("2. VALUE FUNCTION")

tau <- as.numeric(median(gower_dist))

V <- tau - gower_mat
diag(V) <- 0                                   # v_i(i) = 0 by convention

stopifnot(
  "V must be symmetric" = isTRUE(all.equal(V, t(V), check.attributes = FALSE)),
  "diagonal of V must be zero" = all(diag(V) == 0)
)

v_off <- V[upper.tri(V)]
cat(sprintf("tau (median Gower)        : %.4f\n", tau))
cat(sprintf("v_i(j) range              : [%.4f, %.4f]\n", min(v_off), max(v_off)))
cat(sprintf("'friends'  (v > 0)        : %6.2f%%\n", 100 * mean(v_off > 0)))
cat(sprintf("'enemies'  (v < 0)        : %6.2f%%\n", 100 * mean(v_off < 0)))
cat("\nWith tau = 1 every v would be non-negative and the singleton could\n",
    "never be strictly preferred to a coalition (script 05 behaviour).\n", sep = "")

#### 3. Ordinal preference profile for one partition ####################

## Returns each respondent's ranked list over the k+1 feasible coalitions
## (stay, join one of the others, stand alone) plus the diagnostics the
## ranking implies.
## noise_label: label treated as SINGLETONS rather than as a coalition
##              (0L for HDBSCAN, NULL for k-means / PAM).
hedonic_profile <- function(labels, V, noise_label = NULL, eps = 1e-9) {
  
  nU       <- length(labels)
  is_noise <- if (is.null(noise_label)) rep(FALSE, nU) else labels == noise_label
  cl       <- sort(unique(labels[!is_noise]))
  K        <- length(cl)
  stopifnot("fewer than two real coalitions" = K >= 2L)
  
  sizes <- vapply(cl, function(c) sum(labels == c & !is_noise), integer(1))
  
  ## S[i, c] = sum_{j in coalition c} v_i(j)   (self-term is zero)
  S <- vapply(cl, function(c) rowSums(V[, labels == c & !is_noise, drop = FALSE]),
              numeric(nU))
  dim(S) <- c(nU, K)
  
  own_col <- match(labels, cl)                 # NA for noise players
  
  ## Utility of every feasible alternative.
  ## Columns 1..K are the coalitions, column K+1 is the singleton.
  U <- sweep(S, 2, sizes + 1L, "/")            # joining: |C| + 1 denominator
  own_ok <- which(!is_noise)
  U[cbind(own_ok, own_col[own_ok])] <-         # staying: |C| denominator
    S[cbind(own_ok, own_col[own_ok])] / sizes[own_col[own_ok]]
  U <- cbind(U, 0)                             # u_i({i}) = 0
  colnames(U) <- c(paste0("C", cl), "alone")
  
  u_assigned <- numeric(nU)                    # noise players stand alone -> 0
  u_assigned[own_ok] <- U[cbind(own_ok, own_col[own_ok])]
  assigned_col <- ifelse(is_noise, K + 1L, own_col)
  
  ## Ranked list: order of the k+1 alternatives, best first.
  ord <- t(apply(-U, 1, order))
  rnk <- t(apply(-U, 1, rank, ties.method = "min"))
  rank_assigned <- rnk[cbind(seq_len(nU), assigned_col)]
  
  ## --- Def. 3: individual rationality --------------------------------
  ind_rational  <- u_assigned >= -eps
  prefers_alone <- u_assigned < -eps
  
  ## --- Def. 5: Nash blocking -----------------------------------------
  U_alt <- U
  U_alt[cbind(seq_len(nU), assigned_col)] <- -Inf
  best_alt <- apply(U_alt, 1, max)
  gain     <- pmax(best_alt - u_assigned, 0)
  blocker  <- gain > eps
  
  ## where the blockers would move
  best_col <- max.col(U_alt, ties.method = "first")
  move_to  <- ifelse(blocker,
                     ifelse(best_col == K + 1L, "alone",
                            paste0("C", cl[pmin(best_col, K)])),
                     NA_character_)
  
  ## --- Def. 9: top-choice property -----------------------------------
  top_choice <- rank_assigned == 1L
  
  list(
    K = K, cl = cl, sizes = sizes,
    n_noise = sum(is_noise), is_noise = is_noise,
    U = U, order = ord, rank = rnk,
    assigned_col = assigned_col, u_assigned = u_assigned,
    rank_assigned = rank_assigned, mean_rank = mean(rank_assigned),
    ind_rational = ind_rational, n_irrational = sum(!ind_rational),
    prefers_alone = prefers_alone, n_prefers_alone = sum(prefers_alone),
    blocker = blocker, gain = gain, move_to = move_to,
    n_blockers = sum(blocker),
    frac_blockers = mean(blocker),
    frac_blockers_core = if (any(!is_noise)) mean(blocker[!is_noise]) else NA_real_,
    frac_top = mean(top_choice), top_choice_prop = all(top_choice),
    nash_stable = sum(blocker) == 0L,
    total_gain = sum(gain), mean_gain = mean(gain)
  )
}

#### 4. Pairwise exchanges (Karakaya Def. 8, movement type iii) #########

## i in C_a and j in C_b swap coalitions; both entered coalitions keep
## their size, so
##      u_i^swap = ( S[i, b] - v_i(j) ) / |C_b|
##      u_j^swap = ( S[j, a] - v_i(j) ) / |C_a|
## The pair {i, j} strongly Nash blocks if BOTH strictly gain.
exchange_blockers <- function(labels, V, prof, eps = 1e-9) {
  
  cl <- prof$cl; K <- prof$K; sizes <- prof$sizes
  idx_of <- lapply(seq_len(K), function(a) which(labels == cl[a] & !prof$is_noise))
  
  S <- vapply(seq_len(K), function(c) rowSums(V[, idx_of[[c]], drop = FALSE]),
              numeric(length(labels)))
  dim(S) <- c(length(labels), K)
  
  n_pairs  <- 0
  involved <- logical(length(labels))
  
  for (a in seq_len(K - 1L)) for (b in (a + 1L):K) {
    ia <- idx_of[[a]]; jb <- idx_of[[b]]
    if (!length(ia) || !length(jb)) next
    
    Vij <- V[ia, jb, drop = FALSE]                       # |C_a| x |C_b|
    
    Gi <- sweep(-Vij, 1, S[ia, b], "+") / sizes[b]       # gain of i (rows)
    Gi <- sweep(Gi, 1, prof$u_assigned[ia], "-")
    
    Gj <- sweep(-Vij, 2, S[jb, a], "+") / sizes[a]       # gain of j (cols)
    Gj <- sweep(Gj, 2, prof$u_assigned[jb], "-")
    
    hit <- (Gi > eps) & (Gj > eps)
    if (any(hit)) {
      n_pairs      <- n_pairs + sum(hit)
      involved[ia] <- involved[ia] | (rowSums(hit) > 0L)
      involved[jb] <- involved[jb] | (colSums(hit) > 0L)
    }
    rm(Vij, Gi, Gj, hit)
  }
  
  list(n_exchange_pairs = n_pairs,
       n_involved       = sum(involved),
       frac_involved    = mean(involved),
       involved         = involved,
       strongly_nash    = n_pairs == 0 && prof$nash_stable)
}

#### 5. Preference-list renderer (Karakaya Example 1 style) #############

format_pref_list <- function(prof, i, digits = 3, latex = FALSE) {
  lab <- c(paste0("C", prof$cl), "{i}")
  o   <- prof$order[i, ]
  u   <- prof$U[i, o]
  sep <- if (latex) sprintf(" \\succ_{%d} ", i) else sprintf(" >_%d ", i)
  paste0(paste(lab[o], collapse = sep),
         "   [", paste(formatC(u, format = "f", digits = digits),
                       collapse = " > "), "]")
}

## One satisfied respondent, one Nash blocker, one who would rather stand
## alone, one noise point (when present).
illustrate <- function(prof, label) {
  sub_hr(paste("Preference lists:", label))
  picks <- c(
    `top choice` = which(prof$rank_assigned == 1L)[1],
    `blocker`    = which(prof$blocker)[1],
    `loner`      = which(prof$prefers_alone)[1],
    `noise`      = if (any(prof$is_noise)) which(prof$is_noise)[1] else NA_integer_
  )
  picks <- picks[!is.na(picks)]
  if (!length(picks)) { cat("  (no illustrative cases)\n"); return(invisible(NULL)) }
  for (nmk in names(picks)) {
    i <- picks[[nmk]]
    cat(sprintf("  %-11s resp. %5d:  %s\n", nmk, i, format_pref_list(prof, i)))
  }
  invisible(picks)
}

#### 6. Run over all candidate partitions ###############################

hr("3. PER-PARTITION DIAGNOSTICS")

cat(sprintf("%-16s %3s %6s %8s %7s %7s %7s %9s\n",
            "partition", "k", "noise", "blockers", "blk%", "meanRk", "top%", "swaps"))
cat(strrep("-", 72), "\n")

pref_profiles <- setNames(vector("list", length(cluster_cols)), cluster_cols)
exchange_res  <- setNames(vector("list", length(cluster_cols)), cluster_cols)

for (i in seq_along(cluster_cols)) {
  col <- cluster_cols[i]
  lab <- all_clusters_aligned[[col]]
  nl  <- if (is_hdb[i]) 0L else NULL
  
  p <- hedonic_profile(lab, V, noise_label = nl)
  e <- exchange_blockers(lab, V, p)
  
  pref_profiles[[col]] <- p
  exchange_res[[col]]  <- e
  
  cat(sprintf("%-16s %3d %6d %8d %6.2f%% %7.3f %6.1f%% %9s\n",
              col, p$K, p$n_noise, p$n_blockers, 100 * p$frac_blockers,
              p$mean_rank, 100 * p$frac_top,
              format(e$n_exchange_pairs, big.mark = ",")))
}

hedonic_summary <- data.frame(
  partition          = cluster_cols,
  method             = ifelse(is_hdb, "hdbscan", sub("_k?[0-9]+$", "", cluster_cols)),
  k                  = vapply(pref_profiles, function(x) x$K, integer(1)),
  n_noise            = vapply(pref_profiles, function(x) x$n_noise, integer(1)),
  n_blockers         = vapply(pref_profiles, function(x) x$n_blockers, integer(1)),
  frac_blockers      = vapply(pref_profiles, function(x) x$frac_blockers, numeric(1)),
  frac_blockers_core = vapply(pref_profiles, function(x) x$frac_blockers_core, numeric(1)),
  total_gain         = vapply(pref_profiles, function(x) x$total_gain, numeric(1)),
  n_prefers_alone    = vapply(pref_profiles, function(x) x$n_prefers_alone, integer(1)),
  mean_rank          = vapply(pref_profiles, function(x) x$mean_rank, numeric(1)),
  frac_top_choice    = vapply(pref_profiles, function(x) x$frac_top, numeric(1)),
  nash_stable        = vapply(pref_profiles, function(x) x$nash_stable, logical(1)),
  n_exchange_pairs   = vapply(exchange_res, function(x) x$n_exchange_pairs, numeric(1)),
  frac_in_exchange   = vapply(exchange_res, function(x) x$frac_involved, numeric(1)),
  strongly_nash      = vapply(exchange_res, function(x) x$strongly_nash, logical(1)),
  row.names          = NULL,
  stringsAsFactors   = FALSE
)

sub_hr("Individual rationality (Def. 3)")
cat("Respondents strictly preferring to stand alone:\n")
print(hedonic_summary[order(-hedonic_summary$n_prefers_alone),
                      c("partition", "k", "n_prefers_alone")],
      row.names = FALSE)
cat("\nUnder tau = 1 this column would be identically zero.\n")

#### 7. SELECTION: fix k, choose the method #############################

hr("4. SELECTION AT FIXED k")

cat("Deviator shares are comparable only within a k. For each k the\n",
    "methods are ranked by the share of blocking agents; ties are broken\n",
    "by the aggregate deviation gain.\n", sep = "")

k_levels <- sort(unique(hedonic_summary$k))
winners  <- vector("list", length(k_levels))

for (kk in seq_along(k_levels)) {
  k_val <- k_levels[kk]
  d <- hedonic_summary[hedonic_summary$k == k_val, ]
  d <- d[order(d$frac_blockers, d$total_gain), ]
  
  sub_hr(sprintf("k = %d   (%d candidate method%s)",
                 k_val, nrow(d), if (nrow(d) == 1) "" else "s"))
  
  print(data.frame(
    method       = d$method,
    partition    = d$partition,
    blockers     = d$n_blockers,
    `blockers_%` = round(100 * d$frac_blockers, 2),
    agg_gain     = round(d$total_gain, 3),
    mean_rank    = round(d$mean_rank, 3),
    `top_%`      = round(100 * d$frac_top_choice, 1),
    alone        = d$n_prefers_alone,
    swaps        = d$n_exchange_pairs,
    check.names  = FALSE
  ), row.names = FALSE)
  
  if (any(d$n_noise > 0)) {
    cat("  NOTE: rows with noise (HDBSCAN) count unassigned respondents as\n",
        "        singletons, so their blocker share sits on a different\n",
        "        denominator. Core-only share: ",
        sprintf("%.2f%%\n", 100 * d$frac_blockers_core[d$n_noise > 0][1]), sep = "")
  }
  
  best <- d[1, ]
  cat(sprintf("\n  >>> least contested at k = %d : %s (%s), %d blockers (%.2f%%)\n",
              k_val, best$method, best$partition, best$n_blockers,
              100 * best$frac_blockers))
  
  winners[[kk]] <- data.frame(
    k = k_val, method = best$method, partition = best$partition,
    n_blockers = best$n_blockers, frac_blockers = best$frac_blockers,
    total_gain = best$total_gain, mean_rank = best$mean_rank,
    frac_top_choice = best$frac_top_choice,
    n_prefers_alone = best$n_prefers_alone,
    n_exchange_pairs = best$n_exchange_pairs,
    nash_stable = best$nash_stable, strongly_nash = best$strongly_nash,
    row.names = NULL
  )
  
  illustrate(pref_profiles[[best$partition]], best$partition)
}

per_k <- do.call(rbind, winners)

sub_hr("Least-contested method at each k")
print(data.frame(
  k            = per_k$k,
  method       = per_k$method,
  blockers     = per_k$n_blockers,
  `blockers_%` = round(100 * per_k$frac_blockers, 2),
  mean_rank    = round(per_k$mean_rank, 3),
  `top_%`      = round(100 * per_k$frac_top_choice, 1),
  alone        = per_k$n_prefers_alone,
  nash         = per_k$nash_stable,
  strong_nash  = per_k$strongly_nash,
  check.names  = FALSE
), row.names = FALSE)

if (!any(hedonic_summary$nash_stable)) {
  cat("\nNo partition is Nash stable (Def. 5); a fortiori none is strongly\n",
      "Nash stable (Def. 8). Strong Nash stability is the stronger notion,\n",
      "so its failure adds little on its own: report the exchange COUNTS\n",
      "as the extra margin, not the bare fact of failure.\n", sep = "")
}

#### 8. Sensitivity of the per-k selection to tau #######################

hr("5. SENSITIVITY TO tau")

tau_grid <- c(median = as.numeric(median(gower_dist)),
              mean   = as.numeric(mean(gower_dist)),
              one    = 1)

tau_sens <- do.call(rbind, lapply(names(tau_grid), function(tn) {
  Vt <- tau_grid[[tn]] - gower_mat; diag(Vt) <- 0
  do.call(rbind, lapply(seq_along(cluster_cols), function(i) {
    p <- hedonic_profile(all_clusters_aligned[[cluster_cols[i]]], Vt,
                         noise_label = if (is_hdb[i]) 0L else NULL)
    data.frame(tau_name = tn, tau = tau_grid[[tn]],
               partition = cluster_cols[i],
               method = ifelse(is_hdb[i], "hdbscan",
                               sub("_k?[0-9]+$", "", cluster_cols[i])),
               k = p$K, frac_blockers = p$frac_blockers,
               total_gain = p$total_gain, row.names = NULL)
  }))
}))

sub_hr("Per-k winner under each tau")
tau_winners <- do.call(rbind, lapply(
  split(tau_sens, list(tau_sens$tau_name, tau_sens$k), drop = TRUE),
  function(d) {
    d <- d[order(d$frac_blockers, d$total_gain), ]
    data.frame(tau_name = d$tau_name[1], k = d$k[1], winner = d$method[1],
               row.names = NULL)
  }))
print(reshape(tau_winners, idvar = "k", timevar = "tau_name", direction = "wide"),
      row.names = FALSE)

cat("\nIf the winner is stable across the three thresholds, the selection\n",
    "does not rest on the choice of tau.\n", sep = "")

#### 9. Save ############################################################

hr("6. SAVED OBJECTS")

saveRDS(tau,             file.path(out_dir, "tau.rds"))
saveRDS(V,               file.path(out_dir, "value_matrix.rds"))
saveRDS(pref_profiles,   file.path(out_dir, "hedonic_pref_profiles.rds"))
saveRDS(exchange_res,    file.path(out_dir, "hedonic_exchanges.rds"))
saveRDS(hedonic_summary, file.path(out_dir, "hedonic_summary.rds"))
saveRDS(per_k,           file.path(out_dir, "hedonic_per_k.rds"))
saveRDS(tau_sens,        file.path(out_dir, "hedonic_tau_sensitivity.rds"))

cat("Written to", out_dir, ":\n")
cat(paste0("  ", list.files(out_dir), collapse = "\n"), "\n")