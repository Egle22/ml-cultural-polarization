# =============================================================================
# scripts/04_profiles/03_method_overlap.R
# Purpose : How similar are the TOP_N variables selected by the LASSO in the
#           four clustering methods, for k = K_CHOSEN only (run after
#           02_ovr_lasso.R).
# Thesis  : Section 4.2; Table A.7 (variables selected by >= 2 methods, k = 3)
# Input   : results/objects/logit_ovr_ord/logit_lasso_top20_selected_k<K>.rds
# Output  : results/objects/logit_ovr_ord/ (suffix k<K_CHOSEN>):
#   top_overlap_pairwise_k*.rds     pairs of methods, by cluster
#   top_overlap_votes_k*.rds        variable x number of methods choosing it
#   top_overlap_core_k*.rds         core set (all methods)
#   top_overlap_pooled_k*.rds       pairs of methods, clusters pooled
#           results/figures/04_profiles/top_overlap_jaccard_k*.png
#
# -- CHOICES -------------------------------------------------------------------
# [0] The comparison is BY CLUSTER: cluster c of method A against cluster c
#     of method B. This is meaningful because the labels are already aligned
#     (Hungarian) in all_clusters_aligned. A POOLED version (union of the
#     clusters) is also given at the end; it compares the methods as "sets of
#     relevant variables" without relying on the alignment.
#
# [1] Measures: n_shared (intersection) and Jaccard = |A and B| / |A or B|.
#     Baseline: two random sets of size n1 and n2 drawn from P variables are
#     expected to share n1*n2/P. The column exp_shared_random tells whether
#     the observed overlap is informative or what chance would give.
#
# [2] "Core set" = variables selected by ALL methods for that cluster.
#     vote_table reports the full count (1..4 methods), which is the most
#     useful form for the thesis: it shows a gradient, not a dichotomy.
# =============================================================================

# -- Loading -------------------------------------------------------------------
# The selection is always read from disk for the current K_CHOSEN, so that a
# stale object left in the session by a run with another k is never reused.
if (!exists("K_CHOSEN") || !exists("dat_model"))
  source(here::here("scripts", "04_profiles", "01_ovr_framework.R"))
top_lasso_vars <- readRDS(file.path(
  OUT_LOG, sprintf("logit_lasso_top20_selected_k%d.rds", K_CHOSEN)))

tag_k    <- sprintf("k%d", K_CHOSEN)
CLUSTERS <- sort(unique(as.character(top_lasso_vars$cluster)))
P_TOT    <- ncol(dat_model)          # number of candidate variables (chance denominator)

get_top <- function(m, cl) {
  unique(top_lasso_vars$variable[top_lasso_vars$method == m &
                                   as.character(top_lasso_vars$cluster) == cl])
}

jacc <- function(a, b) if (length(union(a, b)) == 0) NA_real_ else
  length(intersect(a, b)) / length(union(a, b))

# =============================================================================
#  (1) PAIRWISE INTERSECTIONS, cluster by cluster
# =============================================================================
mp <- t(utils::combn(METHODS, 2))

pair_overlap <- do.call(rbind, lapply(CLUSTERS, function(cl) {
  do.call(rbind, lapply(seq_len(nrow(mp)), function(i) {
    a <- get_top(mp[i, 1], cl); b <- get_top(mp[i, 2], cl)
    if (length(a) == 0 || length(b) == 0) return(NULL)
    sh <- sort(intersect(a, b))
    data.frame(
      k = K_CHOSEN, cluster = cl,
      method1 = mp[i, 1], method2 = mp[i, 2],
      n1 = length(a), n2 = length(b),
      n_shared = length(sh),
      jaccard  = round(jacc(a, b), 3),
      exp_shared_random = round(length(a) * length(b) / P_TOT, 2),
      shared = paste(sh, collapse = ", "),
      stringsAsFactors = FALSE)
  }))
}))

# =============================================================================
#  (2) INTERSECTION ACROSS ALL METHODS + count of the "votes"
# =============================================================================
vote_table <- top_lasso_vars |>
  mutate(cluster = as.character(cluster)) |>
  group_by(cluster, variable, section_label) |>
  summarise(n_methods = dplyr::n_distinct(method),
            methods   = paste(sort(unique(method)), collapse = " + "),
            mean_rank = round(mean(rank), 1),
            .groups = "drop") |>
  arrange(cluster, desc(n_methods), mean_rank)

core_table <- vote_table |>
  filter(n_methods == length(METHODS)) |>
  dplyr::select(cluster, variable, section_label, mean_rank)

# compact summary: how many variables are shared by 1, 2, 3, 4 methods
vote_summary <- vote_table |>
  dplyr::count(cluster, n_methods) |>
  tidyr::pivot_wider(names_from = n_methods, values_from = n,
                     names_prefix = "n_meth_", values_fill = 0L)

# =============================================================================
#  (3) POOLED VERSION: methods compared ignoring the clusters
# =============================================================================
pooled <- setNames(lapply(METHODS, function(m)
  unique(top_lasso_vars$variable[top_lasso_vars$method == m])), METHODS)

pooled_overlap <- do.call(rbind, lapply(seq_len(nrow(mp)), function(i) {
  a <- pooled[[mp[i, 1]]]; b <- pooled[[mp[i, 2]]]
  data.frame(k = K_CHOSEN, method1 = mp[i, 1], method2 = mp[i, 2],
             n1 = length(a), n2 = length(b),
             n_shared = length(intersect(a, b)),
             jaccard = round(jacc(a, b), 3),
             exp_shared_random = round(length(a) * length(b) / P_TOT, 2),
             stringsAsFactors = FALSE)
}))

pooled_core <- Reduce(intersect, pooled)

# -- Print ---------------------------------------------------------------------
cat(sprintf("\n=== Overlap top-%d, k = %d (%d candidate variables) ===\n",
            max(pair_overlap$n1), K_CHOSEN, P_TOT))

cat("\n-- Pairs of methods, by cluster --\n")
print(pair_overlap[, c("cluster", "method1", "method2", "n_shared",
                       "jaccard", "exp_shared_random")], row.names = FALSE)

cat("\n-- How many variables are shared by how many methods --\n")
print(as.data.frame(vote_summary), row.names = FALSE)

cat("\n-- Core set (all methods) --\n")
if (nrow(core_table) == 0) cat("  no variable common to the 4 methods\n") else
  print(as.data.frame(core_table), row.names = FALSE)

cat("\n-- Pooled (clusters pooled) --\n")
print(pooled_overlap, row.names = FALSE)
cat(sprintf("Pooled core (%d variables): %s\n",
            length(pooled_core), paste(sort(pooled_core), collapse = ", ")))

# -- Jaccard heatmap by cluster ------------------------------------------------
hm <- rbind(
  pair_overlap[, c("cluster", "method1", "method2", "n_shared", "jaccard")],
  setNames(pair_overlap[, c("cluster", "method2", "method1", "n_shared", "jaccard")],
           c("cluster", "method1", "method2", "n_shared", "jaccard")))
hm$method1 <- factor(hm$method1, levels = METHODS)
hm$method2 <- factor(hm$method2, levels = rev(METHODS))

p_hm <- ggplot(hm, aes(method1, method2, fill = jaccard)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%d\n(%.2f)", n_shared, jaccard)), size = 2.8) +
  scale_fill_gradient(low = "grey95", high = "steelblue4", limits = c(0, 1)) +
  facet_wrap(~ paste("Cluster", cluster)) +
  labs(title = sprintf("Shared top-%d variables across clustering methods (k = %d)",
                       max(pair_overlap$n1), K_CHOSEN),
       subtitle = "cell: number of shared variables (Jaccard index)",
       x = NULL, y = NULL, fill = "Jaccard") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        panel.grid = element_blank())

print(p_hm)

plot_dir <- ensure_dir(DIR_FIG, "04_profiles")
ggsave(file.path(plot_dir, sprintf("top_overlap_jaccard_%s.png", tag_k)),
       p_hm, width = 8, height = 6, dpi = 300)

# -- Save ----------------------------------------------------------------------
saveRDS(pair_overlap,   file.path(OUT_LOG, sprintf("top_overlap_pairwise_%s.rds", tag_k)))
saveRDS(vote_table,     file.path(OUT_LOG, sprintf("top_overlap_votes_%s.rds",    tag_k)))
saveRDS(core_table,     file.path(OUT_LOG, sprintf("top_overlap_core_%s.rds",     tag_k)))
saveRDS(pooled_overlap, file.path(OUT_LOG, sprintf("top_overlap_pooled_%s.rds",   tag_k)))

cat("\nDone. Objects in:", OUT_LOG, "| plot in:", plot_dir, "\n")
