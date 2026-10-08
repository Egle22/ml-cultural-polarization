# =============================================================================
# scripts/02_clustering/06_hdbscan_vs_kmeans.R
# Purpose : (1) K-means (k = 6) on the 5D UMAP embedding.
#           (2) HDBSCAN on the 2D UMAP embedding; the number of clusters
#               emerges from density (no k is imposed).
#           (3) On the 2D embedding, compare the stored KM-UMAP labels
#               (all_clusters) with the HDBSCAN labels at the k selected by
#               HDBSCAN: agreement %, scatter and alluvial, Hungarian alignment.
#           (4) Same comparison on the 5D embedding.
#           (5) Mutual-reachability distance histograms.
#           Embeddings are NOT refitted: they are read from results/objects/UMAP/.
# Thesis  : Section 4.1.2; Figure 4.5; Table 4.1
# Input   : results/objects/UMAP/{umap_2d,umap_5d}.rds (from 05_hdbscan.R)
#           results/objects/clusters/all_clusters.rds  (from 04_concordance.R)
# Output  : results/objects/UMAP/kmeans_umap5d_k6.rds
#           results/objects/HDBSCAN/{hdbscan_umap2d_mp25,hdbscan_umap5d_mp25,
#           hdbscan_summary_2d_vs_5d}.rds
#           results/figures/02_clustering/hdbscan_comparison/*.pdf
# Seeds   : SEED = 42
# Note    : all_clusters.rds is built by 04_concordance.R, so run that script
#           first when rebuilding from scratch.
# =============================================================================

source(here::here("R", "00_setup.R"))
library(dbscan)

# -- 0. Configuration ---------------------------------------------------------
DIR_UMAP <- file.path(DIR_OBJ, "UMAP")
DIR_HDB  <- ensure_dir(DIR_OBJ, "HDBSCAN")
DIR_PLT  <- ensure_dir(DIR_FIG, "02_clustering", "hdbscan_comparison", "kmumap2d_vs_hdbscan2d")

SEED        <- 42
K_5D        <- 6      # only k fitted on the 5D embedding
MIN_PTS     <- 25     # same minPts used throughout the HDBSCAN work
KM_NSTART   <- 50
KM_ITER_MAX <- 100

theme_base <- theme_minimal(base_size = 11) +
  theme(plot.title       = element_text(face = "bold", size = 11),
        plot.subtitle    = element_text(color = "grey50", size = 9),
        panel.grid.minor = element_blank())

# ── 1. HELPERS ──────────────────────────────────────────────────────────────

#' Extract a numeric coordinate matrix from whatever the .rds happens to hold
#' (umap object, data.frame with UMAP* columns, or a plain matrix).
as_embedding <- function(obj) {
  if (inherits(obj, "umap"))      obj <- obj$layout
  if (is.data.frame(obj)) {
    keep <- vapply(obj, is.numeric, logical(1))
    obj  <- as.matrix(obj[, keep, drop = FALSE])
  }
  obj <- as.matrix(obj)
  storage.mode(obj) <- "double"
  obj
}

#' Hungarian relabelling of `oth` onto `ref` (assignment, not clustering).
align_hungarian <- function(ref, oth, k_val) {
  tab <- table(factor(ref, levels = as.character(seq_len(k_val))),
               factor(oth, levels = as.character(seq_len(k_val))))
  asg <- clue::solve_LSAP(tab, maximum = TRUE)
  map <- setNames(as.character(seq_len(k_val)), as.character(asg))
  as.character(map[as.character(oth)])
}

# ── 2. LOAD EMBEDDINGS ──────────────────────────────────────────────────────
message("[1] Loading saved embeddings")

emb_2d <- as_embedding(readRDS(file.path(DIR_UMAP, "umap_2d.rds")))
emb_5d <- as_embedding(readRDS(file.path(DIR_UMAP, "umap_5d.rds")))

stopifnot(ncol(emb_2d) == 2, ncol(emb_5d) == 5,
          nrow(emb_2d) == nrow(emb_5d))
n_obs <- nrow(emb_2d)
message(sprintf("    n = %d", n_obs))

all_clusters <- readRDS(file.path(DIR_OBJ, "clusters", "all_clusters.rds"))
stopifnot(nrow(all_clusters) == n_obs)

# ── 3. K-MEANS ON THE 5D EMBEDDING, k = 6 ONLY ─────────────────────────────
message("[2] K-means on the 5D embedding, k = 6")

set.seed(SEED)
km_5d <- stats::kmeans(emb_5d, centers = K_5D,
                       nstart = KM_NSTART, iter.max = KM_ITER_MAX)

lab_km5d <- as.character(km_5d$cluster)
cat("\nCluster sizes, KM-UMAP-5D (k = 6):\n")
print(table(lab_km5d))
print(table(all_clusters$km_umap_6))
# ── 4. HDBSCAN ON THE 2D EMBEDDING ─────────────────────────────────────────
message("[3] HDBSCAN on the 2D embedding")

set.seed(SEED)
hdb_2d <- dbscan::hdbscan(emb_2d, minPts = MIN_PTS)

lab_hdb    <- hdb_2d$cluster                    # 0 = noise
k_opt      <- length(setdiff(unique(lab_hdb), 0))
noise_frac <- mean(lab_hdb == 0)

cat("\n─── HDBSCAN on UMAP-2D (minPts =", MIN_PTS, ") ───\n")
cat("Clusters found (excluding noise): ", k_opt, "\n", sep = "")
cat("Noise fraction: ", sprintf("%.1f%%", 100 * noise_frac), "\n", sep = "")
cat("Cluster sizes (0 = noise):\n")
print(table(lab_hdb))
cat("Cluster stability scores:\n")
print(round(hdb_2d$cluster_scores, 3))

if (k_opt < 2) {
  stop("HDBSCAN returned fewer than 2 clusters on the 2D embedding — ",
       "there is no k at which a comparison against K-means is defined.")
}

# ── 5. SELECT THE MATCHING K-MEANS COLUMN FROM all_clusters ────────────────
message("[4] Comparison at k = ", k_opt)

cand <- c(paste0("km_umap_", k_opt), paste0("km_umap_k", k_opt))
col_km <- cand[cand %in% names(all_clusters)][1]

if (is.na(col_km)) {
  stop("No K-means-on-UMAP-2D column for k = ", k_opt,
       " in all_clusters. Looked for: ", paste(cand, collapse = ", "))
}
message("    using column: ", col_km)

lab_km2d <- as.character(all_clusters[[col_km]])

# ── 6. AGREEMENT (noise points excluded) ───────────────────────────────────
# HDBSCAN assigns a point to no cluster at all when it sits below the density
# threshold. Those points have no K-means counterpart to agree or disagree with,
# so the agreement rate is computed on the clustered subset only. The noise
# share is reported separately and is itself part of the finding.

keep <- lab_hdb != 0

lab_hdb_k  <- as.character(lab_hdb[keep])
lab_km_k   <- lab_km2d[keep]
lab_hdb_al <- align_hungarian(lab_km_k, lab_hdb_k, k_opt)

pct_agree <- round(100 * mean(lab_km_k == lab_hdb_al), 1)
pct_disc  <- round(100 - pct_agree, 1)

cat("\nAgreement on the clustered subset (n = ", sum(keep), "):\n", sep = "")
cat("  Agreement : ", pct_agree, "%\n", sep = "")
cat("  Discordant: ", pct_disc,  "%\n", sep = "")

conf_tab <- table(
  KMeans  = factor(lab_km_k,   levels = as.character(seq_len(k_opt))),
  HDBSCAN = factor(lab_hdb_al, levels = as.character(seq_len(k_opt)))
)
cat("\nConfusion table (after Hungarian alignment):\n")
print(conf_tab)

# ── 7. PLOT DATA ────────────────────────────────────────────────────────────
lev <- as.character(seq_len(k_opt))

df_plot <- data.frame(
  UMAP1      = emb_2d[keep, 1],
  UMAP2      = emb_2d[keep, 2],
  cl_ref     = factor(lab_km_k,   levels = lev),   # KMeans on UMAP-2D
  cl_oth     = factor(lab_hdb_al, levels = lev),   # HDBSCAN on UMAP-2D
  discordant = lab_km_k != lab_hdb_al
)

df_noise <- data.frame(
  UMAP1 = emb_2d[!keep, 1],
  UMAP2 = emb_2d[!keep, 2]
)

ref_lab <- "KM-UMAP-2D"
oth_lab <- "HDBSCAN-2D"

# ── 8. PLOTS ────────────────────────────────────────────────────────────────
p_km <- ggplot(df_plot, aes(UMAP1, UMAP2)) +
  geom_point(data = df_noise, colour = "grey85", alpha = 0.35, size = 0.7) +
  geom_point(data = filter(df_plot, !discordant), aes(color = cl_ref),
             alpha = 0.45, size = 0.9) +
  geom_point(data = filter(df_plot, discordant),
             color = "black", alpha = 0.6, size = 1.2) +
  scale_color_manual(values = CLUSTER_COLORS[seq_len(k_opt)],
                     name = paste0(ref_lab, "\ncluster")) +
  labs(title    = paste0("UMAP-2D space — ", ref_lab, " labels  (K = ", k_opt, ")"),
       subtitle = paste0("Black = discordant with ", oth_lab, "  [", pct_disc,
                         "%]  |  Grey = HDBSCAN noise"),
       x = "UMAP 1", y = "UMAP 2") +
  theme_base

p_hdb <- ggplot(df_plot, aes(UMAP1, UMAP2)) +
  geom_point(data = df_noise, colour = "grey85", alpha = 0.35, size = 0.7) +
  geom_point(data = filter(df_plot, !discordant), aes(color = cl_oth),
             alpha = 0.45, size = 0.9) +
  geom_point(data = filter(df_plot, discordant),
             color = "black", alpha = 0.6, size = 1.2) +
  scale_color_manual(values = CLUSTER_COLORS[seq_len(k_opt)],
                     name = paste0(oth_lab, "\ncluster")) +
  labs(title    = paste0("UMAP-2D space — ", oth_lab, " labels  (K = ", k_opt, ")"),
       subtitle = paste0("Black = discordant with ", ref_lab, "  [", pct_disc,
                         "%]  |  Grey = HDBSCAN noise (",
                         sprintf("%.1f%%", 100 * noise_frac), ")"),
       x = "UMAP 1", y = "UMAP 2") +
  theme_base

df_al <- as.data.frame(table(
  ref = factor(df_plot$cl_ref, levels = lev),
  oth = factor(df_plot$cl_oth, levels = lev)
)) |>
  rename(n = Freq) |>
  filter(n > 0)

p_alluvial <- ggplot(df_al, aes(axis1 = ref, axis2 = oth, y = n)) +
  geom_alluvium(aes(fill = ref), alpha = 0.75, width = 0.3) +
  geom_stratum(aes(fill = after_stat(stratum)),
               width = 0.3, color = "white", linewidth = 0.4) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),
            size = 3.5, fontface = "bold", color = "white") +
  scale_x_discrete(limits = c("ref", "oth"),
                   labels = c(ref_lab, oth_lab), expand = c(0.15, 0.15)) +
  scale_fill_manual(values = CLUSTER_COLORS[seq_len(k_opt)]) +
  labs(title = paste0("Cluster flow: ", ref_lab, " <-> ", oth_lab,
                      "  (K = ", k_opt, ")"),
       y = "Observations (noise excluded)") +
  theme_base +
  theme(legend.position = "none", panel.grid = element_blank())

panel <- (p_km | p_alluvial | p_hdb) +
  plot_annotation(
    title   = paste0(ref_lab, "  vs  ", oth_lab, "  —  K = ", k_opt),
    caption = paste0("Reference: ", ref_lab, "  |  Hungarian alignment  |  ",
                     pct_agree, "% agreement on the clustered subset  |  ",
                     sprintf("%.1f%%", 100 * noise_frac), " HDBSCAN noise"),
    theme = theme(plot.title   = element_text(face = "bold", size = 13, hjust = 0.5),
                  plot.caption = element_text(color = "grey50", size = 8))
  )

# ── 9. CONSOLIDATED SAVES ───────────────────────────────────────────────────
message("[5] Saving")

ggplot2::ggsave(file.path(DIR_PLT, sprintf("panel_kmumap2d_vs_hdbscan2d_k%d.pdf", k_opt)),
                panel, width = 16, height = 5.5, bg = "white")

saveRDS(
  list(k = K_5D, seed = SEED, labels = lab_km5d,
       centers = km_5d$centers, tot_withinss = km_5d$tot.withinss),
  file.path(DIR_UMAP, "kmeans_umap5d_k6.rds")
)

saveRDS(
  list(min_pts        = MIN_PTS,
       hdbscan        = hdb_2d,
       labels         = lab_hdb,
       n_clusters     = k_opt,
       noise_frac     = noise_frac,
       cluster_scores = hdb_2d$cluster_scores,
       km_column_used = col_km,
       pct_agreement  = pct_agree,
       confusion      = conf_tab),
  file.path(DIR_HDB, "hdbscan_umap2d_mp25.rds")
)

cat("\nObjects -> ", DIR_UMAP, " / ", DIR_HDB, "\n", sep = "")
cat("Plots   -> ", DIR_PLT, "\n", sep = "")

# 5d ----
# ═══════════════════════════════════════════════════════════════════════════
#  PART 5D — HDBSCAN ON THE 5D EMBEDDING, COMPARED AGAINST K-MEANS
#
#  Mirrors the 2D comparison above. Two differences, both forced by the data:
#   - there are no stored K-means labels for the 5D embedding in all_clusters,
#     so K-means is fitted here at the k selected by HDBSCAN;
#   - the scatter plots use the first two of the five UMAP dimensions as a
#     display space. The partition itself lives in 5D.
# ═══════════════════════════════════════════════════════════════════════════

DIR_PLT5 <- ensure_dir(DIR_FIG, "02_clustering", "hdbscan_comparison", "kmumap5d_vs_hdbscan5d")

# ── 10. HDBSCAN ON THE 5D EMBEDDING ────────────────────────────────────────
message("[6] HDBSCAN on the 5D embedding")

set.seed(SEED)
hdb_5d <- dbscan::hdbscan(emb_5d, minPts = MIN_PTS)

lab_hdb5    <- hdb_5d$cluster                   # 0 = noise
k_opt5      <- length(setdiff(unique(lab_hdb5), 0))
noise_frac5 <- mean(lab_hdb5 == 0)

cat("\n─── HDBSCAN on UMAP-5D (minPts =", MIN_PTS, ") ───\n")
cat("Clusters found (excluding noise): ", k_opt5, "\n", sep = "")
cat("Noise fraction: ", sprintf("%.1f%%", 100 * noise_frac5), "\n", sep = "")
cat("Cluster sizes (0 = noise):\n")
print(table(lab_hdb5))
cat("Cluster stability scores:\n")
print(round(hdb_5d$cluster_scores, 3))

if (k_opt5 < 2) {
  stop("HDBSCAN returned fewer than 2 clusters on the 5D embedding — ",
       "there is no k at which a comparison against K-means is defined. ",
       "This is itself the substantive result: report the noise fraction.")
}

# ── 11. K-MEANS ON THE 5D EMBEDDING AT k = k_opt5 ──────────────────────────
message("[7] K-means on the 5D embedding at k = ", k_opt5)

if (k_opt5 == K_5D) {
  message("    k matches K_5D — reusing the existing km_5d fit")
  lab_km5 <- lab_km5d
} else {
  set.seed(SEED)
  km_5d_kopt <- stats::kmeans(emb_5d, centers = k_opt5,
                              nstart = KM_NSTART, iter.max = KM_ITER_MAX)
  lab_km5 <- as.character(km_5d_kopt$cluster)
}

cat("\nCluster sizes, KM-UMAP-5D (k = ", k_opt5, "):\n", sep = "")
print(table(lab_km5))

# ── 12. AGREEMENT (noise points excluded) ──────────────────────────────────
# Same convention as the 2D case: HDBSCAN noise points have no K-means
# counterpart to agree or disagree with, so agreement is computed on the
# clustered subset and the noise share is reported separately.

keep5 <- lab_hdb5 != 0

lab_hdb5_k  <- as.character(lab_hdb5[keep5])
lab_km5_k   <- lab_km5[keep5]
lab_hdb5_al <- align_hungarian(lab_km5_k, lab_hdb5_k, k_opt5)

pct_agree5 <- round(100 * mean(lab_km5_k == lab_hdb5_al), 1)
pct_disc5  <- round(100 - pct_agree5, 1)

cat("\nAgreement on the clustered subset (n = ", sum(keep5), "):\n", sep = "")
cat("  Agreement : ", pct_agree5, "%\n", sep = "")
cat("  Discordant: ", pct_disc5,  "%\n", sep = "")

conf_tab5 <- table(
  KMeans  = factor(lab_km5_k,   levels = as.character(seq_len(k_opt5))),
  HDBSCAN = factor(lab_hdb5_al, levels = as.character(seq_len(k_opt5)))
)
cat("\nConfusion table (after Hungarian alignment):\n")
print(conf_tab5)

# ── 13. PLOT DATA ──────────────────────────────────────────────────────────
lev5 <- as.character(seq_len(k_opt5))

df_plot5 <- data.frame(
  UMAP1      = emb_5d[keep5, 1],
  UMAP2      = emb_5d[keep5, 2],
  cl_ref     = factor(lab_km5_k,   levels = lev5),   # KMeans on UMAP-5D
  cl_oth     = factor(lab_hdb5_al, levels = lev5),   # HDBSCAN on UMAP-5D
  discordant = lab_km5_k != lab_hdb5_al
)

df_noise5 <- data.frame(
  UMAP1 = emb_5d[!keep5, 1],
  UMAP2 = emb_5d[!keep5, 2]
)

ref_lab5 <- "KM-UMAP-5D"
oth_lab5 <- "HDBSCAN-5D"

# ── 14. PLOTS ──────────────────────────────────────────────────────────────
p_km5 <- ggplot(df_plot5, aes(UMAP1, UMAP2)) +
  geom_point(data = df_noise5, colour = "grey85", alpha = 0.35, size = 0.7) +
  geom_point(data = filter(df_plot5, !discordant), aes(color = cl_ref),
             alpha = 0.45, size = 0.9) +
  geom_point(data = filter(df_plot5, discordant),
             color = "black", alpha = 0.6, size = 1.2) +
  scale_color_manual(values = CLUSTER_COLORS[seq_len(k_opt5)],
                     name = paste0(ref_lab5, "\ncluster")) +
  labs(title    = paste0("UMAP-5D space — ", ref_lab5, " labels  (K = ", k_opt5, ")"),
       subtitle = paste0("Dims 1-2 of 5 shown  |  Black = discordant with ",
                         oth_lab5, "  [", pct_disc5, "%]  |  Grey = noise"),
       x = "UMAP dim 1", y = "UMAP dim 2") +
  theme_base

p_hdb5 <- ggplot(df_plot5, aes(UMAP1, UMAP2)) +
  geom_point(data = df_noise5, colour = "grey85", alpha = 0.35, size = 0.7) +
  geom_point(data = filter(df_plot5, !discordant), aes(color = cl_oth),
             alpha = 0.45, size = 0.9) +
  geom_point(data = filter(df_plot5, discordant),
             color = "black", alpha = 0.6, size = 1.2) +
  scale_color_manual(values = CLUSTER_COLORS[seq_len(k_opt5)],
                     name = paste0(oth_lab5, "\ncluster")) +
  labs(title    = paste0("UMAP-5D space — ", oth_lab5, " labels  (K = ", k_opt5, ")"),
       subtitle = paste0("Dims 1-2 of 5 shown  |  Black = discordant with ",
                         ref_lab5, "  [", pct_disc5, "%]  |  Grey = noise (",
                         sprintf("%.1f%%", 100 * noise_frac5), ")"),
       x = "UMAP dim 1", y = "UMAP dim 2") +
  theme_base

df_al5 <- as.data.frame(table(
  ref = factor(df_plot5$cl_ref, levels = lev5),
  oth = factor(df_plot5$cl_oth, levels = lev5)
)) |>
  rename(n = Freq) |>
  filter(n > 0)

p_alluvial5 <- ggplot(df_al5, aes(axis1 = ref, axis2 = oth, y = n)) +
  geom_alluvium(aes(fill = ref), alpha = 0.75, width = 0.3) +
  geom_stratum(aes(fill = after_stat(stratum)),
               width = 0.3, color = "white", linewidth = 0.4) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),
            size = 3.5, fontface = "bold", color = "white") +
  scale_x_discrete(limits = c("ref", "oth"),
                   labels = c(ref_lab5, oth_lab5), expand = c(0.15, 0.15)) +
  scale_fill_manual(values = CLUSTER_COLORS[seq_len(k_opt5)]) +
  labs(title = paste0("Cluster flow: ", ref_lab5, " <-> ", oth_lab5,
                      "  (K = ", k_opt5, ")"),
       y = "Observations (noise excluded)") +
  theme_base +
  theme(legend.position = "none", panel.grid = element_blank())

panel5 <- (p_km5 | p_alluvial5 | p_hdb5) +
  plot_annotation(
    title   = paste0(ref_lab5, "  vs  ", oth_lab5, "  —  K = ", k_opt5),
    caption = paste0("Reference: ", ref_lab5, "  |  Hungarian alignment  |  ",
                     pct_agree5, "% agreement on the clustered subset  |  ",
                     sprintf("%.1f%%", 100 * noise_frac5), " HDBSCAN noise"),
    theme = theme(plot.title   = element_text(face = "bold", size = 13, hjust = 0.5),
                  plot.caption = element_text(color = "grey50", size = 8))
  )

# ── 15. CONSOLIDATED SAVES (5D) ────────────────────────────────────────────
message("[8] Saving 5D comparison")

ggplot2::ggsave(file.path(DIR_PLT5, sprintf("panel_kmumap5d_vs_hdbscan5d_k%d.pdf", k_opt5)),
                panel5, width = 16, height = 5.5, bg = "white")

saveRDS(
  list(min_pts        = MIN_PTS,
       hdbscan        = hdb_5d,
       labels         = lab_hdb5,
       n_clusters     = k_opt5,
       noise_frac     = noise_frac5,
       cluster_scores = hdb_5d$cluster_scores,
       kmeans_labels  = lab_km5,
       kmeans_k       = k_opt5,
       pct_agreement  = pct_agree5,
       confusion      = conf_tab5),
  file.path(DIR_HDB, "hdbscan_umap5d_mp25.rds")
)

# ── 16. SIDE-BY-SIDE SUMMARY, 2D vs 5D ─────────────────────────────────────
summary_2d_5d <- data.frame(
  embedding      = c("UMAP 2D", "UMAP 5D"),
  hdbscan_k      = c(k_opt, k_opt5),
  noise_pct      = c(100 * noise_frac, 100 * noise_frac5),
  n_clustered    = c(sum(keep), sum(keep5)),
  agreement_pct  = c(pct_agree, pct_agree5)
)

cat("\n───────────────────────────────────────────────────────────────\n")
cat(" HDBSCAN vs K-means: 2D versus 5D embedding\n")
cat("───────────────────────────────────────────────────────────────\n")
print(summary_2d_5d, row.names = FALSE, digits = 3)

saveRDS(summary_2d_5d, file.path(DIR_HDB, "hdbscan_summary_2d_vs_5d.rds"))

cat("\nPlots -> ", DIR_PLT5, "\n", sep = "")
# Mutual reachability ----
# ── MUTUAL REACHABILITY DISTANCE HISTOGRAM ─────────────────────────────────
# d_mreach(a,b) = max(d_core(a), d_core(b), d(a,b))


EMB      <- emb_5d      # or emb_2d
MIN_PTS  <- 25
TAG      <- "umap5d"    # or "umap2d"

# ── 1. Core distances ───────────────────────────────────────────────────────
core_d <- dbscan::kNNdist(EMB, k = MIN_PTS - 1)
if (is.matrix(core_d)) core_d <- core_d[, ncol(core_d)]

# ── 2. Full mutual reachability matrix ──────────────────────────────────────
D    <- as.matrix(stats::dist(EMB))
Cmax <- pmax(outer(core_d, rep(1, length(core_d))),
             outer(rep(1, length(core_d)), core_d))
Dmr  <- pmax(D, Cmax)

iu       <- upper.tri(Dmr)
mr_pairs <- Dmr[iu]
d_pairs  <- D[iu]

# fraction of pairs whose distance is set by a core distance, not by geometry
core_dominance <- mean(mr_pairs > d_pairs)
cat("Core-dominance fraction: ", sprintf("%.1f%%", 100 * core_dominance), "\n", sep = "")

# ── 3. MST edge weights (the density-relevant subset) ──────────────────────
hdb      <- dbscan::hdbscan(EMB, minPts = MIN_PTS)
mst_w    <- hdb$hc$height      # merge heights = MST edge weights in d_mreach

# ── 4. Plots ────────────────────────────────────────────────────────────────
theme_h <- theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(color = "grey50", size = 9),
        panel.grid.minor = element_blank())

p_all <- ggplot(data.frame(d = mr_pairs), aes(d)) +
  geom_histogram(bins = 80, fill = "#4C7896", colour = "white", linewidth = 0.15) +
  labs(title    = "Mutual reachability distance — all pairs",
       subtitle = sprintf("n pairs = %s  |  core-dominant: %.1f%%",
                          format(length(mr_pairs), big.mark = ","),
                          100 * core_dominance),
       x = expression(d[mreach]), y = "Frequency") +
  theme_h

p_mst <- ggplot(data.frame(d = mst_w), aes(d)) +
  geom_histogram(bins = 60, fill = "#B5705B", colour = "white", linewidth = 0.15) +
  labs(title    = "Mutual reachability distance — MST edges only",
       subtitle = "Merge heights of the single-linkage hierarchy",
       x = expression(d[mreach]), y = "Frequency") +
  theme_h

p_core <- ggplot(data.frame(d = core_d), aes(d)) +
  geom_histogram(bins = 60, fill = "#5A9E6F", colour = "white", linewidth = 0.15) +
  labs(title    = "Core distance",
       subtitle = sprintf("minPts = %d  |  CV = %.2f",
                          MIN_PTS, sd(core_d) / mean(core_d)),
       x = expression(d[core]), y = "Frequency") +
  theme_h

panel_mr <- (p_core | p_mst) / p_all
panel_mr <- (p_all | p_mst) 
ggsave(file.path(DIR_PLT, sprintf("mreach_histograms_%s.pdf", TAG)),
       panel_mr, width = 12, height = 8, bg = "white")
