# =============================================================================
# scripts/02_clustering/07_align_clusters.R
# Purpose : Align the labels of every partition to a reference configuration
#           (KM-PC5) with the Hungarian algorithm, so that "cluster 1" means
#           the same group across methods; recompute the concordance on the
#           aligned labels; save the CANONICAL cluster-label file.
# Thesis  : Sections 3.5, 4.1.3; Tables 4.2, A.3
# Input   : results/objects/clusters/all_clusters.rds       (04_concordance.R)
#           results/objects/HDBSCAN/hdbscan_umap2d_mp25.rds (06_hdbscan_vs_kmeans.R)
#           data/processed/geo_raw.rds                       (respondent IDs)
# Output  : results/objects/clusters/all_clusters_aligned.rds  <- CANONICAL
#             37 columns, one row per respondent (n = 2,579):
#             km_umap_{2..6}, pam_{2..6}, km_pc{2,5,15,75,198}_k{2..6},
#             hdb_umap2_mp25 (raw HDBSCAN labels, 0 = noise), D_INTERVIEW
#           results/objects/{comparison_global_aligned,comparison_cluster_aligned}.rds
#           results/figures/02_clustering/comparison/per_cluster_agreement_k4.pdf
# Note    : no other script may overwrite all_clusters_aligned.rds (CLAUDE.md).
# =============================================================================

source(here::here("R", "00_setup.R"))

out_dir_obj  <- DIR_OBJ
out_dir_comp <- ensure_dir(DIR_FIG, "02_clustering", "comparison")

all_clusters <- readRDS(file.path(DIR_CLUST, "all_clusters.rds"))
k_values     <- c(2, 3, 4, 5, 6, 7, 8, 16)

# -- 1. Global alignment of all partitions to a reference ---------------------
REF_METHOD <- "km_pc5"   # alignment anchor

ref_col_for_k <- function(k_val) {
  if (grepl("^km_pc", REF_METHOD)) paste0(REF_METHOD, "_k", k_val)
  else                              paste0(REF_METHOD, "_", k_val)
}

all_clusters_aligned <- all_clusters

for (k_val in k_values) {
  cols <- get_col_names(k_val)
  cols <- cols[cols %in% colnames(all_clusters)]

  ref_col <- ref_col_for_k(k_val)
  stopifnot(ref_col %in% cols)
  ref_vec <- all_clusters[[ref_col]]

  for (col in cols) {
    if (col == ref_col) {
      # the reference is unchanged (its labels define the scheme)
      all_clusters_aligned[[col]] <- as.integer(all_clusters[[col]])
      next
    }
    all_clusters_aligned[[col]] <- as.integer(
      hungarian_align(ref_vec, all_clusters[[col]], k_val)
    )
  }
  cat(sprintf("k = %d aligned to %s\n", k_val, ref_col))
}

# Check: agreement with the reference now reads on the diagonal
for (k_val in k_values) {
  ref_col <- ref_col_for_k(k_val)
  cols <- get_col_names(k_val)
  cols <- cols[cols %in% colnames(all_clusters_aligned) & cols != ref_col]
  for (col in cols) {
    agr <- mean(all_clusters_aligned[[ref_col]] == all_clusters_aligned[[col]])
    cat(sprintf("k=%2d  %-12s vs %-12s : %.3f\n", k_val, ref_col, col, agr))
  }
}

# -- 2. Concordance on the aligned labels (all k) -----------------------------
conc        <- compute_concordance(all_clusters_aligned, k_values)
global_long <- conc$global_long
clust_long  <- conc$clust_long

saveRDS(global_long, file.path(out_dir_obj, "comparison_global_aligned.rds"))
saveRDS(clust_long,  file.path(out_dir_obj, "comparison_cluster_aligned.rds"))

keep_methods <- c("km_umap", "pam", "km_pc5")
clust_long_5 <- clust_long %>%
  filter(sub("_k?[0-9]+$", "", method_i) %in% keep_methods,
         sub("_k?[0-9]+$", "", method_j) %in% keep_methods)
print(clust_long_5, row.names = FALSE)

# -- 3. Canonical file: thesis columns only -----------------------------------
# The thesis uses k = 2..6 (Table A.3). All PCA dimensions are kept. The
# HDBSCAN partition (UMAP 2D, minPts = 25) and the respondent ID are added.
K_THESIS  <- 2:6
keep_cols <- unlist(lapply(K_THESIS, get_col_names))
keep_cols <- names(all_clusters_aligned)[names(all_clusters_aligned) %in% keep_cols]

hdb_2d  <- readRDS(file.path(DIR_OBJ, "HDBSCAN", "hdbscan_umap2d_mp25.rds"))
geo_raw <- readRDS(file.path(DIR_DATA, "geo_raw.rds"))
stopifnot(length(hdb_2d$labels) == nrow(all_clusters_aligned),
          nrow(geo_raw)         == nrow(all_clusters_aligned))

clusters_canonical <- all_clusters_aligned[, keep_cols]
clusters_canonical[["hdb_umap2_mp25"]] <- as.integer(hdb_2d[["labels"]])
clusters_canonical[["D_INTERVIEW"]]    <- as.integer(geo_raw[["D_INTERVIEW"]])

stopifnot(ncol(clusters_canonical) == 37)
saveRDS(clusters_canonical, PATH_CLUSTERS)
cat("Saved", PATH_CLUSTERS, "\n")

# -- 4. Per-cluster agreement at k = 4 (KM-PC2, KM-PC5, KM-UMAP, PAM) ---------
k_val_plot   <- 4
keep_methods <- c("km_umap", "pam", "km_pc5", "km_pc2")

# column name for a method at k
get_col <- function(method, k) {
  if (method == "pam")     return(paste0("pam_", k))
  if (method == "km_umap") return(paste0("km_umap_", k))
  return(paste0(method, "_k", k))
}

method_label_map <- c(
  km_umap = "KM-UMAP",
  pam     = "PAM",
  km_pc5  = "KM-PC5",
  km_pc2  = "KM-PC2"
)

# agreement for each pair, with Hungarian alignment
pairs <- combn(keep_methods, 2, simplify = FALSE)

bar_rows <- list()

for (p in pairs) {
  m_i <- p[1]; m_j <- p[2]
  col_i <- get_col(m_i, k_val_plot)
  col_j <- get_col(m_j, k_val_plot)

  if (!col_i %in% colnames(all_clusters) ||
      !col_j %in% colnames(all_clusters)) next

  lbl_i <- as.character(all_clusters[[col_i]])
  lbl_j <- as.character(all_clusters[[col_j]])

  # Hungarian alignment: j aligned to i
  lbl_j_al <- hungarian_align(lbl_i, lbl_j, k_val_plot)

  lev <- as.character(seq_len(k_val_plot))

  for (cl in lev) {
    idx <- lbl_i == cl
    n_cl <- sum(idx)
    agr  <- if (n_cl == 0) NA else mean(lbl_j_al[idx] == cl)
    bar_rows[[length(bar_rows) + 1]] <- data.frame(
      method_i = m_i, method_j = m_j,
      label_i  = method_label_map[m_i],
      label_j  = method_label_map[m_j],
      cluster  = cl,
      n        = n_cl,
      agreement = agr,
      stringsAsFactors = FALSE
    )
  }
}

bar_df <- do.call(rbind, bar_rows) %>%
  mutate(
    pair    = paste(label_i, "vs", label_j),
    cluster = factor(as.integer(cluster),
                     levels = seq_len(k_val_plot))
  )

# facet order
pair_order <- c(
  "KM-PC2 vs KM-UMAP",
  "KM-PC2 vs PAM",
  "KM-PC2 vs KM-PC5",
  "KM-PC5 vs KM-UMAP",
  "KM-PC5 vs PAM",
  "KM-UMAP vs PAM"
)

bar_df <- bar_df %>%
  mutate(
    pair = ifelse(pair %in% pair_order, pair,
                  paste(label_j, "vs", label_i)),
    pair = factor(pair, levels = pair_order)
  ) %>%
  filter(!is.na(pair), !is.na(agreement))

pal_k4 <- setNames(cluster_colors[seq_len(k_val_plot)],
                   as.character(seq_len(k_val_plot)))

p_k4 <- ggplot(bar_df, aes(x = cluster, y = agreement, fill = cluster)) +
  geom_col(width = 0.85) +
  geom_text(aes(label = sprintf("%.2f", agreement)),
            vjust = -0.4, size = 3) +
  geom_text(aes(label = paste0("n=", n)),
            vjust = 1.3, size = 2.6, color = "white") +
  facet_wrap(~ pair, nrow = 2) +
  scale_fill_manual(values = pal_k4, name = "Cluster") +
  scale_y_continuous(limits = c(0, 1.12), breaks = seq(0, 1, 0.25)) +
  labs(
    title    = "Per-cluster label concordance at k = 4",
    subtitle = "Row 1: KM-PC2 vs others  |  Row 2: KM-PC5 vs others",
    x        = "Cluster",
    y        = "Per-cluster agreement"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    strip.text         = element_text(face = "bold"),
    legend.position    = "right"
  )

print(p_k4)

ggsave(
  file.path(out_dir_comp, "per_cluster_agreement_k4.pdf"),
  p_k4, width = 14, height = 8, bg = "white"
)
