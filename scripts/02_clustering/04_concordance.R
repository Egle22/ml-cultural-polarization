# =============================================================================
# scripts/02_clustering/04_concordance.R
# Purpose : Concordance between clustering configurations for every k.
#   1. Build all_clusters from the three source objects (UMAP k-means, PAM
#      Gower, PCA k-means on 2, 5, 15, 75, 198 PCs).
#   2. For each pair of configurations and each k: global agreement (after
#      Hungarian alignment) and Cramér's V.
#   3. Per-cluster agreement (how far each cluster of one configuration
#      matches the aligned cluster of the other).
#   4. Global and per-cluster heatmaps; everything saved as .rds and PDF.
# Thesis  : Section 3.5, 4.1.1; Table 4.2
# Input   : results/objects/clusters/umap_df_clusters.rds
#           results/objects/{pam_results_allk,pca_kmeans_objects}.rds
# Output  : results/objects/clusters/all_clusters.rds (raw labels, k = 2..8, 16)
#           results/objects/{comparison_global,comparison_cluster,
#           tabs_agreement_allk,tabs_cramerv_allk}.rds
#           results/figures/02_clustering/comparison/*.pdf
# =============================================================================

source(here::here("R", "00_setup.R"))

out_dir_obj  <- DIR_OBJ
out_dir_comp <- ensure_dir(DIR_FIG, "02_clustering", "comparison")

# -- 1. Source objects --------------------------------------------------------
umap_df     <- readRDS(file.path(DIR_CLUST, "umap_df_clusters.rds"))
results_pam <- readRDS(file.path(out_dir_obj, "pam_results_allk.rds"))
kmeans_pca  <- readRDS(file.path(out_dir_obj, "pca_kmeans_objects.rds"))

# -- 2. Build all_clusters ----------------------------------------------------
k_values  <- c(2, 3, 4, 5, 6, 7, 8, 16)
n_pc_list <- c(2, 5, 15, 75, 198)

n_obs        <- nrow(umap_df)
all_clusters <- data.frame(id = seq_len(n_obs))

# K-means UMAP
for (k in k_values) {
  col_name <- paste0("km_umap_", k)
  all_clusters[[col_name]] <- as.integer(umap_df[[paste0("cluster_k", k)]])
  cat(col_name, "— n:", length(all_clusters[[col_name]]), "\n")
}

# PAM Gower
for (k in k_values) {
  col_name <- paste0("pam_", k)
  all_clusters[[col_name]] <- as.integer(
    results_pam[[as.character(k)]]$clustering
  )
  cat(col_name, "— n:", length(all_clusters[[col_name]]), "\n")
}

# K-means PCA
for (n_pc in n_pc_list) {
  for (k in k_values) {
    col_name <- paste0("km_pc", n_pc, "_k", k)
    all_clusters[[col_name]] <- as.integer(
      kmeans_pca[[paste0("pc", n_pc)]][[paste0("k", k)]]$cluster
    )
    cat(col_name, "— n:", length(all_clusters[[col_name]]), "\n")
  }
}

all_clusters$id <- NULL
str(all_clusters)

saveRDS(all_clusters, file.path(ensure_dir(DIR_CLUST), "all_clusters.rds"))

# -- 3-5. Concordance for each k (helpers in R/alignment.R) -------------------
conc        <- compute_concordance(all_clusters, k_values)
global_long <- conc$global_long
clust_long  <- conc$clust_long

keep_methods <- c("km_umap", "pam", "km_pc5")

clust_long_5 <- clust_long %>%
  filter(sub("_k?[0-9]+$", "", method_i) %in% keep_methods,
         sub("_k?[0-9]+$", "", method_j) %in% keep_methods)

print(clust_long_5, row.names = FALSE)

# -- 6. Printed tables for each k ---------------------------------------------
make_matrix_table <- function(df, metric, k_val) {
  df %>%
    filter(k == k_val) %>%
    select(label_i, label_j, value = all_of(metric)) %>%
    pivot_wider(names_from = label_j, values_from = value) %>%
    rename(Method = label_i)
}

for (k_val in k_values) {
  
  cat(sprintf("\n══════════════════════════════════════\n"))
  cat(sprintf("  k = %d  —  Global Agreement\n", k_val))
  cat(sprintf("══════════════════════════════════════\n"))
  print(make_matrix_table(global_long, "Agreement", k_val), row.names = FALSE)
  
  cat(sprintf("\n══════════════════════════════════════\n"))
  cat(sprintf("  k = %d  —  Cramér's V\n", k_val))
  cat(sprintf("══════════════════════════════════════\n"))
  print(make_matrix_table(global_long, "CramerV", k_val), row.names = FALSE)
  
  cat(sprintf("\n══════════════════════════════════════\n"))
  cat(sprintf("  k = %d  —  Per-cluster Agreement\n", k_val))
  cat(sprintf("══════════════════════════════════════\n"))
  tab_cl <- clust_long %>%
    filter(k == k_val) %>%
    mutate(pair = paste0(label_i, " vs ", label_j)) %>%
    select(pair, cluster, n, agreement) %>%
    arrange(pair, as.integer(cluster))
  print(tab_cl, row.names = FALSE)
}

# -- 7. Global heatmap --------------------------------------------------------
plot_heatmap_global <- function(df_long, metric, k_val) {
  df <- df_long %>%
    filter(k == k_val) %>%
    select(label_i, label_j, value = all_of(metric))
  
  ggplot(df, aes(x = label_j, y = label_i, fill = value)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.2f", value)), size = 3, color = "grey20") +
    scale_fill_gradient2(
      low = "#f7f7f7", mid = "#74add1", high = "#313695",
      midpoint = 0.5, limits = c(0, 1), name = metric
    ) +
    labs(title = sprintf("%s — k = %d", metric, k_val), x = NULL, y = NULL) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x      = element_text(angle = 35, hjust = 1, size = 8),
      axis.text.y      = element_text(size = 8),
      panel.grid       = element_blank(),
      legend.key.width = unit(0.4, "cm")
    )
}

# -- 8. Per-cluster heatmap ---------------------------------------------------
plot_heatmap_cluster <- function(df_long, k_val) {
  df <- df_long %>%
    filter(k == k_val) %>%
    mutate(
      pair    = paste0(label_i, " vs\n", label_j),
      cluster = factor(cluster, levels = sort(unique(as.integer(cluster))))
    )
  
  ggplot(df, aes(x = pair, y = cluster, fill = agreement)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.2f", agreement)),
              size = 2.8, color = "grey20") +
    scale_fill_gradient2(
      low = "#d73027", mid = "#ffffbf", high = "#1a9850",
      midpoint = 0.5, limits = c(0, 1), name = "Agreement"
    ) +
    labs(
      title    = sprintf("Per-cluster agreement — k = %d", k_val),
      subtitle = "Rows = clusters of method_i  |  columns = method pairs",
      x = NULL, y = "Cluster"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x      = element_text(angle = 35, hjust = 1, size = 7),
      axis.text.y      = element_text(size = 9),
      panel.grid       = element_blank(),
      legend.key.width = unit(0.4, "cm")
    )
}

# -- 9. Save plots -------------------------------------------------------------

for (metric in c("Agreement", "CramerV")) {
  plots_hm <- lapply(k_values, function(k_val)
    plot_heatmap_global(global_long, metric, k_val))
  names(plots_hm) <- paste0("k", k_values)
  
  grid_hm <-
    (plots_hm$k2 | plots_hm$k3 | plots_hm$k4 | plots_hm$k5) /
    (plots_hm$k6 | plots_hm$k7 | plots_hm$k8 | plots_hm$k16) +
    plot_annotation(
      title   = sprintf("Pairwise %s — all methods and k values", metric),
      caption = "Methods: KM-UMAP, PAM, KM-PC2, KM-PC5, KM-PC15, KM-PC75, KM-PC198",
      theme   = theme(plot.title = element_text(size = 13, face = "bold"))
    )
  
  ggsave(
    filename = file.path(out_dir_comp,
                         sprintf("comparison_%s_heatmap.pdf", tolower(metric))),
    plot = grid_hm, width = 18, height = 9
  )
  cat("Saved heatmap:", metric, "\n")
}

for (k_val in k_values) {
  p_cl <- plot_heatmap_cluster(clust_long, k_val)
  ggsave(
    filename = file.path(out_dir_comp,
                         sprintf("comparison_cluster_agreement_k%d.pdf", k_val)),
    plot   = p_cl,
    width  = 14,
    height = 3 + k_val * 0.6
  )
  cat("Saved per-cluster agreement k =", k_val, "\n")
}

# -- 10. Save objects and tables ----------------------------------------------
saveRDS(global_long, file.path(out_dir_obj, "comparison_global.rds"))
saveRDS(clust_long,  file.path(out_dir_obj, "comparison_cluster.rds"))

tabs_agr <- lapply(k_values, function(k_val)
  make_matrix_table(global_long, "Agreement", k_val))
tabs_crv <- lapply(k_values, function(k_val)
  make_matrix_table(global_long, "CramerV", k_val))
names(tabs_agr) <- paste0("k", k_values)
names(tabs_crv) <- paste0("k", k_values)

saveRDS(tabs_agr, file.path(out_dir_obj, "tabs_agreement_allk.rds"))
saveRDS(tabs_crv, file.path(out_dir_obj, "tabs_cramerv_allk.rds"))

cat("All objects saved to:", out_dir_obj, "\n")