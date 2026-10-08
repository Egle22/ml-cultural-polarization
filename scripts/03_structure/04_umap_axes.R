# =============================================================================
# scripts/03_structure/04_umap_axes.R
# Purpose : Interpret the two UMAP axes: correlation with the principal
#           components, UMAP coloured by PC1 / PC2 scores, correlation of each
#           item with each axis, random-forest importance of the items for
#           each axis, and comparison with the PCA loadings.
# Thesis  : Section 4.1.1, Appendix A.3.1; Figure 4.3 (a); Table A.5
# Input   : results/objects/clusters/umap_df_clusters.rds (2D UMAP + KM-UMAP labels)
#           results/objects/{pca_scores,pca_model}.rds
#           data/processed/data_imputed_selected.rds
# Output  : results/figures/03_structure/axis_interpretation/*.pdf
# Seeds   : set.seed(42); randomForest(ntree = 500)
# =============================================================================

source(here::here("R", "00_setup.R"))

PLOT_DIR <- ensure_dir(DIR_FIG, "03_structure", "axis_interpretation")

# -- 0. Load objects -----------------------------------------------------------
umap_df      <- readRDS(file.path(DIR_CLUST, "umap_df_clusters.rds"))
pca_scores   <- readRDS(file.path(DIR_OBJ, "pca_scores.rds"))
pca_model    <- readRDS(file.path(DIR_OBJ, "pca_model.rds"))
imputed_data <- readRDS(file.path(DIR_DATA, "data_imputed_selected.rds"))

umap_coords <- umap_df[, c("UMAP1", "UMAP2")]
cluster_k2  <- factor(umap_df$cluster_k2)

# -- 1. Correlation between UMAP and PCs ------------------------------------------
pc_coords <- data.frame(
  PC1 = pca_scores[, 1],
  PC2 = pca_scores[, 2],
  PC3 = pca_scores[, 3],
  PC4 = pca_scores[, 4],
  PC5 = pca_scores[, 5]
)

cor_mat <- cor(cbind(umap_coords, pc_coords))
cat("Correlation between UMAP axes and first 2 PCs:\n")
print(round(cor_mat[1:2, 3:4], 3))

# -- 2. UMAP coloured by PC1 / PC2 ---------------------------------------------
df_plot <- data.frame(
  UMAP1   = umap_coords$UMAP1,
  UMAP2   = umap_coords$UMAP2,
  PC1     = pca_scores[, 1],
  PC2     = pca_scores[, 2],
  cluster = cluster_k2
)

p_pc1 <- ggplot(df_plot, aes(x = UMAP1, y = UMAP2, color = PC1)) +
  geom_point(size = 0.8, alpha = 0.6) +
  scale_color_gradient2(low  = "#2C7BB6", mid = "grey90", high = "#D7191C",
                        midpoint = 0, name = "PC1 score") +
  labs(title = "UMAP 2D — colored by PC1 score") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())

p_pc2 <- ggplot(df_plot, aes(x = UMAP1, y = UMAP2, color = PC2)) +
  geom_point(size = 0.8, alpha = 0.6) +
  scale_color_gradient2(low  = "#2C7BB6", mid = "grey90", high = "#D7191C",
                        midpoint = 0, name = "PC2 score") +
  labs(title = "UMAP 2D — colored by PC2 score") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())

p_cluster <- ggplot(df_plot, aes(x = UMAP1, y = UMAP2, color = cluster)) +
  geom_point(size = 0.8, alpha = 0.6) +
  scale_color_manual(values = c("1" = "#4C7896", "2" = "#B5705B"),
                     name = "Cluster") +
  labs(title = "UMAP 2D — K-means clusters (k=2)") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())

(p_cluster | p_pc1 | p_pc2) +
  plot_annotation(
    title = "UMAP structure vs PCA scores",
    theme = theme(plot.title = element_text(face = "bold", hjust = 0.5))
  )

ggsave(file.path(PLOT_DIR, "umap_vs_pca_scores.pdf"), width = 14, height = 5)


# -- 3. Correlation of UMAP axes with the original items ------------------------
X_num <- imputed_data %>%
  mutate(across(where(is.factor), ~ as.numeric(as.character(.)))) %>%
  as.data.frame()

cor_umap1 <- sapply(X_num, function(v) cor(v, umap_coords$UMAP1,
                                           use = "complete.obs"))
cor_umap2 <- sapply(X_num, function(v) cor(v, umap_coords$UMAP2,
                                           use = "complete.obs"))

cor_vars_umap2 <- data.frame(
  variable = names(cor_umap2),
  cor_UMAP2 = cor_umap2
) %>%
  left_join(section_lookup, by = "variable") %>%
  arrange(desc(abs(cor_UMAP2)))

cat("\nTop 20 variables correlated with UMAP2:\n")
print(head(cor_vars_umap2[, c("variable", "section_label", "cor_UMAP2")], 20))


cor_vars_umap1 <- data.frame(
  variable = names(cor_umap1),
  cor_UMAP1 = cor_umap1
) %>%
  left_join(section_lookup, by = "variable") %>%
  arrange(desc(abs(cor_UMAP1)))

cat("\nTop 20 variables correlated with UMAP1:\n")
print(head(cor_vars_umap1[, c("variable", "section_label", "cor_UMAP1")], 20))


# -- 4. Random forest: variable importance for each UMAP axis --------------------

set.seed(42)

# UMAP2
rf_umap2 <- randomForest(
  x         = X_num,
  y         = umap_coords$UMAP2,
  ntree     = 500,
  importance = TRUE
)

# UMAP1
rf_umap1 <- randomForest(
  x         = X_num,
  y         = umap_coords$UMAP1,
  ntree     = 500,
  importance = TRUE
)

cat("\nRandom Forest R² (UMAP2):", round(1 - rf_umap2$mse[500] / var(umap_coords$UMAP2), 3), "\n")
cat("Random Forest R² (UMAP1):", round(1 - rf_umap1$mse[500] / var(umap_coords$UMAP1), 3), "\n")

# -- 4a. Extract importance and join section_lookup ------------------------------
imp_umap2 <- data.frame(
  variable  = rownames(importance(rf_umap2)),
  IncMSE    = importance(rf_umap2)[, "%IncMSE"],
  IncPurity = importance(rf_umap2)[, "IncNodePurity"]
) %>%
  left_join(section_lookup, by = "variable") %>%
  arrange(desc(IncMSE))

imp_umap1 <- data.frame(
  variable  = rownames(importance(rf_umap1)),
  IncMSE    = importance(rf_umap1)[, "%IncMSE"],
  IncPurity = importance(rf_umap1)[, "IncNodePurity"]
) %>%
  left_join(section_lookup, by = "variable") %>%
  arrange(desc(IncMSE))

cat("\nTop 20 variables by RF importance for UMAP2 (%IncMSE):\n")
print(head(imp_umap2[, c("variable", "section_label", "IncMSE")], 20))

cat("\nTop 20 variables by RF importance for UMAP1 (%IncMSE):\n")
print(head(imp_umap1[, c("variable", "section_label", "IncMSE")], 20))

# -- 4b. Plot UMAP2 ───────────────────────────────────────────────────────────
top_rf_umap2 <- head(imp_umap2, 20)

ggplot(top_rf_umap2, aes(x = reorder(variable, IncMSE),
                         y = IncMSE, fill = section_label)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = section_colors, name = "Section") +
  labs(title    = "Top-20 variables by RF importance for UMAP2",
       subtitle = "Metric: %IncMSE (permutation-based)",
       x = NULL, y = "%IncMSE") +
  theme_minimal(base_size = 10) +
  theme(legend.position  = "bottom",
        panel.grid.minor = element_blank())

ggsave(file.path(PLOT_DIR, "rf_importance_umap2.pdf"),
       width = 8, height = 7)

# -- 4c. Plot UMAP1 ───────────────────────────────────────────────────────────
top_rf_umap1 <- head(imp_umap1, 20)

ggplot(top_rf_umap1, aes(x = reorder(variable, IncMSE),
                         y = IncMSE, fill = section_label)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = section_colors, name = "Section") +
  labs(title    = "Top-20 variables by RF importance for UMAP1",
       subtitle = "Metric: %IncMSE (permutation-based)",
       x = NULL, y = "%IncMSE") +
  theme_minimal(base_size = 10) +
  theme(legend.position  = "bottom",
        panel.grid.minor = element_blank())

ggsave(file.path(PLOT_DIR, "rf_importance_umap1.pdf"),
       width = 8, height = 7)
# PCA vs UMAP ####
# ── Join PC loadings with section_lookup ─────────────────────────────────────

top_pc1_labeled <- pca_model$rotation %>%
  as.data.frame() %>%
  tibble::rownames_to_column("variable") %>%
  dplyr::mutate(abs_loading = abs(PC1)) %>%
  dplyr::arrange(desc(abs_loading)) %>%
  dplyr::select(variable, PC1) %>%
  head(20) %>%
  left_join(section_lookup, by = "variable")

top_pc2_labeled <- pca_model$rotation %>%
  as.data.frame() %>%
  tibble::rownames_to_column("variable") %>%
  dplyr::mutate(abs_loading = abs(PC2)) %>%
  dplyr::arrange(desc(abs_loading)) %>%
  dplyr::select(variable, PC2) %>%
  head(20) %>%
  left_join(section_lookup, by = "variable")

# ── Plots ─────────────────────────────────────────────────────────────────────
p_rf1 <- ggplot(head(imp_umap1, 20),
                aes(x = reorder(variable, IncMSE),
                    y = IncMSE, fill = section_label)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = section_colors, name = "Section") +
  labs(title    = "RF importance — UMAP1",
       subtitle = "Metric: %IncMSE",
       x = NULL, y = "%IncMSE") +
  theme_minimal(base_size = 10) +
  theme(legend.position  = "none",
        panel.grid.minor = element_blank())

p_pc1 <- ggplot(top_pc1_labeled,
                aes(x = reorder(variable, abs(PC1)),
                    y = PC1, fill = section_label)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = section_colors, name = "Section") +
  labs(title    = "PCA loadings — PC1",
       subtitle = "Sorted by |loading|",
       x = NULL, y = "Loading") +
  theme_minimal(base_size = 10) +
  theme(legend.position  = "none",
        panel.grid.minor = element_blank())

p_rf2 <- ggplot(head(imp_umap2, 20),
                aes(x = reorder(variable, IncMSE),
                    y = IncMSE, fill = section_label)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = section_colors, name = "Section") +
  labs(title    = "RF importance — UMAP2",
       subtitle = "Metric: %IncMSE",
       x = NULL, y = "%IncMSE") +
  theme_minimal(base_size = 10) +
  theme(legend.position  = "none",
        panel.grid.minor = element_blank())

p_pc2 <- ggplot(top_pc2_labeled,
                aes(x = reorder(variable, abs(PC2)),
                    y = PC2, fill = section_label)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = section_colors, name = "Section") +
  labs(title    = "PCA loadings — PC2",
       subtitle = "Sorted by |loading|",
       x = NULL, y = "Loading") +
  theme_minimal(base_size = 10) +
  theme(legend.position  = "none",
        panel.grid.minor = element_blank())

# -- Common legend -----------------------------------------------------------------
legend_plot <- ggplot(head(imp_umap1, 20),
                      aes(x = variable, y = IncMSE, fill = section_label)) +
  geom_col() +
  scale_fill_manual(values = section_colors, name = "Section") +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 3))

shared_legend <- cowplot::get_legend(legend_plot)

# -- Final layout ------------------------------------------------------------------

bottom_row    <- (p_rf1 | p_pc2) + plot_annotation(title = "UMAP1 vs PC2")
top_row <- (p_rf2 | p_pc1) + plot_annotation(title = "UMAP2 vs PC1")

final_plot <- cowplot::plot_grid(
  top_row, bottom_row, shared_legend,
  ncol    = 1,
  rel_heights = c(1, 1, 0.15)
)

ggsave(file.path(PLOT_DIR, "rf_vs_pca_loadings.pdf"),
       final_plot, width = 14, height = 16)




