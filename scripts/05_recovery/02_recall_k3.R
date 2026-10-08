# =============================================================================
# scripts/05_recovery/02_recall_k3.R
# Purpose : Recall is computed by 01_recovery_k3.R and saved as recall_k3.rds
#           in each block folder. This script only collects it and writes
#           the table and figure for the thesis. Run after 01_recovery_k3.R.
# Thesis  : Section 4.3.3; Table 4.5 (\label{tab:per-class}), Figure 4.8,
#           Table 4.8 (HDBSCAN recall/precision, console)
# Input   : results/objects/recovery_k3/<block>/recall_k3.rds
# Output  : results/objects/recovery_k3/{recall_k3_all,balanced_accuracy_k3}.rds
#           results/tables/recovery_k3/table_recall_k3.tex
#           results/figures/05_recovery/fig_recall_k3.pdf
# Note    : the short cluster names assume labels aligned to KM-PC5.
# =============================================================================

source(here::here("R", "00_setup.R"))

out_root <- file.path(DIR_OBJ, "recovery_k3")
tex_dir  <- ensure_dir(DIR_TAB, "recovery_k3")
plot_dir <- ensure_dir(DIR_FIG, "05_recovery")

BLOCKS    <- c("pol", "dem", "pol_dem")
BLOCK_LAB <- c(pol = "Political", dem = "Sociodemographic", pol_dem = "Both")

# ── 1. Collect ───────────────────────────────────────────────────────####
recall <- bind_rows(lapply(BLOCKS, function(b) {
  f <- file.path(out_root, b, "recall_k3.rds")
  if (!file.exists(f)) { warning("missing: ", f); return(NULL) }
  readRDS(f)
}))
stopifnot(nrow(recall) > 0)
saveRDS(recall, file.path(out_root, "recall_k3_all.rds"))

recall <- recall %>%
  mutate(block_lab = factor(BLOCK_LAB[block], levels = BLOCK_LAB),
         cluster   = as.character(cluster))

cat("\n──── Recall, LASSO, k = 3, full sample ────\n")
print(recall %>%
        filter(classifier == "LASSO", run_type == "FULL") %>%
        select(block_lab, method, cluster, support, recall, precision, f1) %>%
        arrange(block_lab, method, cluster), n = 60, digits = 3)

# Balanced accuracy and the worst-recovered cluster, per configuration.
bal <- recall %>%
  group_by(block_lab, run_type, method, classifier) %>%
  summarise(balanced_acc = mean(recall, na.rm = TRUE),
            min_recall   = min(recall,  na.rm = TRUE),
            worst        = cluster[which.min(recall)], .groups = "drop")
cat("\n──── Balanced accuracy / worst cluster ────\n")
print(bal %>% filter(classifier == "LASSO") %>% arrange(run_type, block_lab, method),
      n = 60, digits = 3)
saveRDS(bal, file.path(out_root, "balanced_accuracy_k3.rds"))

# ── 2. LaTeX table: recall by cluster and block, FULL runs ───────────####
# Cluster labels for the four k-means/PAM partitions are Hungarian-aligned
# to KM-PC5, so 1/2/3 mean the same thing across those four. HDBSCAN's own
# labels are NOT aligned with them and are excluded from this table.
CL <- c("1" = "Precarious and Disaffected",
        "2" = "Secular and Globally Oriented",
        "3" = "Religious and Order-Oriented")

tab <- recall %>%
  filter(classifier == "LASSO", run_type == "FULL") %>%
  select(block_lab, method, cluster, support, recall) %>%
  pivot_wider(names_from = block_lab, values_from = recall) %>%
  arrange(method, cluster)

lines <- c(
  "\\begin{table}[H]", "\\centering",
  paste0("\\caption{Per-cluster recall at $k=3$ on the full sample, under ",
         "LASSO-penalised multinomial regression. Recall is the share of a ",
         "cluster's true members correctly identified. Cluster labels are ",
         "Hungarian-aligned to KM-PC5, so they denote the same three positions ",
         "across configurations.}"),
  "\\label{tab:per-class}", "\\small",
  "\\begin{tabular}{llcccc}", "\\toprule",
  "Configuration & Cluster & $n$ & Political & Sociodem. & Both \\\\", "\\midrule")

prev <- NULL
for (i in seq_len(nrow(tab))) {
  r <- tab[i, ]
  if (!is.null(prev) && r$method != prev) lines <- c(lines, "\\addlinespace")
  lines <- c(lines, sprintf("%s & %s & %d & %.3f & %.3f & %.3f \\\\",
                            if (identical(r$method, prev)) "" else r$method,
                            CL[[r$cluster]], r$support,
                            r$Political, r$Sociodemographic, r$Both))
  prev <- r$method
}
lines <- c(lines, "\\bottomrule", "\\end{tabular}", "\\end{table}")
writeLines(lines, file.path(tex_dir, "table_recall_k3.tex"))
cat("\nLaTeX table written to", file.path(tex_dir, "table_recall_k3.tex"), "\n")

# ── 3. Figure ────────────────────────────────────────────────────────####
fig <- recall %>%
  filter(classifier == "LASSO", run_type == "FULL") %>%
  mutate(cluster = factor(cluster, levels = names(CL),
                          labels = gsub(" and ", "\nand ", CL))) %>%
  ggplot(aes(cluster, recall, fill = block_lab)) +
  geom_hline(yintercept = 1/3, linetype = 2, colour = "grey55", linewidth = 0.4) +
  geom_col(position = position_dodge(0.8), width = 0.72) +
  facet_wrap(~ method, nrow = 1) +
  scale_fill_manual(values = c(Political = "#1F3B57",
                               Sociodemographic = "#E15759",
                               Both = "#59A14F")) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  labs(x = NULL, y = "Recall", fill = NULL,
       title = "Which cluster is recovered, and from which block",
       subtitle = "k = 3, LASSO, out of sample. Dashed line: recall under uniform guessing.") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold", size = 9),
        legend.position = "bottom")
ggsave(file.path(plot_dir, "fig_recall_k3.pdf"), fig, width = 11, height = 4.4)
cat("Figure written to", file.path(plot_dir, "fig_recall_k3.pdf"), "\n")

# ── 4. HDBSCAN, reported separately ──────────────────────────────────####
# Its three clusters are not the same objects as the aligned 1/2/3 above,
# so recall is comparable within this table only.
cat("\n──── HDBSCAN recall (own labels, assigned respondents only) ────\n")
print(recall %>%
        filter(classifier == "LASSO", run_type == "HDB") %>%
        select(block_lab, cluster, support, recall, precision) %>%
        arrange(block_lab, cluster), n = 20, digits = 3)