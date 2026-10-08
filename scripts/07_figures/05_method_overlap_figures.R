## =============================================================================
## scripts/07_figures/05_method_overlap_figures.R
##   FIG 1 - between-method agreement as a function of k (without Cramer's V)
##   FIG 2 - cluster shares by partition, with and without HDBSCAN
##   FIG 3 - cross-method overlap of the LASSO variables (+ code to generate it)
## Thesis  : Tables 4.2 / A.3 / A.7 shown as figures (supporting material)
## Input   : results/objects/clusters/all_clusters_aligned.rds
##           (FIG 1 uses the agreement values of the concordance tables; FIG 3
##           falls back to the k = 3 vote counts of Table A.7)
## Output  : results/figures/07_figures/method_overlap/*.pdf
## =============================================================================

source(here::here("R", "00_setup.R"))
library(tibble)

#### 0. Paths and constants ####################################################

path_obj      <- DIR_OBJ
path_clusters <- DIR_CLUST
path_plots    <- ensure_dir(DIR_FIG, "07_figures", "method_overlap")

NOISE_COLOR <- "grey75"   # for HDBSCAN noise, when present as cluster 0

theme_thesis <- theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.3),
    strip.text         = element_text(face = "bold", size = 9),
    plot.title         = element_text(face = "bold", size = 12),
    plot.subtitle      = element_text(colour = "grey35", size = 9),
    plot.caption       = element_text(colour = "grey45", size = 7.5, hjust = 0),
    legend.position    = "bottom",
    legend.key.height  = unit(0.4, "cm")
  )
theme_set(theme_thesis)


#### 1. FIG 1 - Between-method agreement =======================================
## Values from the tables tab:concordance_low / tab:concordance_high.
## NB: "agreement", not "accuracy": there is no ground truth.

agr <- tribble(
  ~pair,          ~k, ~agreement,
  "PC2 vs PC5",     2, 99.1,  "PC2 vs PC5",     3, 91.8,
  "PC2 vs PC5",     4, 62.3,  "PC2 vs PC5",     5, 71.5,
  "PC2 vs PC5",     6, 62.7,  "PC2 vs PC5",     7, 52.5,
  "PC2 vs PC5",     8, 55.7,  "PC2 vs PC5",    16, 45.5,

  "UMAP vs PC5",    2, 92.3,  "UMAP vs PC5",    3, 85.4,
  "UMAP vs PC5",    4, 58.0,  "UMAP vs PC5",    5, 73.9,
  "UMAP vs PC5",    6, 66.0,  "UMAP vs PC5",    7, 56.2,
  "UMAP vs PC5",    8, 49.7,  "UMAP vs PC5",   16, 42.4,

  "PC2 vs PAM",     2, 86.6,  "PC2 vs PAM",     3, 68.4,
  "PC2 vs PAM",     4, 64.1,  "PC2 vs PAM",     5, 51.8,
  "PC2 vs PAM",     6, 54.6,  "PC2 vs PAM",     7, 47.1,
  "PC2 vs PAM",     8, 41.2,  "PC2 vs PAM",    16, 24.8,

  "PC5 vs PAM",     2, 86.7,  "PC5 vs PAM",     3, 68.5,
  "PC5 vs PAM",     4, 55.7,  "PC5 vs PAM",     5, 53.8,
  "PC5 vs PAM",     6, 45.6,  "PC5 vs PAM",     7, 44.0,
  "PC5 vs PAM",     8, 45.7,  "PC5 vs PAM",    16, 27.1
)

K_MAX    <- 6                      # set to 16 for the full version
agr_plot <- agr %>% filter(k <= K_MAX)

## Baseline: expected agreement between two independent balanced partitions = 1/k
baseline <- tibble(k = sort(unique(agr_plot$k))) %>% mutate(exp_agr = 100 / k)

## Labels at the end of the lines instead of a legend: more readable with 4 series
lab_end <- agr_plot %>% filter(k == max(k))

pair_cols <- c("PC2 vs PC5"  = "#4E79A7",
               "UMAP vs PC5" = "#59A14F",
               "PC2 vs PAM"  = "#E15759",
               "PC5 vs PAM"  = "#F28E2B")

fig1 <- ggplot(agr_plot, aes(k, agreement, colour = pair)) +
  geom_line(data = baseline, aes(k, exp_agr), inherit.aes = FALSE,
            linetype = 2, colour = "grey55", linewidth = 0.4) +
  annotate("text", x = K_MAX, y = 100 / K_MAX - 4,
           label = "chance agreement (1/k)",
           hjust = 1, size = 2.7, colour = "grey45") +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.1) +
  geom_text(data = lab_end, aes(label = pair),
            hjust = -0.12, size = 3, show.legend = FALSE) +
  scale_colour_manual(values = pair_cols) +
  scale_x_continuous(breaks = sort(unique(agr_plot$k)),
                     expand = expansion(mult = c(0.03, 0.30))) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25),
                     labels = function(x) paste0(x, "%")) +
  guides(colour = "none") +
  labs(x = "Number of clusters (k)",
       y = "Observations with the same label",
       title = "Between-method agreement declines as k grows",
       subtitle = "Labels aligned with the Hungarian algorithm",
       caption = paste("Agreement is downward-biased for k >= 4 because label",
                       "alignment is imperfect; the ordering across k is",
                       "nonetheless informative."))

print(fig1)
ggsave(file.path(path_plots, "fig_agreement_vs_k.pdf"), fig1,
       width = 6.5, height = 4.2)


#### 2. FIG 2 - Cluster shares =================================================
## Computed directly from all_clusters_aligned.rds: no numbers typed by hand.

all_cl <- readRDS(PATH_CLUSTERS)
n_tot  <- nrow(all_cl)

## Column-name parser: handles km_pc2_k3, pam_3, km_umap_5, hdb_umap2_mp25 ...
parse_partition <- function(nm) {
  is_hdb <- grepl("^hdb", nm)
  kv <- suppressWarnings(as.integer(sub(".*?_k?(\\d+)$", "\\1", nm)))
  method <- dplyr::case_when(
    grepl("^km_pc2",  nm) ~ "KM-PC2",
    grepl("^km_pc5",  nm) ~ "KM-PC5",
    grepl("^km_umap", nm) ~ "KM-UMAP",
    grepl("^pam",     nm) ~ "PAM",
    is_hdb                ~ "HDBSCAN",
    TRUE                  ~ NA_character_
  )
  tibble(partition = nm, method = method, k = kv, is_hdb = is_hdb)
}

## Partitions to show
keep <- grep("^(km_pc2|km_pc5|km_umap|pam)_k?[2-5]$|^hdb_", names(all_cl),
             value = TRUE)
stopifnot(length(keep) > 0)

shares <- lapply(keep, function(nm) {
  v  <- all_cl[[nm]]
  tb <- as.data.frame(table(cluster = v), stringsAsFactors = FALSE)
  meta <- parse_partition(nm)
  ## for HDBSCAN, k = number of non-noise clusters
  if (meta$is_hdb) meta$k <- sum(as.character(tb$cluster) != "0")
  tibble(
    partition = nm,
    method    = meta$method,
    k         = meta$k,
    cluster   = as.character(tb$cluster),
    n         = tb$Freq
  )
}) %>%
  bind_rows() %>%
  group_by(partition) %>%
  mutate(share = n / sum(n)) %>%      # over the assigned n
  ungroup() %>%
  mutate(
    method  = factor(method, levels = c("KM-PC2", "KM-PC5", "KM-UMAP",
                                        "PAM", "HDBSCAN")),
    cluster = factor(cluster, levels = as.character(0:16))
  )

## For HDBSCAN also compute the noise share on the FULL n (n = 2579), not on
## the assigned subset only.
noise_share <- shares %>%
  filter(method == "HDBSCAN", cluster == "0") %>%
  mutate(noise_pct = 100 * n / n_tot) %>%
  dplyr::select(partition, noise_pct)
print(noise_share)

## Palette: CLUSTER_COLORS + grey for the noise
pal <- c("0" = NOISE_COLOR, CLUSTER_COLORS)

make_shares_plot <- function(dat, title, subtitle = NULL) {
  ggplot(dat, aes(x = factor(k), y = share, fill = cluster)) +
    geom_col(width = 0.72, colour = "white", linewidth = 0.4) +
    geom_text(aes(label = ifelse(share >= 0.05,
                                 sprintf("%.0f%%", 100 * share), "")),
              position = position_stack(vjust = 0.5),
              size = 2.5, colour = "white", fontface = "bold") +
    facet_grid(~ method, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = pal, drop = TRUE,
                      breaks = as.character(1:6)) +
    scale_y_continuous(labels = function(x) paste0(100 * x, "%"),
                       expand = expansion(mult = c(0, 0.02))) +
    labs(x = "Number of clusters (k)", y = "Share of respondents",
         fill = "Cluster", title = title, subtitle = subtitle) +
    theme(panel.grid.major.y = element_blank())
}

## --- 2a: without HDBSCAN ------------------------------------------------------
fig2a <- make_shares_plot(
  filter(shares, method != "HDBSCAN"),
  title    = "Cluster shares across partitions",
  subtitle = "Partition-based methods split the sample into comparable blocks"
)
print(fig2a)
ggsave(file.path(path_plots, "fig_shares_nohdb.pdf"), fig2a,
       width = 8, height = 4)

## --- 2b: with HDBSCAN (noise in grey) -----------------------------------------
fig2b <- make_shares_plot(
  shares,
  title    = "Cluster shares across partitions, including HDBSCAN",
  subtitle = "Grey: HDBSCAN noise points (shares computed on assigned units)"
) +
  labs(caption = paste0(
    "HDBSCAN noise as a share of the full sample (n = ", n_tot, "): ",
    paste(sprintf("%s = %.1f%%", noise_share$partition,
                  noise_share$noise_pct), collapse = "; "), "."))
print(fig2b)
ggsave(file.path(path_plots, "fig_shares_hdb.pdf"), fig2b,
       width = 9, height = 4.2)


#### 3. FIG 3 - Cross-method overlap of the LASSO variables ====================
## 3.1 DATA GENERATION ---------------------------------------------------------
## The vote table exists only for k = 3 (04_profiles/03_method_overlap.R). This
## block regenerates it for any k with a simpler all-dummy one-vs-rest LASSO
## (5 folds, deviance). Run it once, save, then move to the plots.

if (FALSE) {

  imputed_labeled <- readRDS(file.path(DIR_DATA, "imputed_labeled.rds"))
  vars_analysis   <- setdiff(names(imputed_labeled), c(Political, Demographics))
  all_cl          <- readRDS(PATH_CLUSTERS)
  stopifnot(nrow(imputed_labeled) == nrow(all_cl))

  X_df <- imputed_labeled[, vars_analysis, drop = FALSE]

  ## CRITICAL: UNORDERED factors, otherwise model.matrix produces polynomial
  ## contrasts and the LASSO becomes uninterpretable.
  X_df <- dplyr::mutate(X_df, dplyr::across(dplyr::everything(),
                                            ~ factor(as.character(.x))))
  stopifnot(!any(vapply(X_df, is.ordered, logical(1))))

  mm_full <- model.matrix(~ ., data = X_df)
  asg     <- attr(mm_full, "assign")
  lbls    <- attr(terms(~ ., data = X_df), "term.labels")
  ## column -> original variable via `assign` (NOT via regex on the names:
  ## Q1821 vs Q18210 would be ambiguous)
  col_var <- c("(Intercept)", lbls)[asg + 1L]
  keep_c  <- col_var != "(Intercept)"
  X       <- mm_full[, keep_c, drop = FALSE]
  col_var <- col_var[keep_c]

  methods_map <- c(PC2 = "km_pc2", PC5 = "km_pc5",
                   UMAP = "km_umap", PAM = "pam")
  KS <- 2:6
  TOP_N <- 20        # as in the thesis: the 20 strongest regressors

  sel <- list()
  for (kv in KS) {
    sel[[as.character(kv)]] <- list()
    for (cl in seq_len(kv)) {
      per_method <- list()
      for (mn in names(methods_map)) {
        cand <- grep(paste0("^", methods_map[[mn]], "_k?", kv, "$"),
                     names(all_cl), value = TRUE)
        if (!length(cand)) next
        y <- as.integer(all_cl[[cand[1]]] == cl)
        if (length(unique(y)) < 2) next

        set.seed(42)
        fit <- cv.glmnet(X, y, family = "binomial", alpha = 1,
                         type.measure = "deviance", nfolds = 5)
        co  <- coef(fit, s = "lambda.1se")
        nz  <- which(as.numeric(co)[-1] != 0)          # intercept removed
        if (!length(nz)) next
        ord <- nz[order(abs(as.numeric(co)[-1][nz]), decreasing = TRUE)]
        per_method[[mn]] <- unique(col_var[ord])[seq_len(min(TOP_N, length(ord)))]
      }
      sel[[as.character(kv)]][[paste0("Cluster ", cl)]] <- per_method
    }
    message("k = ", kv, " done")
  }
  saveRDS(sel, file.path(path_obj, "lasso_selection_by_k.rds"))
}

## 3.2 VOTE COUNTS -------------------------------------------------------------

vote_counts <- function(method_lists) {
  method_lists <- method_lists[lengths(method_lists) > 0]
  if (!length(method_lists)) return(integer(0))
  tb <- table(unlist(lapply(method_lists, unique), use.names = FALSE))
  vapply(seq_along(method_lists), function(m) sum(tb == m), integer(1))
}

overlap_table <- function(sel) {
  do.call(rbind, lapply(names(sel), function(kk) {
    do.call(rbind, lapply(names(sel[[kk]]), function(cl) {
      cnt <- vote_counts(sel[[kk]][[cl]])
      if (!length(cnt)) return(NULL)
      data.frame(k = as.integer(kk), cluster = cl,
                 m = seq_along(cnt), n_vars = cnt,
                 stringsAsFactors = FALSE)
    }))
  }))
}

## If block 3.1 has been run, use instead:
##   ov <- as_tibble(overlap_table(readRDS(file.path(path_obj, "lasso_selection_by_k.rds"))))

## Default: k = 3 only, from the table tab:top-overlap-votes-k3 (Table A.7)
ov <- tribble(
  ~k, ~cluster,     ~m, ~n_vars,
  3, "Cluster 1",   4,  1,   3, "Cluster 1",   3,  8,
  3, "Cluster 1",   2, 12,   3, "Cluster 1",   1, 28,
  3, "Cluster 2",   4,  2,   3, "Cluster 2",   3,  9,
  3, "Cluster 2",   2,  8,   3, "Cluster 2",   1, 29,
  3, "Cluster 3",   4,  5,   3, "Cluster 3",   3,  4,
  3, "Cluster 3",   2,  9,   3, "Cluster 3",   1, 30
)

## 3.3 PLOTS -------------------------------------------------------------------

vote_pal <- c("all 4 methods" = "#2C3E50", "3 methods" = "#5D7A94",
              "2 methods" = "#9FB3C4", "1 method only" = "grey85")

fig3a <- ov %>%
  mutate(m = factor(m, levels = 4:1, labels = names(vote_pal))) %>%
  ggplot(aes(x = cluster, y = n_vars, fill = m)) +
  geom_col(width = 0.68, colour = "white", linewidth = 0.4) +
  facet_wrap(~ k, labeller = label_both, scales = "free_x", nrow = 1) +
  scale_fill_manual(values = vote_pal) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "Variables selected", fill = NULL,
       title = "Methods rarely agree on which variables matter",
       subtitle = "One-vs-rest LASSO, 20 strongest regressors per method")
print(fig3a)
ggsave(file.path(path_plots, "fig_lasso_votes.pdf"), fig3a,
       width = 7, height = 4)

fig3b <- ov %>%
  group_by(k, cluster) %>%
  summarise(`>= 2 methods` = sum(n_vars[m >= 2]),
            `all 4 methods` = sum(n_vars[m == 4]), .groups = "drop") %>%
  group_by(k) %>%
  summarise(across(everything(), mean), .groups = "drop") %>%
  dplyr::select(-cluster) %>%
  pivot_longer(-k, names_to = "level", values_to = "n_vars") %>%
  ggplot(aes(k, n_vars, colour = level, shape = level)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.4) +
  scale_colour_manual(values = c("all 4 methods" = "#2C3E50",
                                 ">= 2 methods"  = "#9FB3C4")) +
  scale_x_continuous(breaks = sort(unique(ov$k))) +
  labs(x = "Number of clusters (k)", y = "Variables per cluster (mean)",
       colour = NULL, shape = NULL,
       title = "Cross-method overlap of selected variables")
print(fig3b)
