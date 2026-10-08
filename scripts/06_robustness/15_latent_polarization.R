################################################################################
# scripts/06_robustness/15_latent_polarization.R
# Thesis  : not reported (Desmet-style latent value polarisation, robustness)
# Output  : results/objects/polarization/lvp_all_methods_vs_Q240.rds, figures
#
#     LVP(P) = ( dbar_all - dbar_within(P) ) / dbar_all
#
# Share of total value heterogeneity accounted for by a partition P.
# Computed for EVERY clustering method stored in all_clusters_aligned.rds
# (k <= K_MAX) and for the benchmark partitions built from Q240.
#
# Two distance spaces:
#   "gower" -> pre-computed Gower distance on the 198 value items (reference)
#   "pc2"   -> Euclidean distance on the first 2 PCs (robustness only: it is the
#              space KMeans_2pc optimises in, so it favours that method)
#
# NOTE: LVP is invariant to label permutation (aligned vs raw give the same value)
#       and increases mechanically with k -> compare only within a given k.
# NOTE: dataset fully imputed; no NA handling.
################################################################################

#### 00. Setup #################################################################

rm(list = intersect(ls(), c("wvs", "clusters", "dist_list", "partitions",
                            "lvp_results", "p_lvp", "D_gower", "D_pc2")))

source(here::here("R", "00_setup.R"))

obj_dir  <- DIR_OBJ
data_dir <- DIR_DATA
out_dir  <- ensure_dir(DIR_OBJ, "polarization")
plot_dir <- ensure_dir(DIR_FIG, "06_robustness", "polarization")

K_MAX <- 5   # maximum number of groups considered

#### 01. Data ##################################################################

wvs        <- readRDS(file.path(data_dir, "imputed_labeled.rds"))
clusters   <- readRDS(PATH_CLUSTERS)
gower_dist <- readRDS(file.path(DIR_CLUST, "gower_dist.rds"))
pca_scores <- readRDS(file.path(obj_dir,  "pca_scores.rds"))

n <- nrow(wvs)

if ("id" %in% names(wvs) && "id" %in% names(clusters)) {
  clusters <- clusters[match(wvs$id, clusters$id), , drop = FALSE]
  stopifnot(identical(as.character(wvs$id), as.character(clusters$id)))
} else {
  stopifnot(nrow(clusters) == n)
}

#### 02. Distance matrices #####################################################

D_gower <- as.matrix(gower_dist)

pc2   <- as.matrix(pca_scores)[, 1:2, drop = FALSE]
D_pc2 <- as.matrix(dist(pc2, method = "euclidean"))

dist_list <- list(gower = D_gower, pc2 = D_pc2)

for (nm in names(dist_list)) {
  stopifnot(is.matrix(dist_list[[nm]]),
            nrow(dist_list[[nm]]) == n, ncol(dist_list[[nm]]) == n,
            isTRUE(all.equal(unname(diag(dist_list[[nm]])), rep(0, n))))
}
rm(gower_dist, pca_scores)

#### 03. Core function #########################################################

lvp_from_dist <- function(D, g) {
  stopifnot(is.matrix(D), nrow(D) == length(g), !anyNA(g))
  
  g <- as.integer(as.factor(g))
  n <- length(g)
  
  total_sum   <- sum(D) / 2
  total_pairs <- n * (n - 1) / 2
  dbar_all    <- total_sum / total_pairs
  
  within_sum <- 0; within_pairs <- 0
  for (k in sort(unique(g))) {
    idx <- which(g == k); nk <- length(idx)
    if (nk < 2L) next
    within_sum   <- within_sum   + sum(D[idx, idx]) / 2
    within_pairs <- within_pairs + nk * (nk - 1) / 2
  }
  dbar_within  <- within_sum / within_pairs
  dbar_between <- (total_sum - within_sum) / (total_pairs - within_pairs)
  
  data.frame(
    K            = length(unique(g)),
    n            = n,
    dbar_all     = dbar_all,
    dbar_within  = dbar_within,
    dbar_between = dbar_between,
    LVP          = (dbar_all - dbar_within) / dbar_all,
    sep_ratio    = dbar_between / dbar_all,
    stringsAsFactors = FALSE
  )
}

#### 04. Partitions: ALL clustering methods ####################################

partitions <- list()

# --- (A) every method x every k found in the object ---------------------------
clu_cols <- setdiff(names(clusters), "id")
clu_cols <- grep("_k\\d+(_aligned)?$", clu_cols, value = TRUE)
stopifnot(length(clu_cols) > 0)

clu_info <- data.frame(
  col     = clu_cols,
  method  = sub("_k\\d+(_aligned)?$", "", clu_cols),
  k       = as.integer(sub("^.*_k(\\d+)(_aligned)?$", "\\1", clu_cols)),
  aligned = grepl("_aligned$", clu_cols),
  stringsAsFactors = FALSE
)
clu_info <- clu_info[clu_info$k <= K_MAX, , drop = FALSE]

# keep the aligned column when both versions exist for the same method x k
clu_info <- clu_info[order(clu_info$method, clu_info$k, -clu_info$aligned), ]
clu_info <- clu_info[!duplicated(paste(clu_info$method, clu_info$k)), ]

message("Partitions found:\n",
        paste(sprintf("  %-18s k = %s", unique(clu_info$method),
                      tapply(clu_info$k, clu_info$method,
                             function(x) paste(sort(x), collapse = ", "))[unique(clu_info$method)]),
              collapse = "\n"))

for (i in seq_len(nrow(clu_info))) {
  partitions[[clu_info$col[i]]] <- list(
    partition = clu_info$method[i],
    k_nominal = clu_info$k[i],
    scheme    = NA_character_,
    g         = clusters[[clu_info$col[i]]]
  )
}

# --- (B) benchmark partitions from Q240 ---------------------------------------
q240_num <- as.integer(as.character(wvs$Q240))
stopifnot(!anyNA(q240_num), all(q240_num %in% 1:10))

q240_schemes <- list(
  "2" = c(5, 10),            # 1-5 | 6-10
  "3" = c(3, 6, 10),         # 1-3 | 4-6 | 7-10
  "4" = c(3, 5, 7, 10),      # 1-3 | 4-5 | 6-7 | 8-10
  "5" = c(2, 4, 6, 8, 10)    # 1-2 | 3-4 | 5-6 | 7-8 | 9-10
)

cut_q240 <- function(x, upper) {
  brk <- c(0, upper)
  lab <- paste0(head(brk, -1) + 1, "-", upper)
  cut(x, breaks = brk, labels = lab, include.lowest = TRUE, right = TRUE)
}

for (kk in names(q240_schemes)) {
  g <- cut_q240(q240_num, q240_schemes[[kk]])
  stopifnot(!anyNA(g))
  partitions[[paste0("Q240_k", kk)]] <- list(
    partition = "Q240",
    k_nominal = as.integer(kk),
    scheme    = paste(levels(g), collapse = " | "),
    g         = g
  )
}

message("Total partitions to evaluate: ", length(partitions))

#### 05. Compute LVP: every partition x every distance space ###################

lvp_results <- do.call(rbind, lapply(names(dist_list), function(space) {
  D <- dist_list[[space]]
  do.call(rbind, lapply(partitions, function(p) {
    data.frame(dist_space = space,
               partition  = p$partition,
               k_nominal  = p$k_nominal,
               lvp_from_dist(D, p$g),
               scheme     = p$scheme,
               row.names  = NULL)
  }))
}))

lvp_results <- lvp_results[order(lvp_results$dist_space,
                                 lvp_results$k_nominal,
                                 lvp_results$partition), ]
rownames(lvp_results) <- NULL

print(lvp_results[, c("dist_space", "partition", "k_nominal", "K",
                      "dbar_within", "LVP", "sep_ratio")], digits = 4)

#### 06. Plot ##################################################################

# --- readable method labels ---------------------------------------------------
pretty_map <- c(
  "km_pc2"       = "k-means (2 PCs)",
  "km_pc5"       = "k-means (5 PCs)",
  "pam"    = "PAM (Gower)",
  "km_umap"  = "k-means (UMAP)",
  "Q240"         = "Q240 (left-right)"
)

pretty_lab <- function(x) ifelse(x %in% names(pretty_map),
                                 pretty_map[x], gsub("_", " ", x))

lvp_results$method_lab <- pretty_lab(lvp_results$partition)

# Q240 is the benchmark: always last in the legend, black and dashed.
meth_levels <- c(sort(setdiff(unique(lvp_results$method_lab), "Q240 (left-right)")),
                 "Q240 (left-right)")
lvp_results$method_lab <- factor(lvp_results$method_lab, levels = meth_levels)

n_meth <- length(meth_levels) - 1L
pal    <- setNames(c(hcl.colors(n_meth, palette = "Dark 3"), "black"), meth_levels)
ltys   <- setNames(c(rep("solid", n_meth), "22"), meth_levels)
shp    <- setNames(c(rep_len(c(16, 17, 15, 18, 8, 4, 3), n_meth), 21), meth_levels)
lwds   <- setNames(c(rep(0.6, n_meth), 0.9), meth_levels)

facet_lab <- c(gower = "Gower distance (198 value items)",
               pc2   = "Euclidean distance (first 2 PCs)")

p_lvp <- ggplot(lvp_results,
                aes(x = k_nominal, y = LVP,
                    colour = method_lab, shape = method_lab,
                    linetype = method_lab, linewidth = method_lab)) +
  geom_line() +
  geom_point(size = 2.3, fill = "white") +
  facet_wrap(~ dist_space, scales = "free_y",
             labeller = labeller(dist_space = facet_lab)) +
  scale_colour_manual(values = pal) +
  scale_linetype_manual(values = ltys) +
  scale_shape_manual(values = shp) +
  scale_linewidth_manual(values = lwds) +
  scale_x_continuous(breaks = sort(unique(lvp_results$k_nominal))) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Number of groups (k)",
       y = "Latent value-based polarisation",
       colour = "Partition", shape = "Partition",
       linetype = "Partition", linewidth = "Partition") +
  theme_minimal(base_size = 11) +
  theme(legend.position  = "bottom",
        legend.key.width = unit(1.4, "lines"),
        panel.grid.minor = element_blank(),
        strip.text       = element_text(face = "bold")) +
  guides(colour   = guide_legend(nrow = 2, byrow = TRUE),
         shape    = guide_legend(nrow = 2, byrow = TRUE),
         linetype = guide_legend(nrow = 2, byrow = TRUE),
         linewidth = guide_legend(nrow = 2, byrow = TRUE))

print(p_lvp)

# Reference-space-only version, for the main text
p_lvp_gower <- p_lvp %+% subset(lvp_results, dist_space == "gower")

#### 07. Save ##################################################################

saveRDS(lvp_results, file.path(out_dir, "lvp_all_methods_vs_Q240.rds"))
write.csv(lvp_results, file.path(out_dir, "lvp_all_methods_vs_Q240.csv"),
          row.names = FALSE)
ggsave(file.path(plot_dir, "lvp_all_methods_vs_Q240.pdf"),
       p_lvp, width = 8.5, height = 4.8)
ggsave(file.path(plot_dir, "lvp_all_methods_gower.pdf"),
       p_lvp_gower, width = 5.5, height = 4.8)