# =============================================================================
# scripts/06_robustness/13_pca_minmax_desmet.R
# Purpose : Comparison with Desmet et al. (2025), who rescale every item to
#           [0, 1] (min-max) instead of standardising it (z-score). The
#           thesis always uses the z-score PCA; this script only documents
#           how the min-max normalisation changes the picture:
#   1. Min-max PCA and its top-10 loadings on PC1 / PC2.
#   2. Density of the min-max PCA scores (HDR and "paper" rendering).
#   3. Diagnostics: variance explained, Procrustes distance to the z-score
#      configuration.
#   4. Five-dimensional (min-max) structure seen in two dimensions.
# Thesis  : Section 2.2.2 (normalisation); not a thesis figure
# Input   : data/processed/data_imputed_selected.rds
#           results/objects/pca_scores.rds (z-score PCA, for sign alignment)
# Output  : results/objects/{pca_minmax,pc5_knn_density}.rds
#           results/figures/06_robustness/pca_minmax/*.pdf
# =============================================================================

source(here::here("R", "00_setup.R"))

OBJ_DIR  <- DIR_OBJ
PLOT_DIR <- ensure_dir(DIR_FIG, "06_robustness", "pca_minmax")

#### 1. Min-max PCA on the item codes and top loadings ####

# Desmet et al. keep the political items and exclude only the demographics.
# Note: data_imputed_selected.rds already contains the value block only, so
# here only the value items enter.
full_data <- readRDS(file.path(DIR_DATA, "data_imputed_selected.rds"))

w_data <- full_data %>%
  dplyr::select(-dplyr::any_of(Demographics)) %>%
  mutate(across(where(is.factor), ~ as.numeric(as.character(.))))

cat("Number of questions used:", ncol(w_data), "\n")

# Min-max normalisation to [0, 1] (Desmet et al. style)
w_data_01 <- w_data %>%
  mutate(across(everything(), function(col) {
    rng <- range(col, na.rm = TRUE)
    if (diff(rng) == 0) return(rep(0, length(col)))  # constant -> 0
    (col - rng[1]) / (rng[2] - rng[1])
  }))

# PCA on the covariance matrix (centre yes, scale NO)
pca_min_max <- prcomp(w_data_01, center = TRUE, scale. = FALSE)

# variance explained: compare with their 16.4% for the first 2 PCs
var_exp <- pca_min_max$sdev^2 / sum(pca_min_max$sdev^2)
cat(sprintf("Var. explained PC1 = %.1f%%, PC2 = %.1f%%, PC1+PC2 = %.1f%%\n",
            var_exp[1]*100, var_exp[2]*100, sum(var_exp[1:2])*100))

loadings <- as.data.frame(pca_min_max$rotation[, 1:2])
loadings$variable <- rownames(loadings)

top_pc1 <- loadings[order(-abs(loadings$PC1)), ][1:10, c("variable", "PC1")]
top_pc2 <- loadings[order(-abs(loadings$PC2)), ][1:10, c("variable", "PC2")]

cat("\n=== TOP 10 loadings PC1 (min-max [0,1]) ===\n")
print(top_pc1, row.names = FALSE)

cat("\n=== TOP 10 loadings PC2 (min-max [0,1]) ===\n")
print(top_pc2, row.names = FALSE)

#### 2. Min-max PCA with dummy-coded nominal items ####

dat <- readRDS(file.path(DIR_DATA, "data_imputed_selected.rds"))

D <- dat
stopifnot(!anyNA(D))

# If TRUE, each question contributes equally regardless of how many dummy
# columns it generates. Document whichever you pick.
EQUAL_QUESTION_WEIGHT <- FALSE

minmax <- function(v) {
  r <- range(v)
  if (diff(r) == 0) return(rep(0, length(v)))   # constant column -> no signal
  (v - r[1]) / diff(r)
}

# Block-by-block construction. This avoids the Q182 / Q1821 regex collision
# entirely, because provenance is carried by the loop, not by the column name.
blocks <- vector("list", ncol(D))
names(blocks) <- names(D)

for (v in names(D)) {
  col <- D[[v]]
  if (is.numeric(col) || is.ordered(col)) {
    B <- matrix(minmax(as.numeric(col)), ncol = 1,
                dimnames = list(NULL, v))
  } else {
    col <- droplevels(as.factor(col))
    B <- stats::model.matrix(~ col - 1)          # already 0/1
    colnames(B) <- paste(v, levels(col), sep = "__")
  }
  if (EQUAL_QUESTION_WEIGHT) B <- B / sqrt(ncol(B))
  blocks[[v]] <- B
}

X_mm <- do.call(cbind, blocks)
stopifnot(is.matrix(X_mm), nrow(X_mm) == nrow(D),
          all(X_mm >= 0 - 1e-12), all(X_mm <= 1 + 1e-12))

# Centre but DO NOT scale: the [0,1] rescaling IS the normalisation.
pca_mm <- stats::prcomp(X_mm, center = TRUE, scale. = FALSE)

# PCA signs are arbitrary. Align to the z-score solution so the two figures
# are visually comparable rather than mirrored.
pca_scores_z <- readRDS(file.path(OBJ_DIR, "pca_scores.rds"))
scores_mm <- pca_mm$x[, 1:5, drop = FALSE]
for (j in 1:2) {
  if (stats::cor(scores_mm[, j], pca_scores_z[, j]) < 0)
    scores_mm[, j] <- -scores_mm[, j]
}

pts_mm <- data.frame(X = scores_mm[, 1], Y = scores_mm[, 2])
dd_mm  <- density_hdr(pts_mm$X, pts_mm$Y)

p_mm       <- plot_hdr(dd_mm, pts_mm, "PC1 (min-max)", "PC2 (min-max)", "hdr")
p_mm_paper <- plot_hdr(dd_mm, pts_mm, "PC1 (min-max)", "PC2 (min-max)", "paper")

print(p_mm); print(p_mm_paper)

ggsave(file.path(PLOT_DIR, "pca_minmax_density_hdr.pdf"),   p_mm,
       width = 7, height = 6, device = cairo_pdf)
ggsave(file.path(PLOT_DIR, "pca_minmax_density_paper.pdf"), p_mm_paper,
       width = 7, height = 6, device = cairo_pdf)

#### 3. Diagnostics: is the difference normalisation or rendering? ####

# (a) variance explained
ve_mm <- pca_mm$sdev^2 / sum(pca_mm$sdev^2)
cat(sprintf("min-max PC1 %.3f  PC2 %.3f  cum5 %.3f\n",
            ve_mm[1], ve_mm[2], sum(ve_mm[1:5])))

# (b) top loadings of PC1
head(sort(abs(pca_mm$rotation[, 1]), decreasing = TRUE), 20)

# (c) how similar are the two 2D configurations, up to rotation/reflection?
proc <- vegan::procrustes(pca_scores_z[, 1:2], scores_mm[, 1:2],
                          scale = TRUE, symmetric = TRUE)
cat(sprintf("Procrustes SS = %.4f  (0 = identical configuration)\n", proc$ss))

saveRDS(list(pca = pca_mm, scores = scores_mm,
             equal_weight = EQUAL_QUESTION_WEIGHT),
        file.path(OBJ_DIR, "pca_minmax.rds"))


#### 4. Five-dimensional structure, seen in two dimensions ####

rm(list = ls(pattern = "^(S5|dens5|pts5|p_5d|p_pairs)$"))

S5 <- scores_mm[, 1:5, drop = FALSE]
stopifnot(ncol(S5) == 5, !anyNA(S5))

# (a) kNN density estimate in the FULL 5D space, then shown on the PC1-PC2 plane.
d      <- 5L
k      <- max(10L, round(sqrt(nrow(S5))))
r_k    <- dbscan::kNNdist(S5, k = k)[, k]
vol_d  <- pi^(d / 2) / gamma(d / 2 + 1)
dens5  <- k / (nrow(S5) * vol_d * r_k^d)
stopifnot(all(is.finite(dens5)), all(dens5 > 0))

pts5 <- data.frame(X = S5[, 1], Y = S5[, 2], dens5 = dens5)

p_5d <- ggplot(pts5, aes(X, Y, colour = dens5)) +
  geom_point(size = 0.9, alpha = 0.85) +
  scale_colour_viridis_c(name = "5D density\n(kNN, log)",
                         trans = "log10", option = "magma") +
  coord_equal() + labs(x = "PC1 (min-max)", y = "PC2 (min-max)") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

print(p_5d)
ggsave(file.path(PLOT_DIR, "pc5_density_on_pc12.pdf"), p_5d,
       width = 7, height = 6, device = cairo_pdf)

# (b) all ten planes of the 5D subspace, each with its own HDR contours.
pairs_df <- do.call(rbind, lapply(utils::combn(5, 2, simplify = FALSE), function(ij) {
  i <- ij[1]; j <- ij[2]
  dd <- density_hdr(S5[, i], S5[, j])
  cbind(dd$grid,
        panel = paste0("PC", i, " vs PC", j),
        lev50 = dd$levels[["50%"]],
        lev90 = dd$levels[["90%"]])
}))

pts_pairs <- do.call(rbind, lapply(utils::combn(5, 2, simplify = FALSE), function(ij) {
  data.frame(X = S5[, ij[1]], Y = S5[, ij[2]],
             panel = paste0("PC", ij[1], " vs PC", ij[2]))
}))

p_pairs <- ggplot() +
  geom_point(data = pts_pairs, aes(X, Y),
             colour = "grey60", size = 0.3, alpha = 0.25) +
  geom_contour(data = pairs_df, aes(X, Y, z = z, group = panel),
               breaks = NULL, colour = "#2171B5", linewidth = 0.4, bins = 6) +
  facet_wrap(~ panel, scales = "free", ncol = 4) +
  coord_equal() + labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank())

print(p_pairs)
ggsave(file.path(PLOT_DIR, "pc5_pairwise_density.pdf"), p_pairs,
       width = 11, height = 7, device = cairo_pdf)

saveRDS(list(dens5 = dens5, k = k, d = d),
        file.path(OBJ_DIR, "pc5_knn_density.rds"))

