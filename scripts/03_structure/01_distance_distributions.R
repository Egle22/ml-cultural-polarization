# =============================================================================
# scripts/03_structure/01_distance_distributions.R
# Purpose : Distributions of pairwise distances between respondents under the
#           three representations: Euclidean distance in PCA space, Euclidean
#           distance on the 2D UMAP embedding, Gower distance on the full
#           data. All distances are normalised to [0, 1].
# Thesis  : Section 4.1.1; Figure 4.1
# Input   : results/objects/pca_scores.rds
#           results/objects/umapd_2d.rds
#           results/objects/clusters/gower_dist.rds
# Output  : results/figures/03_structure/distance_distributions{,_horizontal}.pdf
# Note    : KNOWN LIMITATION (not fixed, see README): the panel titled
#           "Euclidean distance on 2 PCs" is computed with dist(pca_scores) on
#           all 198 principal components, i.e. it is the full-space Euclidean
#           distance.
# =============================================================================

source(here::here("R", "00_setup.R"))

plot_dir <- ensure_dir(DIR_FIG, "03_structure")

#### 1. Data ####

gower_dist <- readRDS(file.path(DIR_CLUST, "gower_dist.rds"))
pca_scores <- readRDS(file.path(DIR_OBJ, "pca_scores.rds"))
umapd_2d   <- readRDS(file.path(DIR_OBJ, "umapd_2d.rds"))

#### 2. Normalised pairwise distances ####

D <- dist(pca_scores, method = "euclidean")
D <- D / max(D, na.rm = TRUE)
# D: either a 'dist' object or a full square matrix of pairwise distances
d_vec <- if (inherits(D, "dist")) as.numeric(D) else D[lower.tri(D)]

stopifnot(length(d_vec) == choose(attr(as.dist(D), "Size"), 2), !anyNA(d_vec))

g_gower     <- as.numeric(gower_dist)
d_euclidean <- as.numeric(D)
u_euclidean <- as.numeric(dist(umapd_2d, method = "euclidean"))
u_euclidean <- u_euclidean / max(u_euclidean)

# Remove invalid values
d_euclidean <- d_euclidean[is.finite(d_euclidean)]

#### 3. Histograms ####

# One histogram with its mean as a red vertical line.
hist_with_mean <- function(x, main, xlab) {
  hist(x, breaks = 65, main = main,
       xlab = xlab, ylab = "Frequency",
       xlim = c(0, 1), col = "grey80", border = "black")
  abline(v = mean(x, na.rm = TRUE), col = "red", lwd = 2)
}

# The three panels with a common legend, stacked vertically or side by side.
draw_distance_panels <- function(horizontal = FALSE) {
  if (horizontal) {
    layout(matrix(c(1, 1, 1,
                    2, 3, 4), nrow = 2, byrow = TRUE),
           heights = c(1, 6))
  } else {
    layout(matrix(c(1, 2, 3, 4), ncol = 1), heights = c(1, 6, 6, 6))
  }

  # common legend
  par(mar = c(0, 0, 0, 0))
  plot.new()
  legend("center", legend = "Mean", col = "red", lwd = 2, lty = 1,
         bty = "n", horiz = horizontal)

  par(mar = if (horizontal) c(4.5, 4.5, 3, 1) else c(4, 5, 3, 2))
  hist_with_mean(d_euclidean, "Euclidean distance on 2 PCs",
                 "Normalized Euclidean distance")
  hist_with_mean(u_euclidean, "Euclidean distance on UMAP embedding",
                 "Normalized Euclidean distance")
  hist_with_mean(g_gower, "Gower distance on full dataset",
                 "Gower distance")

  # restore the standard layout
  layout(1)
}

# Vertical version
draw_distance_panels(horizontal = FALSE)
pdf(file.path(plot_dir, "distance_distributions.pdf"), width = 7, height = 12)
draw_distance_panels(horizontal = FALSE)
dev.off()

# Horizontal version
draw_distance_panels(horizontal = TRUE)
pdf(file.path(plot_dir, "distance_distributions_horizontal.pdf"), width = 14, height = 5)
draw_distance_panels(horizontal = TRUE)
dev.off()
