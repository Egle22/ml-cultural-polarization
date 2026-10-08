# =============================================================================
# scripts/03_structure/02_mds_gower.R
# Purpose : Classical MDS (principal coordinates) of the Gower dissimilarity,
#           distortion diagnostics, scree plot, and the relative density
#           surface on the first two MDS axes.
# Thesis  : Section 4.1.1; Figure 4.2 (MDS panel)
# Input   : results/objects/clusters/gower_dist.rds
# Output  : results/objects/mds_gower.rds
#           results/figures/03_structure/mds_gower_scree.pdf
#           results/figures/03_structure/mds_gower_density_surface.html
# =============================================================================

source(here::here("R", "00_setup.R"))

OBJ_DIR  <- DIR_OBJ
PLOT_DIR <- ensure_dir(DIR_FIG, "03_structure")

#### 1. Load the Gower distance ####

gower_dist <- readRDS(file.path(DIR_CLUST, "gower_dist.rds"))
if (!inherits(gower_dist, "dist")) gower_dist <- stats::as.dist(gower_dist)
stopifnot(inherits(gower_dist, "dist"))

n_obs <- attr(gower_dist, "Size")
cat(sprintf("n = %d  |  %s pairwise distances  |  range [%.4f, %.4f]\n",
            n_obs, format(length(gower_dist), big.mark = ","),
            min(gower_dist), max(gower_dist)))
stopifnot(n_obs > 0, !anyNA(gower_dist))

#### 2. Classical MDS (principal coordinates analysis) ####

# k = 5 so the same object serves both the 2D figure and any 5D diagnostic.
# eig = TRUE is essential: Gower dissimilarity is not Euclidean, so the
# eigenvalue spectrum is what tells you how much the projection distorts.
K_MDS <- 5L
mds_obj <- stats::cmdscale(gower_dist, k = K_MDS, eig = TRUE)

mds_pts <- mds_obj$points
stopifnot(is.matrix(mds_pts), nrow(mds_pts) == n_obs,
          ncol(mds_pts) == K_MDS, !anyNA(mds_pts))

#### 3. Diagnostics: how badly does the Euclidean embedding distort? ####

ev <- mds_obj$eig

n_neg     <- sum(ev < 0)
neg_ratio <- if (n_neg > 0) abs(min(ev)) / max(ev) else 0

# Two conventions. Mardia's |lambda| criterion is optimistic because it treats
# negative eigenvalues as if they carried information; the positive-only
# version is the conservative one. Report both.
gof_abs <- cumsum(abs(ev))[1:2] / sum(abs(ev))
ev_pos  <- ev[ev > 0]
gof_pos <- cumsum(ev_pos)[1:2] / sum(ev_pos)

cat(sprintf("Negative eigenvalues: %d  |  |min neg| / max pos = %.3f\n",
            n_neg, neg_ratio))
cat(sprintf("Axes 1-2, |lambda| criterion : %.1f%%\n", gof_abs[2] * 100))
cat(sprintf("Axes 1-2, positive-only      : %.1f%%\n", gof_pos[2] * 100))
if (neg_ratio > 0.2)
  message("Non-Euclidean distortion is substantial; run the non-metric check ",
          "and report it in the appendix.")

# Direct fidelity check: correlation between original and embedded distances.
d_emb2 <- stats::dist(mds_pts[, 1:2])
d_emb5 <- stats::dist(mds_pts[, 1:5])
cat(sprintf("cor(Gower, MDS-2D) = %.3f   cor(Gower, MDS-5D) = %.3f\n",
            stats::cor(gower_dist, d_emb2),
            stats::cor(gower_dist, d_emb5)))

#### 4. Scree plot and save ####

scree <- data.frame(axis = 1:min(20, length(ev)),
                    eig  = ev[1:min(20, length(ev))])

p_scree <- ggplot(scree, aes(axis, eig)) +
  geom_col(fill = "#2171B5") +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  labs(x = "MDS axis", y = "Eigenvalue") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

print(p_scree)
ggsave(file.path(PLOT_DIR, "mds_gower_scree.pdf"), p_scree,
       width = 6, height = 4, device = cairo_pdf)

saveRDS(mds_obj, file.path(OBJ_DIR, "mds_gower.rds"))
message("Saved: ", file.path(OBJ_DIR, "mds_gower.rds"))


#### 5. Density surface on the MDS plane ####

pts_mds <- data.frame(X = mds_pts[, 1], Y = mds_pts[, 2])
dd_mds  <- density_hdr(pts_mds$X, pts_mds$Y)
kde_mds <- dd_mds$kde

# Relative height, peak = 100. Removes the units question and shows only
# whether a second local maximum exists and how tall it is relative to the first.
z_rel_mds <- t(kde_mds$z) / max(kde_mds$z) * 100

surf_mds <- plotly::plot_ly(x = kde_mds$x, y = kde_mds$y, z = z_rel_mds) |>
  plotly::add_surface(
    colorscale = "Viridis",
    contours = list(z = list(show = TRUE, usecolormap = TRUE,
                             project = list(z = TRUE))),
    colorbar = list(title = "Relative\ndensity\n(peak = 100)")) |>
  plotly::layout(
    title = "Classical MDS on Gower dissimilarity — density surface",
    scene = list(
      xaxis = list(title = "MDS 1"),
      yaxis = list(title = "MDS 2"),
      zaxis = list(title = "Relative density (peak = 100)"),
      aspectmode  = "manual",
      aspectratio = list(x = 1, y = 1, z = 0.6),
      camera = list(eye = list(x = 1.6, y = -1.6, z = 0.9))))

surf_mds

htmlwidgets::saveWidget(
  plotly::partial_bundle(surf_mds),
  file.path(PLOT_DIR, "mds_gower_density_surface.html"),
  selfcontained = TRUE)
