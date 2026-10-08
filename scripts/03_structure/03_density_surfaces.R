# =============================================================================
# scripts/03_structure/03_density_surfaces.R
# Purpose : Kernel density surfaces (with highest-density-region contours)
#           over two-dimensional representations of the respondents:
#           z-score PCA (PC1-PC2) and the 2D UMAP embedding. The MDS panel is
#           built in 02_mds_gower.R. The min-max PCA used to compare with
#           Desmet et al. is in 06_robustness/13_pca_minmax_desmet.R.
# Thesis  : Section 4.1.1; Figure 4.2
# Input   : results/objects/pca_scores.rds (z-score PCA)
#           results/objects/UMAP/umap_2d.rds
# Output  : results/objects/pca_density_2d.rds
#           results/figures/03_structure/{pca_density_hdr_contours,
#           pca_density_filled,umap2d_density_hdr}.pdf
#           results/figures/03_structure/umap2d_density_surface.html
# Note    : deterministic (no sampling); bandwidths from the normal-reference rule.
# =============================================================================

source(here::here("R", "00_setup.R"))

OBJ_DIR  <- DIR_OBJ
PLOT_DIR <- ensure_dir(DIR_FIG, "03_structure")

#### 1. Data ####

pca_scores <- readRDS(file.path(OBJ_DIR, "pca_scores.rds"))
pc2 <- as.data.frame(pca_scores)[, 1:2]
names(pc2) <- c("PC1", "PC2")
stopifnot(nrow(pc2) > 0, all(vapply(pc2, is.numeric, logical(1))),
          !anyNA(pc2))

#### 2. Kernel density estimate on a grid ####

# Bandwidth: normal-reference rule per axis (Silverman). Kept explicit so the
# choice is documentable in the thesis appendix rather than hidden in a default.
h_vec <- c(MASS::bandwidth.nrd(pc2$PC1), MASS::bandwidth.nrd(pc2$PC2))

# Grid extended 10% beyond the data range so outer contours are not clipped.
pad <- function(v, f = 0.10) range(v) + c(-1, 1) * diff(range(v)) * f
kde <- MASS::kde2d(pc2$PC1, pc2$PC2,
                   h = h_vec, n = 400,
                   lims = c(pad(pc2$PC1), pad(pc2$PC2)))

dens_grid <- expand.grid(PC1 = kde$x, PC2 = kde$y)
dens_grid$z <- as.vector(kde$z)

#### 3. HDR levels: contours enclosing a given share of respondents ####

# Density value at each observation (nearest grid node), then the (1 - p)
# quantile of those values is the level enclosing p of the sample.
ix <- findInterval(pc2$PC1, kde$x, all.inside = TRUE)
iy <- findInterval(pc2$PC2, kde$y, all.inside = TRUE)
dens_at_point <- kde$z[cbind(ix, iy)]

hdr_probs  <- c(0.25, 0.50, 0.75, 0.90)
hdr_levels <- stats::quantile(dens_at_point, probs = 1 - hdr_probs, names = FALSE)
names(hdr_levels) <- paste0(hdr_probs * 100, "%")
stopifnot(!any(duplicated(round(hdr_levels, 12))))

# Levels are decreasing in density; factor() sorts increasing, so labels must
# be reversed to match the plotting order.
hdr_labels <- rev(names(hdr_levels))
hdr_cols   <- rev(c("#08306B", "#2171B5", "#6BAED6", "#C6DBEF"))

#### 4. Contour plot over the scatter ####

p_contour <- ggplot() +
  geom_point(data = pc2, aes(PC1, PC2),
             colour = "grey55", size = 0.55, alpha = 0.30) +
  geom_contour(data = dens_grid,
               aes(PC1, PC2, z = z, colour = factor(after_stat(level))),
               breaks = sort(hdr_levels), linewidth = 0.75) +
  scale_colour_manual(name = "Share of\nrespondents",
                      values = stats::setNames(hdr_cols, sort(hdr_levels)),
                      labels = hdr_labels) +
  coord_equal() +
  labs(x = "PC1", y = "PC2") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        legend.position  = "right")

print(p_contour)

ggsave(file.path(PLOT_DIR, "pca_density_hdr_contours.pdf"),
       p_contour, width = 7, height = 6, device = cairo_pdf)

#### 5. Filled version (the 3D surface seen from above) ####

p_filled <- ggplot(dens_grid, aes(PC1, PC2)) +
  geom_raster(aes(fill = z), interpolate = TRUE) +
  geom_contour(aes(z = z), breaks = sort(hdr_levels),
               colour = "white", linewidth = 0.4) +
  scale_fill_viridis_c(name = "Density", option = "magma") +
  coord_equal(expand = FALSE) +
  labs(x = "PC1", y = "PC2") +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank())

print(p_filled)

ggsave(file.path(PLOT_DIR, "pca_density_filled.pdf"),
       p_filled, width = 7, height = 6, device = cairo_pdf)

#### 6. Save objects ####

saveRDS(list(kde = kde, hdr_levels = hdr_levels, bandwidth = h_vec),
        file.path(OBJ_DIR, "pca_density_2d.rds"))

plotly::plot_ly(x = kde$x, y = kde$y, z = t(kde$z)) |>
  plotly::add_surface(colorscale = "Viridis") |>
  plotly::layout(scene = list(xaxis = list(title = "PC1"),
                              yaxis = list(title = "PC2"),
                              zaxis = list(title = "Density")))


#### 7. UMAP 2D density ####

rm(list = ls(pattern = "^(pts_um|dd_um|p_um)$"))

umap_2d <- readRDS(file.path(OBJ_DIR, "UMAP", "umap_2d.rds"))
emb <- if (is.list(umap_2d) && !is.null(umap_2d$layout)) umap_2d$layout else umap_2d
emb <- as.matrix(emb)
stopifnot(ncol(emb) == 2, nrow(emb) == nrow(pc2), !anyNA(emb))

pts_um <- data.frame(X = emb[, 1], Y = emb[, 2])
dd_um  <- density_hdr(pts_um$X, pts_um$Y)

p_um <- plot_hdr(dd_um, pts_um, "UMAP 1", "UMAP 2", "hdr")
print(p_um)

ggsave(file.path(PLOT_DIR, "umap2d_density_hdr.pdf"), p_um,
       width = 7, height = 6, device = cairo_pdf)

#### 8. UMAP 2D — interactive 3D density surface ####

rm(list = ls(pattern = "^(kde_um|z_rel_um|z_cnt_um|surf_um)$"))

stopifnot(exists("dd_um"))   # from block 7
kde_um <- dd_um$kde
n_um   <- nrow(pts_um)

# Two readable vertical scales. Relative height is the one to use in the
# thesis: it removes the units question and shows only what matters, i.e.
# whether a second local maximum exists and how tall it is.
z_rel_um <- t(kde_um$z) / max(kde_um$z) * 100     # peak = 100
z_cnt_um <- t(kde_um$z) * n_um                    # respondents per unit area

# NOTE the t(): kde2d returns z[i, j] with i indexing x, while plotly expects
# rows indexed by y. Without the transpose the surface is rotated by 90 degrees.

surf_um <- plotly::plot_ly(x = kde_um$x, y = kde_um$y, z = z_rel_um) |>
  plotly::add_surface(
    colorscale = "Viridis",
    colorbar = list(title = "Relative\ndensity\n(peak = 100)")) |>
  plotly::layout(
    #title = "UMAP 2D — kernel density surface",
    scene = list(
      xaxis = list(title = "UMAP 1"),
      yaxis = list(title = "UMAP 2"),
      zaxis = list(title = "Relative density (peak = 100)"),
      aspectmode = "manual",
      aspectratio = list(x = 1, y = 1, z = 0.6),
      camera = list(eye = list(x = 1.6, y = -1.6, z = 0.9))))

surf_um

htmlwidgets::saveWidget(
  plotly::partial_bundle(surf_um),
  file.path(PLOT_DIR, "umap2d_density_surface.html"),
  selfcontained = TRUE)
