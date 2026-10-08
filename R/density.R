#### 99_density_helpers.R — reusable 2D density machinery ####

# Deterministic. No sampling, no jitter, no subsetting. See notes below.

density_hdr <- function(x, y,
                        n_grid    = 400,
                        pad_frac  = 0.10,
                        bw_factor = 1.0,
                        probs     = c(0.25, 0.50, 0.75, 0.90)) {
  
  stopifnot(is.numeric(x), is.numeric(y), length(x) == length(y),
            !anyNA(x), !anyNA(y), all(probs > 0), all(probs < 1))
  
  h <- c(MASS::bandwidth.nrd(x), MASS::bandwidth.nrd(y)) * bw_factor
  stopifnot(all(h > 0))
  
  pad <- function(v) range(v) + c(-1, 1) * diff(range(v)) * pad_frac
  kde <- MASS::kde2d(x, y, h = h, n = n_grid, lims = c(pad(x), pad(y)))
  
  grid <- expand.grid(X = kde$x, Y = kde$y)
  grid$z <- as.vector(kde$z)
  
  # Density at each observation (nearest grid node) -> HDR levels.
  ix <- findInterval(x, kde$x, all.inside = TRUE)
  iy <- findInterval(y, kde$y, all.inside = TRUE)
  d_obs <- kde$z[cbind(ix, iy)]
  
  lev <- stats::quantile(d_obs, probs = 1 - probs, names = FALSE)
  names(lev) <- paste0(probs * 100, "%")
  stopifnot(!any(duplicated(round(lev, 12))))
  
  list(kde = kde, grid = grid, levels = lev, bandwidth = h,
       dens_at_point = d_obs)
}

plot_hdr <- function(dd, pts, xlab = "Dim 1", ylab = "Dim 2",
                     style = c("hdr", "paper")) {
  
  style <- match.arg(style)
  lev_sorted <- sort(dd$levels)
  lab_sorted <- rev(names(dd$levels))
  cols       <- rev(c("#08306B", "#2171B5", "#6BAED6", "#C6DBEF"))
  
  if (style == "hdr") {
    ggplot() +
      geom_point(data = pts, aes(X, Y),
                 colour = "grey55", size = 0.55, alpha = 0.30) +
      geom_contour(data = dd$grid,
                   aes(X, Y, z = z, colour = factor(after_stat(level))),
                   breaks = lev_sorted, linewidth = 0.75) +
      scale_colour_manual(name = "Share of\nrespondents",
                          values = stats::setNames(cols, lev_sorted),
                          labels = lab_sorted) +
      coord_equal() + labs(x = xlab, y = ylab) +
      theme_minimal(base_size = 12) +
      theme(panel.grid.minor = element_blank())
  } else {
    # Deliberately mimics the paper's rendering: many equally spaced filled
    # bands + red points. Used ONLY to separate rendering from normalisation.
    ggplot() +
      geom_contour_filled(data = dd$grid, aes(X, Y, z = z), bins = 14) +
      geom_point(data = pts, aes(X, Y), colour = "red", size = 0.7) +
      scale_fill_brewer(palette = "Blues", guide = "none") +
      coord_equal(expand = FALSE) + labs(x = xlab, y = ylab) +
      theme_minimal(base_size = 12) +
      theme(panel.grid = element_blank())
  }
}