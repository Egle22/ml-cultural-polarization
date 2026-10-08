# ═══════════════════════════════════════════════════════════════════════
#  scripts/06_robustness/06_hdbscan_mst_diagnostics.R
#  HDBSCAN with full mutual-reachability / MST diagnostics — WVS7 USA
#
#  What this script produces
#    P1  core distance d_core distribution                (density contrast)
#    P2  raw d vs mutual-reachability d_mreach            (metric inflation)
#    P3  MST edge-weight profile + true bifurcations      (WHERE the cuts are)
#    P4  HDBSCAN condensed tree                           (the dendrogram)
#    P5  single-linkage dendrogram on d_mreach            (shape only)
#    G   3D force-directed graph, neo4j-style, standalone HTML
#
#  Key objects saved for reuse
#    mst_graph         igraph — FULL MST on d_mreach. This is the object
#                      that GENERATES the entire HDBSCAN hierarchy.
#    mst_graph_final   igraph — same MST with between-cluster and
#                      noise-incident edges deleted: its connected
#                      components ARE the extracted clusters.
#    mst_cuts          data.frame — every cut, in order, with the sizes of
#                      the two resulting components and a flag marking
#                      TRUE bifurcations (both sides >= MIN_PTS).
#
#  Thesis  : Appendix A.4.2 (why so many noise points)
#  Input   : results/objects/clusters/gower_dist.rds, results/objects/pca_scores.rds
#            or the standardised value items (INPUT_MODE below)
#  Output  : results/objects/HDBSCAN/{mst_graph,mst_graph_final,mst_cuts,
#            hdbscan_obj,coredist,mrdist}_<tag>.rds; figures and HTML graphs
#            in results/figures/06_robustness/hdbscan/
# ═══════════════════════════════════════════════════════════════════════

source(here::here("R", "00_setup.R"))

suppressPackageStartupMessages({
  library(dbscan)
  library(igraph)
  library(plotly)
  library(htmlwidgets)
  library(ggplot2)
  library(patchwork)
})

# ═══════════════════════════════════════════════════════════════════════
#  0. INPUT CONFIGURATION — THIS IS THE ONLY BLOCK YOU EDIT
# ═══════════════════════════════════════════════════════════════════════

INPUT_MODE  <- "umap"    # "gower" | "pca" | "umap" | "raw" | "custom"
N_DIM       <- 2L         # embedding dimensionality; IGNORED when "gower"
MIN_PTS     <- 25L        # core-distance k (= MinPts-1 NN) AND min cluster size
LAYOUT_MODE <- "force"    # "force" = neo4j-style FR | "embedding" = MDS coords
CUSTOM_PATH <- ""         # only used when INPUT_MODE == "custom"
SEED        <- 42L

#  INPUT_MODE reference
#    "gower"   results/objects/clusters/gower_dist.rds -> dist used AS IS (recommended)
#    "pca"     results/objects/pca_scores.rds          -> first N_DIM columns, Euclidean
#    "umap"    computed here, min_dist = 0      -> N_DIM components, Euclidean
#    "raw"     imputed_labeled.rds              -> one-hot + z-score, Euclidean
#    "custom"  CUSTOM_PATH                      -> dist OR matrix, auto-detected

set.seed(SEED)

RUN_TAG <- paste0(INPUT_MODE,
                  if (INPUT_MODE == "gower") "" else as.character(N_DIM),
                  "_mp", MIN_PTS)

# ── 0b. Paths ─────────────────────────────────────────────────────── #####
# The HTML graphs draw one node per respondent: they are written to
# results/figures (not tracked by git), never to docs/ (GitHub Pages).
OUT_OBJ <- ensure_dir(DIR_OBJ, "HDBSCAN")
OUT_PLT <- ensure_dir(DIR_FIG, "06_robustness", "hdbscan")
OUT_WEB <- ensure_dir(DIR_FIG, "06_robustness", "hdbscan", "html")

# ── 0c. Stale-object guards (you run scripts in chunks) ───────────── #####
rm(list = intersect(
  ls(),
  c("D", "Dm", "Dmr", "core_d", "hdb", "mst_edges", "mst_cuts",
    "mst_graph", "mst_graph_final", "coords3d", "labels_final")
))

# ── 0d. Object resolver (clusters/ takes priority over Objects/) ──── #####
resolve_object <- function(fname) {
  candidates <- c(
    file.path(DIR_CLUST, fname),
    file.path(DIR_OBJ, fname),
    file.path(DIR_DATA, fname)
  )
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0L) {
    stop("Object not found: ", fname, "\n  looked in:\n  ",
         paste(candidates, collapse = "\n  "), call. = FALSE)
  }
  hit[1L]
}

# ═══════════════════════════════════════════════════════════════════════
#  1. INPUT DISPATCH  ->  a single `dist` object D
# ═══════════════════════════════════════════════════════════════════════

message("[1] Building input (mode = ", INPUT_MODE, ")")

D <- switch(
  INPUT_MODE,
  
  "gower" = {
    d <- readRDS(resolve_object("gower_dist.rds"))
    stopifnot(inherits(d, "dist"))
    d
  },
  
  "pca" = {
    sc <- readRDS(resolve_object("pca_scores.rds"))
    sc <- as.matrix(sc)
    stopifnot(ncol(sc) >= N_DIM)
    dist(sc[, seq_len(N_DIM), drop = FALSE])
  },
  
  "umap" = {
    library(umap)
    X <- readRDS(resolve_object("imputed_labeled.rds"))
    X <- model.matrix(~ . - 1, data = as.data.frame(X))
    X <- scale(X)                                   # z-score, not min-max
    X <- X[, apply(X, 2, function(v) all(is.finite(v))), drop = FALSE]
    cfg <- umap::umap.defaults
    cfg$n_neighbors  <- 15
    cfg$min_dist     <-1e-4         # 0, NOT 0.1 — required for density clustering
    cfg$metric       <- "euclidean"
    cfg$n_epochs     <- 500
    cfg$n_components <- N_DIM
    cfg$random_state <- SEED
    dist(umap::umap(X, config = cfg)$layout)
  },
  
  "raw" = {
    X <- readRDS(resolve_object("imputed_labeled.rds"))
    X <- model.matrix(~ . - 1, data = as.data.frame(X))
    X <- scale(X)
    X <- X[, apply(X, 2, function(v) all(is.finite(v))), drop = FALSE]
    dist(X)
  },
  
  "custom" = {
    stopifnot(nzchar(CUSTOM_PATH), file.exists(CUSTOM_PATH))
    obj <- readRDS(CUSTOM_PATH)
    if (inherits(obj, "dist")) obj else dist(as.matrix(obj))
  },
  
  stop("Unknown INPUT_MODE: ", INPUT_MODE)
)

stopifnot(inherits(D, "dist"))
n_obs <- attr(D, "Size")
message("    n = ", n_obs, " observations")

# ═══════════════════════════════════════════════════════════════════════
#  2. CORE DISTANCE + MUTUAL REACHABILITY DISTANCE
#     dbscan::coredist() = distance to the (MinPts - 1)-th neighbor
#     dbscan::mrdist()   = max(coredist(a), coredist(b), dist(a, b))
# ═══════════════════════════════════════════════════════════════════════

message("[2] Core distance and mutual reachability distance")

core_d <- dbscan::coredist(D, minPts = MIN_PTS)
Dmr    <- dbscan::mrdist(D, minPts = MIN_PTS, coredist = core_d)
stopifnot(length(core_d) == n_obs, inherits(Dmr, "dist"))

Dm  <- as.matrix(D)     # n x n raw distances
Dmm <- as.matrix(Dmr)   # n x n mutual reachability distances

# ── Core-dominance fraction: share of pairs where d_mreach is set by the
#    SPARSITY of the endpoints rather than by their actual distance.
#    High value  -> the metric has been flattened -> weak density contrast.
core_dom <- 0; n_pairs <- 0
for (i in seq_len(n_obs - 1L)) {
  j        <- (i + 1L):n_obs
  cm       <- pmax(core_d[i], core_d[j])
  core_dom <- core_dom + sum(cm > Dm[i, j])
  n_pairs  <- n_pairs + length(j)
}
core_dominance <- core_dom / n_pairs

message(sprintf("    core-dominance fraction = %.1f%%", 100 * core_dominance))
message(sprintf("    d_core: mean = %.4f  sd = %.4f  CV = %.3f",
                mean(core_d), sd(core_d), sd(core_d) / mean(core_d)))

# ═══════════════════════════════════════════════════════════════════════
#  3. MINIMUM SPANNING TREE ON d_mreach  (vectorized Prim, O(n^2))
#     This tree generates the ENTIRE HDBSCAN hierarchy.
# ═══════════════════════════════════════════════════════════════════════

message("[3] Minimum spanning tree on d_mreach")

prim_mst <- function(M) {
  n       <- nrow(M)
  in_tree <- logical(n); in_tree[1L] <- TRUE
  best_d  <- M[1L, ]; best_d[1L] <- Inf
  best_p  <- rep(1L, n)
  from <- integer(n - 1L); to <- integer(n - 1L); w <- numeric(n - 1L)
  for (i in seq_len(n - 1L)) {
    j <- which.min(best_d)
    from[i] <- best_p[j]; to[i] <- j; w[i] <- best_d[j]
    in_tree[j] <- TRUE
    best_d[j]  <- Inf
    upd <- (M[j, ] < best_d) & !in_tree
    best_d[upd] <- M[j, upd]
    best_p[upd] <- j
  }
  data.frame(from = from, to = to, weight = w)
}

mst_edges <- prim_mst(Dmm)
stopifnot(nrow(mst_edges) == n_obs - 1L)

# ═══════════════════════════════════════════════════════════════════════
#  4. HDBSCAN EXTRACTION (authoritative labels come from the package)
# ═══════════════════════════════════════════════════════════════════════

message("[4] HDBSCAN extraction (EOM)")

hdb <- dbscan::hdbscan(D, minPts = MIN_PTS)
labels_final <- hdb$cluster                       # 0 = noise
n_clusters   <- length(unique(labels_final[labels_final != 0L]))
noise_frac   <- mean(labels_final == 0L)

message(sprintf("    %d clusters | %.1f%% noise", n_clusters, 100 * noise_frac))

# ── Validation: single-linkage merge heights MUST equal the sorted MST
#    edge weights. If this fails, our MST and dbscan's disagree.
hc_h  <- sort(hdb$hc$height)
mst_h <- sort(mst_edges$weight)
mst_check <- isTRUE(all.equal(mst_h, hc_h, tolerance = 1e-8))
if (!mst_check) {
  warning("MST weights do not match hdbscan$hc$height — inspect before using.")
} else {
  message("    OK: MST weights == hclust merge heights")
}

# ═══════════════════════════════════════════════════════════════════════
#  5. CUT SEQUENCE — which cuts are TRUE bifurcations?
#     Union-find over MST edges in ASCENDING weight order. Each merge,
#     read backwards, is exactly one cut of the tree; the two component
#     sizes at merge time are the sizes the cut would produce.
# ═══════════════════════════════════════════════════════════════════════

message("[5] Cut sequence (union-find)")

cut_sequence <- function(edges, n) {
  edges  <- edges[order(edges$weight), , drop = FALSE]
  parent <- seq_len(n)
  sz     <- rep(1L, n)
  find   <- function(i) { while (parent[i] != i) i <- parent[i]; i }
  sa <- sb <- integer(nrow(edges))
  for (e in seq_len(nrow(edges))) {
    a <- find(edges$from[e]); b <- find(edges$to[e])
    x <- sz[a]; y <- sz[b]
    if (x < y) { tmp <- a; a <- b; b <- tmp; tmp <- x; x <- y; y <- tmp }
    parent[b] <- a
    sz[a]     <- x + y
    sa[e] <- x; sb[e] <- y
  }
  data.frame(weight = edges$weight, size_a = sa, size_b = sb,
             from = edges$from, to = edges$to)
}

mst_cuts <- cut_sequence(mst_edges, n_obs)
mst_cuts <- mst_cuts[order(-mst_cuts$weight), , drop = FALSE]   # cut order
mst_cuts$cut_rank       <- seq_len(nrow(mst_cuts))
mst_cuts$lambda         <- 1 / mst_cuts$weight
mst_cuts$min_side       <- pmin(mst_cuts$size_a, mst_cuts$size_b)
mst_cuts$is_bifurcation <- mst_cuts$min_side >= MIN_PTS
rownames(mst_cuts) <- NULL

n_bifurc <- sum(mst_cuts$is_bifurcation)
message(sprintf("    %d cuts total, of which %d TRUE bifurcations (both sides >= %d)",
                nrow(mst_cuts), n_bifurc, MIN_PTS))
if (n_bifurc > 0L) {
  print(utils::head(
    mst_cuts[mst_cuts$is_bifurcation,
             c("cut_rank", "weight", "lambda", "size_a", "size_b")], 10))
}

# ═══════════════════════════════════════════════════════════════════════
#  6. GRAPH OBJECTS  (full MST + final within-cluster forest)
# ═══════════════════════════════════════════════════════════════════════

message("[6] Building igraph objects")

mst_graph <- igraph::graph_from_data_frame(
  d = mst_edges, directed = FALSE,
  vertices = data.frame(name = as.character(seq_len(n_obs)))
)
igraph::E(mst_graph)$weight  <- mst_edges$weight
igraph::E(mst_graph)$lambda  <- 1 / mst_edges$weight
igraph::V(mst_graph)$cluster <- labels_final
igraph::V(mst_graph)$prob    <- hdb$membership_prob
igraph::V(mst_graph)$glosh   <- hdb$outlier_scores
igraph::V(mst_graph)$coredist <- core_d

# Final graph: delete edges that cross clusters or touch noise.
# Connected components of the result ARE the extracted clusters.
ends_mat  <- igraph::ends(mst_graph, igraph::E(mst_graph), names = FALSE)
cl_from   <- labels_final[ends_mat[, 1]]
cl_to     <- labels_final[ends_mat[, 2]]
drop_edge <- (cl_from != cl_to) | (cl_from == 0L) | (cl_to == 0L)

mst_graph_final <- igraph::delete_edges(mst_graph, igraph::E(mst_graph)[drop_edge])
message(sprintf("    final graph: %d of %d MST edges retained",
                igraph::ecount(mst_graph_final), igraph::ecount(mst_graph)))

# ── 3D coordinates ────────────────────────────────────────────────── #####
message("[6b] 3D layout (", LAYOUT_MODE, ")")
set.seed(SEED)
coords3d <- if (LAYOUT_MODE == "force") {
  igraph::layout_with_fr(mst_graph, dim = 3, niter = 500)
} else {
  cmdscale(D, k = 3)                       # classical MDS, works for Gower too
}
colnames(coords3d) <- c("x", "y", "z")
igraph::V(mst_graph)$x <- igraph::V(mst_graph_final)$x <- coords3d[, 1]
igraph::V(mst_graph)$y <- igraph::V(mst_graph_final)$y <- coords3d[, 2]
igraph::V(mst_graph)$z <- igraph::V(mst_graph_final)$z <- coords3d[, 3]

# ═══════════════════════════════════════════════════════════════════════
#  7. PLOTS
# ═══════════════════════════════════════════════════════════════════════

message("[7] Plots")

noise_col <- "grey75"
pal <- if (exists("cluster_colors")) {
  grDevices::colorRampPalette(unname(cluster_colors))(max(1L, n_clusters))
} else {
  grDevices::colorRampPalette(c("#1D9E75", "#7F77DD", "#D85A30",
                                "#378ADD", "#EF9F27"))(max(1L, n_clusters))
}
names(pal) <- as.character(seq_len(max(1L, n_clusters)))
pal_full <- c("0" = noise_col, pal)

# ── P1. core distance distribution ────────────────────────────────── #####
p1 <- ggplot2::ggplot(data.frame(core_d = core_d), ggplot2::aes(core_d)) +
  ggplot2::geom_histogram(bins = 60, fill = "#378ADD", colour = "white") +
  ggplot2::geom_vline(xintercept = mean(core_d), linetype = 2) +
  ggplot2::labs(
    title    = sprintf("Core distance distribution (minPts = %d)", MIN_PTS),
    subtitle = sprintf("mean = %.4f | sd = %.4f | CV = %.3f  (low CV = no density contrast)",
                       mean(core_d), sd(core_d), sd(core_d) / mean(core_d)),
    x = expression(d[core]), y = "Count") +
  ggplot2::theme_minimal(base_size = 11)

# ── P2. raw d vs d_mreach ─────────────────────────────────────────── #####
set.seed(SEED)
idx  <- sample(length(D), min(3e5, length(D)))
p2df <- rbind(
  data.frame(metric = "d (raw)",     value = as.vector(D)[idx]),
  data.frame(metric = "d_mreach",    value = as.vector(Dmr)[idx])
)
p2 <- ggplot2::ggplot(p2df, ggplot2::aes(value, fill = metric, colour = metric)) +
  ggplot2::geom_density(alpha = 0.35) +
  ggplot2::scale_fill_manual(values = c("d (raw)" = "#888780", "d_mreach" = "#D85A30")) +
  ggplot2::scale_colour_manual(values = c("d (raw)" = "#5F5E5A", "d_mreach" = "#993C1D")) +
  ggplot2::labs(
    title    = "Metric inflation: raw distance vs mutual reachability",
    subtitle = sprintf("core-dominance fraction = %.1f%% of pairs (max set by sparsity, not by distance)",
                       100 * core_dominance),
    x = "Distance", y = "Density", fill = NULL, colour = NULL) +
  ggplot2::theme_minimal(base_size = 11)

# ── P3. MST edge-weight profile + true bifurcations ───────────────── #####
p3 <- ggplot2::ggplot(mst_cuts, ggplot2::aes(cut_rank, weight)) +
  ggplot2::geom_line(colour = "#5F5E5A", linewidth = 0.4) +
  ggplot2::geom_point(
    data = mst_cuts[mst_cuts$is_bifurcation, , drop = FALSE],
    ggplot2::aes(cut_rank, weight), colour = "#D85A30", size = 2) +
  ggplot2::labs(
    title    = "MST edge-weight profile (cut sequence, heaviest first)",
    subtitle = sprintf("%d true bifurcations (orange) out of %d cuts | a KNEE here = real structure, a smooth curve = continuum",
                       n_bifurc, nrow(mst_cuts)),
    x = "Cut rank (order of edge removal)", y = expression(d[mreach]~of~removed~edge)) +
  ggplot2::theme_minimal(base_size = 11)

p_diag <- (p1 | p2) / p3 +
  patchwork::plot_annotation(
    title = sprintf("HDBSCAN density diagnostics — %s", RUN_TAG),
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(size = 13, face = "bold")))

# ── P4/P5. dendrograms (base R — note: do NOT pass `main` to plot.hdbscan) ####
dendro_file <- file.path(OUT_PLT, sprintf("hdbscan_dendrograms_%s.pdf", RUN_TAG))
pdf(dendro_file, width = 12, height = 9)
op <- par(mfrow = c(2, 1), mar = c(3, 4, 3, 1))
plot(hdb, show_flat = TRUE)
title(main = sprintf("Condensed tree — %s (minPts = %d)", RUN_TAG, MIN_PTS))
plot(hdb$hc, labels = FALSE, hang = -1, xlab = "", sub = "",
     main = sprintf("Single-linkage dendrogram on d_mreach — %s", RUN_TAG))
par(op)
dev.off()

# ═══════════════════════════════════════════════════════════════════════
#  8. 3D GRAPH (neo4j style) -> standalone HTML
# ═══════════════════════════════════════════════════════════════════════

message("[8] 3D graph")

build_graph3d <- function(g, title_txt) {
  vdf <- data.frame(
    x = igraph::V(g)$x, y = igraph::V(g)$y, z = igraph::V(g)$z,
    cluster = factor(igraph::V(g)$cluster),
    prob    = round(igraph::V(g)$prob, 3),
    glosh   = round(igraph::V(g)$glosh, 3),
    id      = seq_len(igraph::vcount(g))
  )
  em <- igraph::ends(g, igraph::E(g), names = FALSE)
  
  # one trace for all edges, NA separators between segments
  ex <- as.vector(rbind(vdf$x[em[, 1]], vdf$x[em[, 2]], NA))
  ey <- as.vector(rbind(vdf$y[em[, 1]], vdf$y[em[, 2]], NA))
  ez <- as.vector(rbind(vdf$z[em[, 1]], vdf$z[em[, 2]], NA))
  
  p <- plotly::plot_ly() |>
    plotly::add_trace(
      x = ex, y = ey, z = ez, type = "scatter3d", mode = "lines",
      line = list(color = "rgba(150,150,150,0.25)", width = 1),
      hoverinfo = "none", showlegend = FALSE, name = "MST edges") |>
    plotly::add_trace(
      x = vdf$x, y = vdf$y, z = vdf$z, type = "scatter3d", mode = "markers",
      color = vdf$cluster, colors = pal_full,
      marker = list(size = 3, opacity = 0.9,
                    line = list(width = 0.5, color = "rgba(0,0,0,0.3)")),
      text = paste0("id: ", vdf$id,
                    "<br>cluster: ", vdf$cluster,
                    "<br>P(member): ", vdf$prob,
                    "<br>GLOSH: ", vdf$glosh),
      hoverinfo = "text", name = ~vdf$cluster)
  
  ax <- list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE,
             title = "", showbackground = FALSE)
  plotly::layout(
    p,
    title  = list(text = title_txt, font = list(color = "#EEEEEE", size = 15)),
    scene  = list(xaxis = ax, yaxis = ax, zaxis = ax,
                  bgcolor = "#111318"),
    paper_bgcolor = "#111318",
    legend = list(font = list(color = "#EEEEEE")))
}

g3d_full  <- build_graph3d(
  mst_graph,
  sprintf("MST on d_mreach (generates the hierarchy) — %s", RUN_TAG))

g3d_final <- build_graph3d(
  mst_graph_final,
  sprintf("Final HDBSCAN graph: %d clusters, %.0f%% noise — %s",
          n_clusters, 100 * noise_frac, RUN_TAG))

html_full  <- file.path(OUT_WEB, sprintf("mst_full_%s.html",  RUN_TAG))
html_final <- file.path(OUT_WEB, sprintf("mst_final_%s.html", RUN_TAG))

save_html <- function(widget, path) {
  ok <- tryCatch({
    htmlwidgets::saveWidget(widget, file = normalizePath(path, mustWork = FALSE),
                            selfcontained = TRUE, title = basename(path))
    TRUE
  }, error = function(e) FALSE)
  if (!ok) {
    htmlwidgets::saveWidget(widget, file = normalizePath(path, mustWork = FALSE),
                            selfcontained = FALSE, title = basename(path))
    message("    (pandoc unavailable -> saved with a _files/ dependency folder)")
  }
}

# ═══════════════════════════════════════════════════════════════════════
#  9. CONSOLIDATED SAVES  (all writes happen here, after all computation)
# ═══════════════════════════════════════════════════════════════════════

message("[9] Saving")

ggplot2::ggsave(file.path(OUT_PLT, sprintf("hdbscan_diagnostics_%s.pdf", RUN_TAG)),
                p_diag, width = 13, height = 9)

save_html(g3d_full,  html_full)
save_html(g3d_final, html_final)

saveRDS(mst_graph,       file.path(OUT_OBJ, sprintf("mst_graph_%s.rds",       RUN_TAG)))
saveRDS(mst_graph_final, file.path(OUT_OBJ, sprintf("mst_graph_final_%s.rds", RUN_TAG)))
saveRDS(mst_cuts,        file.path(OUT_OBJ, sprintf("mst_cuts_%s.rds",        RUN_TAG)))
saveRDS(hdb,             file.path(OUT_OBJ, sprintf("hdbscan_obj_%s.rds",     RUN_TAG)))
saveRDS(core_d,          file.path(OUT_OBJ, sprintf("coredist_%s.rds",        RUN_TAG)))
saveRDS(Dmr,             file.path(OUT_OBJ, sprintf("mrdist_%s.rds",          RUN_TAG)))

saveRDS(
  list(
    run_tag          = RUN_TAG,
    input_mode       = INPUT_MODE,
    n_dim            = if (INPUT_MODE == "gower") NA_integer_ else N_DIM,
    min_pts          = MIN_PTS,
    n_obs            = n_obs,
    n_clusters       = n_clusters,
    noise_frac       = noise_frac,
    core_dominance   = core_dominance,
    coredist_cv      = sd(core_d) / mean(core_d),
    n_cuts           = nrow(mst_cuts),
    n_bifurcations   = n_bifurc,
    mst_check_passed = mst_check,
    cluster_scores   = hdb$cluster_scores
  ),
  file.path(OUT_OBJ, sprintf("hdbscan_summary_%s.rds", RUN_TAG))
)

message("Done. Objects -> ", OUT_OBJ)
message("      Plots   -> ", OUT_PLT)
message("      HTML    -> ", OUT_WEB)

# ═══════════════════════════════════════════════════════════════════════
#  RELOAD BLOCK (run this instead of the script to reuse the saved graph)
#
#    RUN_TAG   <- "gower_mp25"
#    OUT_OBJ   <- file.path(DIR_OBJ, "HDBSCAN")
#    mst_graph       <- readRDS(file.path(OUT_OBJ, paste0("mst_graph_",       RUN_TAG, ".rds")))
#    mst_graph_final <- readRDS(file.path(OUT_OBJ, paste0("mst_graph_final_", RUN_TAG, ".rds")))
#    mst_cuts        <- readRDS(file.path(OUT_OBJ, paste0("mst_cuts_",        RUN_TAG, ".rds")))
#
#  Everything you need is on the graph: V()$cluster, V()$prob, V()$glosh,
#  V()$coredist, V()$x/y/z, E()$weight, E()$lambda.
# ═══════════════════════════════════════════════════════════════════════

#Plots 
# ═══════════════════════════════════════════════════════════════════════
#  DISPLAY ALL PLOTS IN R — run after the script and the histogram block
# ═══════════════════════════════════════════════════════════════════════

# ── ggplot / patchwork: explicit print() needed ────────────────────
print(p_hist_compare)   # P0 — confronto Gower / d_core / d_mreach
print(p_diag)           # P1+P2+P3 — density diagnostics + MST cut profile

# ── dendrogrammi base R (condensed tree + single-linkage) ──────────
op <- par(mfrow = c(2, 1), mar = c(3, 4, 3, 1))
plot(hdb, show_flat = TRUE)
title(main = sprintf("Condensed tree — %s (minPts = %d)", RUN_TAG, MIN_PTS))
plot(hdb$hc, labels = FALSE, hang = -1, xlab = "", sub = "",
     main = sprintf("Single-linkage dendrogram — %s", RUN_TAG))
par(op)

# ── 3D graph — opens in the RStudio Viewer pane ─────────────────────
g3d_full     # full MST (also shows the edges towards noise)
g3d_final    # final graph with the extracted clusters only (for slides)

htmlwidgets::saveWidget(g3d_final, tempfile(fileext = ".html"), selfcontained = TRUE) |>
  browseURL()
