# ═══════════════════════════════════════════════════════════════════════
#  scripts/06_robustness/09_hdbscan_pruning_animation.R
#
#  Same working plot as file 16 (single SVG trace + ids + manual legend, no
#  scattergl), but with a clearer default colouring.
#
#  COLOR_MODE:
#    "fate"  (default, legible) -- colour a point by its FINAL cluster label;
#            grey it out only once its current component drops below minPts.
#            Frame 0 shows the 6 clusters at once; then points switch off as
#            they fall out (centre first, cores last). The frame at which a
#            point greys = its lambda_death. This is the static peeling map,
#            animated.
#    "bloom" (original)         -- colour by the PLURALITY final label of the
#            current component. With 69% noise the whole cloud starts and ends
#            grey and clusters only flash mid-timeline. Correct, but hard to
#            read. Kept for reference.
# ═══════════════════════════════════════════════════════════════════════

source(here::here("R", "00_setup.R"))

suppressPackageStartupMessages({ library(dbscan); library(plotly); library(htmlwidgets) })

# ── CONFIG ─────────────────────────────────────────────────────────────
ANALYSIS_EMB_PATH <- file.path(DIR_OBJ, "UMAP", "umap_5d.rds")
CANVAS_EMB_PATH   <- file.path(DIR_OBJ, "UMAP", "umap_2d.rds")
CANVAS_MODE     <- "separate"     # "separate" = umap_2d ; "project2" = first 2 of 5D
COLOR_MODE      <- "fate"         # "fate" (clear) or "bloom" (original)
MIN_PTS         <- 25
N_FRAMES        <- 60
POINT_SUBSAMPLE <- NULL
DOCS_DIR <- ensure_dir(DIR_FIG, "06_robustness", "hdbscan", "html")
dir.create(DOCS_DIR, recursive = TRUE, showWarnings = FALSE)
noise_col <- "grey80"

# ── LOAD + MST (5D) ────────────────────────────────────────────────────
Xa <- as.matrix(readRDS(ANALYSIS_EMB_PATH)); storage.mode(Xa) <- "double"; n <- nrow(Xa)
XY <- if (CANVAS_MODE == "separate") as.matrix(readRDS(CANVAS_EMB_PATH))[, 1:2, drop = FALSE] else
  Xa[, 1:2, drop = FALSE]
colnames(XY) <- c("UMAP1", "UMAP2")

knn    <- dbscan::kNN(Xa, k = MIN_PTS - 1L); d_core <- knn$dist[, MIN_PTS - 1L]
D  <- as.matrix(dist(Xa)); Dm <- pmax(D, outer(d_core, d_core, pmax)); diag(Dm) <- 0
in_tree <- logical(n); in_tree[1L] <- TRUE
best_w <- Dm[1L, ]; best_w[1L] <- Inf; best_from <- rep(1L, n)
mst <- matrix(NA_real_, n - 1L, 3L, dimnames = list(NULL, c("from", "to", "weight")))
for (i in seq_len(n - 1L)) {
  j <- which.min(ifelse(in_tree, Inf, best_w))
  mst[i, ] <- c(best_from[j], j, best_w[j]); in_tree[j] <- TRUE
  upd <- (!in_tree) & (Dm[j, ] < best_w); best_w[upd] <- Dm[j, upd]; best_from[upd] <- j
}
i_mst <- as.integer(mst[, "from"]); j_mst <- as.integer(mst[, "to"]); w_mst <- mst[, "weight"]

ref <- dbscan::hdbscan(Xa, minPts = MIN_PTS)$cluster           # dbscan's own labels
K   <- max(ref)

components_at <- function(eps) {
  parent <- seq_len(n); sz <- rep(1L, n)
  find <- function(i) { while (parent[i] != i) i <- parent[i]; i }
  for (e in which(w_mst <= eps)) {
    ra <- find(i_mst[e]); rb <- find(j_mst[e])
    if (ra != rb) { if (sz[ra] < sz[rb]) { t<-ra;ra<-rb;rb<-t }; parent[rb] <- ra; sz[ra] <- sz[ra] + sz[rb] }
  }
  vapply(seq_len(n), find, integer(1))
}

base_cols <- grDevices::hcl.colors(max(3L, K), "Dark 3")[seq_len(K)]
pal  <- c(setNames(base_cols, paste0("C", seq_len(K))), noise = noise_col)
lvls <- c(paste0("C", seq_len(K)), "noise")
use_idx  <- if (is.null(POINT_SUBSAMPLE)) seq_len(n) else sort(sample(seq_len(n), min(POINT_SUBSAMPLE, n)))
eps_grid <- rev(exp(seq(log(as.numeric(quantile(w_mst, 0.001))), log(max(w_mst)), length.out = N_FRAMES)))

# ── BUILD FRAMES ───────────────────────────────────────────────────────
frames <- vector("list", N_FRAMES)
for (f in seq_len(N_FRAMES)) {
  comp <- components_at(eps_grid[f]); csz <- ave(comp, comp, FUN = length)
  
  if (COLOR_MODE == "fate") {
    # colour by FINAL label; grey once the point's component is < minPts
    collab <- ifelse(ref != 0L & csz >= MIN_PTS, paste0("C", ref), "noise")
  } else {                                       # "bloom" : plurality label
    plur    <- tapply(ref, comp, function(v) which.max(tabulate(v + 1L)) - 1L)
    plur_pt <- plur[as.character(comp)]
    collab  <- ifelse(csz >= MIN_PTS & plur_pt != 0L, paste0("C", plur_pt), "noise")
  }
  
  frames[[f]] <- data.frame(
    id = use_idx, UMAP1 = XY[use_idx, 1], UMAP2 = XY[use_idx, 2],
    frame = factor(sprintf("%03d", f), levels = sprintf("%03d", seq_len(N_FRAMES))),
    col   = unname(pal[collab[use_idx]]),
    hover = sprintf("comp %d | final %s | lambda %.3f",
                    csz[use_idx], ifelse(ref[use_idx] == 0, "noise", paste0("C", ref[use_idx])),
                    1 / eps_grid[f]))
}
anim_df <- do.call(rbind, frames)

# ── PLOT: one animated trace + manual legend ───────────────────────────
p <- plot_ly()
p <- add_trace(p, data = anim_df, x = ~UMAP1, y = ~UMAP2, frame = ~frame, ids = ~id,
               type = "scatter", mode = "markers",
               marker = list(size = 5, opacity = 0.75, color = ~col, line = list(width = 0)),
               text = ~hover, hoverinfo = "text", showlegend = FALSE)
for (lab in lvls)
  p <- add_trace(p, x = c(NA), y = c(NA), type = "scatter", mode = "markers", name = lab,
                 marker = list(size = 9, color = unname(pal[lab])),
                 hoverinfo = "skip", showlegend = TRUE, inherit = FALSE)
p <- layout(p,
            title = list(text = sprintf("HDBSCAN pruning on 5D UMAP (%s mode, minPts = %d)", COLOR_MODE, MIN_PTS)),
            xaxis = list(title = "UMAP1", zeroline = FALSE),
            yaxis = list(title = "UMAP2", zeroline = FALSE, scaleanchor = "x"),
            legend = list(title = list(text = "final cluster")))
p <- animation_opts(p, frame = 130, transition = 0, redraw = TRUE)
p <- animation_slider(p, currentvalue = list(prefix = "frame (low lambda -> high lambda): "))

htmlwidgets::saveWidget(p, file.path(DOCS_DIR,
                                     sprintf("pruning_animation_%s_umap5_minpts25.html", COLOR_MODE)),
                        selfcontained = TRUE, title = "HDBSCAN pruning animation")
message(sprintf("Done (%s mode). HTML in results/figures/06_robustness/hdbscan/html/", COLOR_MODE))
