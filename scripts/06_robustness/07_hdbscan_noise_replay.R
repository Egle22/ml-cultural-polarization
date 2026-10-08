# ═══════════════════════════════════════════════════════════════════════
#  scripts/06_robustness/07_hdbscan_noise_replay.R
#
#  ONE unified script covering:
#    (A) the MST -> condensed-tree machinery, in a SINGLE replay pass
#    (B) lambda_birth / lambda_death per cluster, stability, EOM selection,
#        and the reconstructed labels cross-checked against dbscan::hdbscan
#    (C) the diagnostic plots (lambda_death distribution, cluster lifespans,
#        deaths-per-iteration, stability bars, peeling map, final labels)
#    (D) an INTERACTIVE plotly animation of the pruning: as edges are cut
#        (lambda rises), clusters "bloom" and noise spreads, drawn on the
#        2D canvas for the 5D-UMAP clustering.
#
#  ── SECTION MAP (delete any OPTIONAL block you do not need) ───────────
#    0  config
#    1  load embeddings (analysis = 5D, canvas = 2D)
#    2  core distance, mutual reachability, MST         (on the 5D space)
#    3  SINGLE replay pass -> events + condensed tree
#    4  stability S(C), EOM selection, reconstructed labels
#    5  cross-check vs dbscan::hdbscan  (source of truth for point labels)
#    6  prints: verdict, selection audit, radial test
#    7  static plots                                    (OPTIONAL, cheap)
#    8  INTERACTIVE animation  (the headline)
#    9  consolidated saves
#
#  Thesis  : Appendix A.4.2 (why so many noise points)
#  Input   : results/objects/UMAP/{umap_5d,umap_2d}.rds (from 02_clustering/05_hdbscan.R)
#  Output  : results/objects/HDBSCAN/*.rds; figures and HTML animation in
#            results/figures/06_robustness/hdbscan/ (one point per respondent:
#            never published in docs/)
# ═══════════════════════════════════════════════════════════════════════

source(here::here("R", "00_setup.R"))

suppressPackageStartupMessages({
  library(dbscan)
  library(ggplot2)
  library(patchwork)
  library(plotly)
  library(htmlwidgets)
})

# ── 0. CONFIG ──────────────────────────────────────────────────────────
# Two embeddings are needed: the 5D one where clustering happens, and the
# 2D one used only as a drawing canvas. Adjust the paths to your filenames.
ANALYSIS_EMB_PATH <- file.path(DIR_OBJ, "UMAP", "umap_5d.rds")  # 5D
CANVAS_EMB_PATH   <- file.path(DIR_OBJ, "UMAP", "umap_2d.rds")  # 2D

# CANVAS_MODE:
#   "separate"  -> draw on the independent 2D UMAP (your convention for all
#                  other plots). Faithful to the pipeline, but edges/clusters
#                  may look scattered: 5D neighbours are not 2D neighbours.
#   "project2"  -> draw on the first 2 of the 5 analysis dimensions. Only 2 of
#                  5 dims, but the same space the MST lives in, so splits look
#                  spatially cleaner. Does NOT need CANVAS_EMB_PATH.
CANVAS_MODE <- "separate"

MIN_PTS  <- 25          # core-distance k / min cluster size (as in script 11)
N_FRAMES <- 60          # animation frames (log-spaced in lambda)
POINT_SUBSAMPLE <- NULL # set e.g. 1500 if the HTML is too heavy; NULL = all
N_DIM    <- 5
RUN_TAG  <- paste0("umap", N_DIM, "_minpts", MIN_PTS)

OBJ_DIR  <- ensure_dir(DIR_OBJ, "HDBSCAN")
PLOT_DIR <- ensure_dir(DIR_FIG, "06_robustness", "hdbscan")
DOCS_DIR <- ensure_dir(DIR_FIG, "06_robustness", "hdbscan", "html")
for (d in c(OBJ_DIR, PLOT_DIR, DOCS_DIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

noise_col <- "grey80"

# ── STALE-OBJECT GUARD ─────────────────────────────────────────────────
rm(list = intersect(ls(), c("Xa", "XY", "D", "Dm", "d_core", "mst_edges",
                            "events", "nodes_df", "diag_pt", "anim_df", "ref")))

# ── 1. LOAD EMBEDDINGS ─────────────────────────────────────────────────
stopifnot(file.exists(ANALYSIS_EMB_PATH))
Xa <- readRDS(ANALYSIS_EMB_PATH)
if (is.data.frame(Xa)) Xa <- as.matrix(Xa)
storage.mode(Xa) <- "double"
stopifnot(is.matrix(Xa), ncol(Xa) == N_DIM, nrow(Xa) > MIN_PTS)
n <- nrow(Xa)

if (CANVAS_MODE == "separate") {
  stopifnot(file.exists(CANVAS_EMB_PATH))
  XY <- readRDS(CANVAS_EMB_PATH)
  if (is.data.frame(XY)) XY <- as.matrix(XY)
  stopifnot(nrow(XY) == n, ncol(XY) >= 2)
  XY <- XY[, 1:2, drop = FALSE]
} else {
  XY <- Xa[, 1:2, drop = FALSE]
}
colnames(XY) <- c("UMAP1", "UMAP2")
message(sprintf("Analysis: n = %d, d = %d | canvas mode = %s", n, ncol(Xa), CANVAS_MODE))

# ── 2. CORE DISTANCE, MUTUAL REACHABILITY, MST (all on the 5D space) ────
# core distance = distance to the (minPts-1)-th neighbour EXCLUDING self.
knn    <- dbscan::kNN(Xa, k = MIN_PTS - 1L)
d_core <- knn$dist[, MIN_PTS - 1L]

D  <- as.matrix(dist(Xa))
Dm <- pmax(D, outer(d_core, d_core, pmax))
diag(Dm) <- 0

# Prim MST on mutual reachability (O(n^2), memory-light)
in_tree   <- logical(n); in_tree[1L] <- TRUE
best_w    <- Dm[1L, ];   best_w[1L]  <- Inf
best_from <- rep(1L, n)
mst_edges <- matrix(NA_real_, nrow = n - 1L, ncol = 3L,
                    dimnames = list(NULL, c("from", "to", "weight")))
for (i in seq_len(n - 1L)) {
  j <- which.min(ifelse(in_tree, Inf, best_w))
  mst_edges[i, ] <- c(best_from[j], j, best_w[j])
  in_tree[j] <- TRUE
  upd <- (!in_tree) & (Dm[j, ] < best_w)
  best_w[upd]    <- Dm[j, upd]
  best_from[upd] <- j
}
stopifnot(all(is.finite(mst_edges[, "weight"])))

i_mst <- as.integer(mst_edges[, "from"])
j_mst <- as.integer(mst_edges[, "to"])
w_mst <- mst_edges[, "weight"]

# core-dominance on the ONLY pairs HDBSCAN uses (the MST edges)
core_dom_mst <- mean(w_mst > D[cbind(i_mst, j_mst)] + 1e-12)
message(sprintf("core-dominance on MST edges = %.3f", core_dom_mst))

# ── 3. SINGLE REPLAY PASS: events + condensed tree ─────────────────────
# Kruskal in INCREASING weight = pruning in DECREASING weight = TOP-DOWN in
# lambda. A merge of components (a, b) at weight w IS the split, at
# lambda = 1/w, of a parent into two children of size a, b.
#   both a,b >= minPts   -> TRUE BIFURCATION (two clusters born)
#   one side  <  minPts  -> FALL-OUT   (the small side leaves its cluster)
#   both <minPts, sum>=mP -> DISSOLUTION (a leaf cluster ends; all points exit)
ord       <- order(w_mst, decreasing = FALSE)
parent_uf <- seq_len(n); csize <- rep(1L, n)
members   <- split(seq_len(n), seq_len(n))
node_of   <- rep(NA_integer_, n)
find2 <- function(i) { while (parent_uf[i] != i) i <- parent_uf[i]; i }

MAXN       <- n
nd_birth   <- rep(NA_real_,    MAXN); nd_death  <- rep(NA_real_,    MAXN)
nd_size    <- rep(NA_integer_, MAXN); nd_parent <- rep(NA_integer_, MAXN)
nd_child1  <- rep(NA_integer_, MAXN); nd_child2 <- rep(NA_integer_, MAXN)
nd_sumfall <- rep(0,           MAXN); nd_nfall  <- rep(0L,          MAXN)
n_nodes    <- 0L

exit_node <- rep(NA_integer_, n)
lam_exit  <- rep(NA_real_,    n)          # per-point lambda_death

events <- data.frame(step = seq_len(n - 1L), weight = NA_real_, lambda = NA_real_,
                     size_small = NA_integer_, size_large = NA_integer_,
                     n_big = NA_integer_, max_comp = NA_integer_,
                     n_death = 0L, n_born = 0L)
n_big <- 0L; max_comp_running <- 1L

for (s in seq_len(n - 1L)) {
  w   <- w_mst[ord[s]]
  lam <- if (w > 0) 1 / w else Inf
  ra  <- find2(i_mst[ord[s]]); rb <- find2(j_mst[ord[s]])
  a   <- csize[ra]; b <- csize[rb]; newsz <- a + b
  big_a <- a >= MIN_PTS; big_b <- b >= MIN_PTS
  new_node <- NA_integer_
  
  if (newsz >= MIN_PTS) {
    if (big_a && big_b) {                       # TRUE BIFURCATION
      n_nodes <- n_nodes + 1L; new_node <- n_nodes
      ca <- node_of[ra]; cb <- node_of[rb]
      nd_death[new_node] <- lam
      nd_child1[new_node] <- ca; nd_child2[new_node] <- cb
      nd_parent[ca] <- new_node; nd_parent[cb] <- new_node
      nd_birth[ca]  <- lam; nd_birth[cb] <- lam
      nd_size[ca]   <- a;   nd_size[cb]  <- b
      events$n_born[s] <- 2L
    } else if (big_a || big_b) {                # FALL-OUT
      if (big_a) { keep <- ra; drp <- rb; kd <- b } else { keep <- rb; drp <- ra; kd <- a }
      nd <- node_of[keep]; pts <- members[[drp]]
      exit_node[pts] <- nd; lam_exit[pts] <- lam
      nd_sumfall[nd] <- nd_sumfall[nd] + kd * lam
      nd_nfall[nd]   <- nd_nfall[nd] + kd
      new_node <- nd; events$n_death[s] <- kd
    } else {                                    # DISSOLUTION (leaf ends)
      n_nodes <- n_nodes + 1L; new_node <- n_nodes
      nd_death[new_node] <- lam
      pts <- c(members[[ra]], members[[rb]])
      exit_node[pts] <- new_node; lam_exit[pts] <- lam
      nd_sumfall[new_node] <- newsz * lam; nd_nfall[new_node] <- newsz
      events$n_death[s] <- newsz
    }
  }
  
  # bookkeeping for the "broom" test
  was_big <- big_a + big_b
  n_big   <- n_big - was_big + as.integer(newsz >= MIN_PTS)
  max_comp_running <- max(max_comp_running, newsz)
  
  # union by size, carrying the node id onto the surviving root
  if (a < b) { t<-ra; ra<-rb; rb<-t; t<-a; a<-b; b<-t }
  parent_uf[rb] <- ra; csize[ra] <- newsz
  members[[ra]] <- c(members[[ra]], members[[rb]]); members[[rb]] <- integer(0)
  node_of[ra]   <- new_node
  
  events[s, c("weight","lambda","size_small","size_large","n_big","max_comp")] <-
    list(w, lam, min(a, b), max(a, b), n_big, max_comp_running)
}
events$true_bifurcation <- events$n_born == 2L
events$rank             <- nrow(events) - events$step + 1L   # rank 1 = first cut

root_id <- n_nodes
stopifnot(node_of[find2(1L)] == root_id)
nd_birth[root_id] <- 0; nd_size[root_id] <- n
lam_fallout <- lam_exit                                       # alias

# invariants
stopifnot(sum(events$n_death) == n)
stopifnot(sum(events$n_born) + 1L == n_nodes)
stopifnot(!anyNA(exit_node), !anyNA(lam_exit))

# ── 4. STABILITY + EOM SELECTION + RECONSTRUCTED LABELS ────────────────
nd_ids      <- seq_len(n_nodes)
children_of <- lapply(nd_ids, function(i) { ch <- c(nd_child1[i], nd_child2[i]); ch[!is.na(ch)] })

# S(C) = sum_{x in C} (lambda_death(x) - lambda_birth(C))
S <- rep(NA_real_, n_nodes)
for (nd in nd_ids) {
  b0 <- nd_birth[nd]
  val <- nd_sumfall[nd] - nd_nfall[nd] * b0
  for (ch in children_of[[nd]]) val <- val + (nd_death[nd] - b0) * nd_size[ch]
  S[nd] <- val
}

descendants <- function(i) {
  out <- integer(0); stk <- children_of[[i]]
  while (length(stk)) { v <- stk[1L]; stk <- stk[-1L]; out <- c(out, v); stk <- c(stk, children_of[[v]]) }
  out
}
selected <- rep(TRUE, n_nodes); selected[root_id] <- FALSE
S_hat <- S
for (nd in setdiff(nd_ids, root_id)) {          # children before parents (id order)
  ch <- children_of[[nd]]; if (!length(ch)) next
  if (sum(S_hat[ch]) > S_hat[nd]) { selected[nd] <- FALSE; S_hat[nd] <- sum(S_hat[ch]) }
  else                            { selected[descendants(nd)] <- FALSE }
}
sel_ids <- which(selected)

nearest_sel <- function(nd) { cur <- nd; while (!is.na(cur)) { if (selected[cur]) return(cur); cur <- nd_parent[cur] }; NA_integer_ }
sel_of_node <- vapply(nd_ids, nearest_sel, integer(1))
lab_rec     <- ifelse(is.na(sel_of_node[exit_node]), 0L,
                      match(sel_of_node[exit_node], sel_ids))

# ── 5. dbscan::hdbscan  (SOURCE OF TRUTH for point labels) ─────────────
hdb <- dbscan::hdbscan(Xa, minPts = MIN_PTS)
ref <- hdb$cluster                              # 0 = noise, 1..K = clusters
K   <- max(ref)

cat("\n══════════════ CONDENSED TREE vs dbscan ══════════════\n")
cat(sprintf("nodes: %d (%d leaves, %d internal) | true bifurcations: %d\n",
            n_nodes, sum(lengths(children_of) == 0L),
            sum(lengths(children_of) > 0L), sum(events$true_bifurcation)))
cat(sprintf("reconstructed: %d clusters, %.1f%% noise\n",
            length(sel_ids), 100 * mean(lab_rec == 0L)))
cat(sprintf("dbscan       : %d clusters, %.1f%% noise\n", K, 100 * mean(ref == 0L)))
cat("cross-tab (rows = reconstructed, cols = dbscan):\n")
print(table(reconstructed = lab_rec, dbscan = ref))
cat("══════════════════════════════════════════════════════\n\n")

# ── 6. PRINTS: verdict, selection audit, radial test ───────────────────
if (max(events$n_big) <= 1L) {
  cat("VERDICT: broom. The root never bifurcates -> noise by construction.\n\n")
} else {
  cat(sprintf("VERDICT: not a broom. %d true bifurcations, up to %d eligible components\n",
              sum(events$true_bifurcation), max(events$n_big)))
  cat("coexisting. Noise comes from fall-out; EOM keeps only the stable subset.\n\n")
}

# selection audit from dbscan (size, stability, mean membership prob)
audit <- data.frame(
  cluster    = as.integer(names(hdb$cluster_scores)),
  stability  = as.numeric(hdb$cluster_scores),
  size       = as.integer(table(factor(ref[ref != 0], levels = seq_len(K)))),
  mean_prob  = tapply(hdb$membership_prob, ref, mean)[as.character(seq_len(K))]
)
audit <- audit[order(-audit$stability), ]
cat("Selection audit (dbscan selected clusters):\n")
print(audit, row.names = FALSE, digits = 4)

# reconstructed condensed-tree nodes (full picture, selected + not)
lam_min <- min(events$lambda[is.finite(events$lambda)])
nodes_df <- data.frame(node = nd_ids, lam_birth = nd_birth[nd_ids],
                       lam_death = nd_death[nd_ids], size_birth = nd_size[nd_ids],
                       n_fallout = nd_nfall[nd_ids], stability = S,
                       selected = selected, is_root = nd_ids == root_id)
nodes_df$birth_plot <- pmax(nodes_df$lam_birth, lam_min)
nodes_df$ypos       <- rank(nodes_df$birth_plot, ties.method = "first")

# radial "periphery-first" test, done in the 5D space, per cluster
rad_local <- rep(NA_real_, n)
for (k in seq_len(K)) {
  idx <- which(ref == k); cen <- colMeans(Xa[idx, , drop = FALSE])
  rad_local[idx] <- sqrt(rowSums(sweep(Xa[idx, , drop = FALSE], 2, cen, "-")^2))
}
rho_local <- sapply(seq_len(K), function(k)
  suppressWarnings(cor(rad_local[ref == k], lam_fallout[ref == k], method = "spearman")))
names(rho_local) <- paste0("cluster_", seq_len(K))
cat("\nSpearman(radius from OWN 5D centroid, lambda_death), per cluster:\n")
print(round(rho_local, 3))
cat("Strongly negative => each cluster is peeled from its periphery inward.\n\n")

# ── 7. STATIC PLOTS (OPTIONAL) ─────────────────────────────────────────
thm <- theme_minimal(base_size = 10) + theme(panel.grid.minor = element_blank(),
                                             legend.position = "bottom")
diag_pt <- data.frame(UMAP1 = XY[, 1], UMAP2 = XY[, 2], lam_death = lam_fallout,
                      ref = factor(ref), fate = factor(ifelse(ref == 0L, "noise", "clustered"),
                                                       levels = c("noise", "clustered")))
fate_cols <- c(noise = "grey65", clustered = "#B22222")
bif_lam   <- events$lambda[events$true_bifurcation]

# Q1 — distribution of lambda_death (one value per point), split by fate
q1 <- ggplot(diag_pt[is.finite(diag_pt$lam_death), ], aes(lam_death, fill = fate)) +
  geom_histogram(bins = 60) +
  geom_vline(xintercept = bif_lam, linetype = "dotted", colour = "steelblue", linewidth = 0.3) +
  scale_x_log10() + scale_fill_manual(values = fate_cols, name = NULL) +
  labs(title = expression("Q1. distribution of " * lambda[death]),
       subtitle = "dotted = the true bifurcations. Noise dies early (low lambda); clustered points die late.",
       x = expression(lambda[death] ~ "(log)"), y = "points") + thm

# Q2 — cluster lifespans: lambda_birth -> lambda_death, EOM selection in red
q2 <- ggplot(nodes_df) +
  geom_segment(aes(x = birth_plot, xend = lam_death, y = ypos, yend = ypos,
                   colour = selected), linewidth = 1.2) +
  geom_point(aes(x = birth_plot, y = ypos), size = 0.8, colour = "black") +
  scale_x_log10(expand = expansion(mult = c(0.02, 0.12))) +
  scale_colour_manual(values = c(`FALSE` = "grey70", `TRUE` = "#B22222"),
                      labels = c("not selected", "selected by EOM"), name = NULL) +
  labs(title = "Q2. cluster lifespans (condensed tree)",
       subtitle = "each segment = one node. Colored length ~ how the animation blooms/greys.",
       x = expression(lambda ~ "(log)"), y = "node (ordered by birth)") + thm

# Q3 — deaths per pruning iteration: where do most exits happen
q3 <- ggplot(events, aes(rank, n_death)) + geom_col(width = 1, fill = "grey35") +
  scale_y_sqrt() +
  labs(title = "Q3. deaths per pruning iteration",
       subtitle = "sqrt scale. Floor at 1 = point-by-point peeling; spikes = leaf clusters dissolving.",
       x = "pruning rank (1 = heaviest edge cut first)", y = "deaths (sqrt)") + thm

# Q4 — stability bars, selected vs not (root excluded: it is never eligible)
q4 <- ggplot(subset(nodes_df, !is_root),
             aes(reorder(factor(node), -stability), stability, fill = selected)) +
  geom_col() +
  scale_fill_manual(values = c(`FALSE` = "grey75", `TRUE` = "#B22222"),
                    labels = c("not selected", "selected"), name = NULL) +
  labs(title = "Q4. node stability S(C)",
       subtitle = "EOM keeps a node only if it beats the sum of its selected children.",
       x = "condensed-tree node", y = "stability") +
  thm + theme(axis.text.x = element_blank())

panel_lambda <- (q1 | q2) / (q3 | q4)
ggsave(file.path(PLOT_DIR, paste0("lambda_selection_", RUN_TAG, ".pdf")),
       panel_lambda, width = 13, height = 9)

# peeling map + final labels on the canvas
m1 <- ggplot(diag_pt, aes(UMAP1, UMAP2, colour = log10(lam_death))) +
  geom_point(size = 0.6) + scale_colour_viridis_c(option = "magma") + coord_equal() +
  labs(title = "peeling: lambda_death per point", colour = expression(log[10](lambda[death]))) + thm
m2 <- ggplot(diag_pt, aes(UMAP1, UMAP2, colour = ref)) +
  geom_point(size = 0.6) + coord_equal() +
  scale_colour_manual(values = c("0" = noise_col,
                                 setNames(scales::hue_pal()(K), as.character(seq_len(K)))),
                      name = "cluster") +
  labs(title = "final HDBSCAN labels (0 = noise)") + thm
ggsave(file.path(PLOT_DIR, paste0("maps_", RUN_TAG, ".pdf")), (m1 | m2), width = 12, height = 5.5)

# ── 8. INTERACTIVE ANIMATION (the headline) ────────────────────────────
# frame f: keep MST edges with weight <= eps_f (i.e. cut everything heavier).
# As eps falls (lambda rises), the tree fragments. A point is coloured by the
# FINAL cluster that is the PLURALITY of its current component, and grey while
# its component is < minPts or plurality-noise. Clusters therefore BLOOM at
# their birth lambda and GREY at their death -- a spatial rendering of Q2.
eps_grid <- rev(exp(seq(log(as.numeric(quantile(w_mst, 0.001))),
                        log(max(w_mst)), length.out = N_FRAMES)))

components_at <- function(eps) {
  parent <- seq_len(n); sz <- rep(1L, n)
  find <- function(i) { while (parent[i] != i) i <- parent[i]; i }
  for (e in which(w_mst <= eps)) {
    ra <- find(i_mst[e]); rb <- find(j_mst[e])
    if (ra != rb) { if (sz[ra] < sz[rb]) { t<-ra;ra<-rb;rb<-t }
      parent[rb] <- ra; sz[ra] <- sz[ra] + sz[rb] }
  }
  vapply(seq_len(n), find, integer(1))
}

# palette: one colour per final cluster + grey for noise/unresolved
base_cols <- if (K <= 8) RColorBrewer::brewer.pal(max(3L, K), "Dark2")[seq_len(K)] else
  grDevices::hcl.colors(K, "Dark 3")
pal <- c(setNames(base_cols, paste0("C", seq_len(K))), noise = noise_col)
lvls <- c(paste0("C", seq_len(K)), "noise")

use_idx <- if (is.null(POINT_SUBSAMPLE)) seq_len(n) else
  sort(sample(seq_len(n), min(POINT_SUBSAMPLE, n)))

frames <- vector("list", N_FRAMES)
for (f in seq_len(N_FRAMES)) {
  comp <- components_at(eps_grid[f])
  csz  <- ave(comp, comp, FUN = length)
  plur <- tapply(ref, comp, function(v) { tb <- tabulate(v + 1L); (which.max(tb) - 1L) })
  plur_pt <- plur[as.character(comp)]
  collab <- ifelse(csz >= MIN_PTS & plur_pt != 0L, paste0("C", plur_pt), "noise")
  frames[[f]] <- data.frame(
    UMAP1 = XY[use_idx, 1], UMAP2 = XY[use_idx, 2],
    frame = factor(sprintf("%03d", f), levels = sprintf("%03d", seq_len(N_FRAMES))),
    color_label = factor(collab[use_idx], levels = lvls),
    hover = sprintf("comp size %d | final %s | lambda %.3f",
                    csz[use_idx], ifelse(ref[use_idx] == 0, "noise", paste0("C", ref[use_idx])),
                    1 / eps_grid[f]))
}
anim_df <- do.call(rbind, frames)

anim <- plot_ly(anim_df, x = ~UMAP1, y = ~UMAP2, frame = ~frame,
                color = ~color_label, colors = pal, type = "scattergl", mode = "markers",
                marker = list(size = 5, opacity = 0.75, line = list(width = 0)),
                text = ~hover, hoverinfo = "text") %>%
  layout(title = list(text = sprintf("HDBSCAN pruning on 5D UMAP (drawn on 2D, minPts = %d)", MIN_PTS)),
         xaxis = list(title = "UMAP1", zeroline = FALSE),
         yaxis = list(title = "UMAP2", zeroline = FALSE, scaleanchor = "x"),
         legend = list(title = list(text = "plurality label"))) %>%
  animation_opts(frame = 130, transition = 0, redraw = TRUE) %>%
  animation_slider(currentvalue = list(prefix = "frame (low lambda -> high lambda): "))

htmlwidgets::saveWidget(anim,
                        file.path(DOCS_DIR, paste0("pruning_animation_", RUN_TAG, ".html")),
                        selfcontained = TRUE, title = "HDBSCAN pruning animation")

# ── 9. CONSOLIDATED SAVES ──────────────────────────────────────────────
saveRDS(mst_edges, file.path(OBJ_DIR, paste0("mst_edges_",       RUN_TAG, ".rds")))
saveRDS(events,    file.path(OBJ_DIR, paste0("mst_events_",      RUN_TAG, ".rds")))
saveRDS(nodes_df,  file.path(OBJ_DIR, paste0("condensed_nodes_", RUN_TAG, ".rds")))
saveRDS(data.frame(id = seq_len(n), lam_death = lam_fallout, ref = ref, lab_rec = lab_rec),
        file.path(OBJ_DIR, paste0("point_labels_", RUN_TAG, ".rds")))

message("Done. HTML and plots in results/figures/06_robustness/hdbscan, objects in results/objects/HDBSCAN.")