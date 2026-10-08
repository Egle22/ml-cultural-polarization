# ═══════════════════════════════════════════════════════════════════════
#  scripts/06_robustness/08_hdbscan_peeling_maps.R
#
#  PART 1  the peeling + labels figure you liked, for d = 2, PLUS the two
#          diagnostics that answer "why do the central points die":
#            - a canvas coloured by d_core  (should be bright in the centre)
#            - a d_core vs lambda_death scatter (near-perfect inverse)
#          Rationale: d_core(x) is a FLOOR on every mutual-reachability edge
#          of x, so a high-d_core point has only heavy edges -> cut first ->
#          dies early. Central-in-canvas = low-density void = high d_core.
#
#  PART 2  the same peeling map, IMPROVED: clustering still on the 5D UMAP,
#          but points drawn in an interactive 3D UMAP space, coloured by
#          lambda_death (peeling) and by final label. Exported as self-
#          contained HTML (results/figures/06_robustness/hdbscan/html/).
#
#  Thesis  : Appendix A.4.1 (geometry of the unassigned respondents)
# ═══════════════════════════════════════════════════════════════════════

source(here::here("R", "00_setup.R"))

suppressPackageStartupMessages({
  library(dbscan)
  library(ggplot2)
  library(patchwork)
  library(plotly)
  library(htmlwidgets)
})

# ── CONFIG ─────────────────────────────────────────────────────────────
EMB2D_PATH   <- file.path(DIR_OBJ, "UMAP", "umap_2d.rds")   # PART 1
EMB5D_PATH   <- file.path(DIR_OBJ, "UMAP", "umap_5d.rds")   # PART 2 clustering
EMB3D_PATH   <- file.path(DIR_OBJ, "UMAP", "umap_3d.rds")   # PART 2 canvas
SRC_DATA     <- file.path(DIR_DATA, "data_imputed_selected.rds")  # only if 3D must be built
MIN_PTS      <- 25
SEED         <- 42

OBJ_DIR  <- ensure_dir(DIR_OBJ, "HDBSCAN")
PLOT_DIR <- ensure_dir(DIR_FIG, "06_robustness", "hdbscan")
DOCS_DIR <- ensure_dir(DIR_FIG, "06_robustness", "hdbscan", "html")
for (d in c(OBJ_DIR, PLOT_DIR, DOCS_DIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

noise_col <- "grey80"

# ── SHARED HELPER: lambda_death per point on any embedding ─────────────
# Returns per-point lambda_death (1 / weight of the edge that isolates the
# point's sub-minPts fragment), core distance, and the dbscan labels.
compute_peeling <- function(emb, minPts) {
  emb <- as.matrix(emb); storage.mode(emb) <- "double"
  n <- nrow(emb)
  knn    <- dbscan::kNN(emb, k = minPts - 1L)
  d_core <- knn$dist[, minPts - 1L]
  D  <- as.matrix(dist(emb))
  Dm <- pmax(D, outer(d_core, d_core, pmax)); diag(Dm) <- 0

  # Prim MST on mutual reachability
  in_tree <- logical(n); in_tree[1L] <- TRUE
  best_w  <- Dm[1L, ]; best_w[1L] <- Inf; best_from <- rep(1L, n)
  mst <- matrix(NA_real_, n - 1L, 3L, dimnames = list(NULL, c("from", "to", "weight")))
  for (i in seq_len(n - 1L)) {
    j <- which.min(ifelse(in_tree, Inf, best_w))
    mst[i, ] <- c(best_from[j], j, best_w[j]); in_tree[j] <- TRUE
    upd <- (!in_tree) & (Dm[j, ] < best_w); best_w[upd] <- Dm[j, upd]; best_from[upd] <- j
  }

  # replay merges in increasing weight = pruning in decreasing weight
  ord <- order(mst[, "weight"])
  parent <- seq_len(n); csize <- rep(1L, n)
  members <- split(seq_len(n), seq_len(n)); lam_death <- rep(NA_real_, n)
  find <- function(i) { while (parent[i] != i) i <- parent[i]; i }
  for (s in ord) {
    w <- mst[s, "weight"]; lam <- if (w > 0) 1 / w else Inf
    ra <- find(as.integer(mst[s, "from"])); rb <- find(as.integer(mst[s, "to"]))
    a <- csize[ra]; b <- csize[rb]
    if (a + b >= minPts) {                     # a point dies when its <minPts
      if (a < minPts) lam_death[members[[ra]]] <- lam   # fragment joins a
      if (b < minPts) lam_death[members[[rb]]] <- lam   # >=minPts component
    }
    if (a < b) { t <- ra; ra <- rb; rb <- t }
    parent[rb] <- ra; csize[ra] <- a + b
    members[[ra]] <- c(members[[ra]], members[[rb]]); members[[rb]] <- integer(0)
  }
  stopifnot(!anyNA(lam_death))

  hdb <- dbscan::hdbscan(emb, minPts = minPts)
  list(lam_death = lam_death, d_core = d_core,
       ref = hdb$cluster, K = max(hdb$cluster), prob = hdb$membership_prob)
}

label_scale <- function(K) scale_colour_manual(
  values = c("0" = noise_col, setNames(scales::hue_pal()(K), as.character(seq_len(K)))),
  name = "cluster")

thm <- theme_minimal(base_size = 10) + theme(panel.grid.minor = element_blank())

# ═══════════════════════════════════════════════════════════════════════
#  PART 1 — d = 2 : peeling + labels + the "why the centre dies" diagnostics
# ═══════════════════════════════════════════════════════════════════════
stopifnot(file.exists(EMB2D_PATH))
XY <- as.matrix(readRDS(EMB2D_PATH))[, 1:2, drop = FALSE]
colnames(XY) <- c("UMAP1", "UMAP2")
pe2 <- compute_peeling(XY, MIN_PTS)

df2 <- data.frame(UMAP1 = XY[, 1], UMAP2 = XY[, 2],
                  lam_death = pe2$lam_death, d_core = pe2$d_core,
                  ref = factor(pe2$ref))

rho_dcore <- cor(pe2$d_core, pe2$lam_death, method = "spearman")
cat(sprintf("d = 2 : %d clusters, %.1f%% noise | Spearman(d_core, lambda_death) = %+.3f\n",
            pe2$K, 100 * mean(pe2$ref == 0), rho_dcore))

g1 <- ggplot(df2, aes(UMAP1, UMAP2, colour = log10(lam_death))) +
  geom_point(size = 0.6) + scale_colour_viridis_c(option = "magma") + coord_equal() +
  labs(title = "peeling: lambda_death per point",
       subtitle = "dark = dies early (low lambda_death)",
       colour = expression(log[10](lambda[death]))) + thm

g2 <- ggplot(df2, aes(UMAP1, UMAP2, colour = ref)) +
  geom_point(size = 0.6) + coord_equal() + label_scale(pe2$K) +
  labs(title = "final HDBSCAN labels (0 = noise)") + thm

# THE DIAGNOSTIC: d_core map. Centre should be BRIGHT here (high d_core),
# i.e. the inverse of the dark centre in g1 -> that is the whole explanation.
g3 <- ggplot(df2, aes(UMAP1, UMAP2, colour = log10(d_core))) +
  geom_point(size = 0.6) + scale_colour_viridis_c(option = "viridis") + coord_equal() +
  labs(title = "d_core per point (the reason)",
       subtitle = "bright centre = sparse neighbourhood = heavy edges = dies first",
       colour = expression(log[10](d[core]))) + thm

g4 <- ggplot(df2, aes(d_core, lam_death,
                      colour = factor(ifelse(ref == 0, "noise", "clustered"),
                                      levels = c("noise", "clustered")))) +
  geom_point(size = 0.5, alpha = 0.4) + scale_x_log10() + scale_y_log10() +
  scale_colour_manual(values = c(noise = "grey60", clustered = "#B22222"), name = NULL) +
  labs(title = sprintf("d_core sets lambda_death  (rho = %+.3f)", rho_dcore),
       subtitle = "d_core is a floor on every reachability edge, so it almost fixes the death time",
       x = expression(d[core] ~ "(log)"), y = expression(lambda[death] ~ "(log)")) +
  thm + theme(legend.position = "bottom")

panel2d <- (g1 | g2) / (g3 | g4)
ggsave(file.path(PLOT_DIR, "peeling_maps_umap2_minpts25.pdf"), panel2d, width = 12, height = 10)
saveRDS(df2, file.path(OBJ_DIR, "peeling_df_umap2_minpts25.rds"))

# ═══════════════════════════════════════════════════════════════════════
#  PART 2 — clustering on 5D, drawn in an interactive 3D UMAP canvas
# ═══════════════════════════════════════════════════════════════════════
stopifnot(file.exists(EMB5D_PATH))
E5 <- as.matrix(readRDS(EMB5D_PATH))
stopifnot(ncol(E5) == 5)
pe5 <- compute_peeling(E5, MIN_PTS)                 # lambda_death from the 5D space
cat(sprintf("d = 5 : %d clusters, %.1f%% noise\n", pe5$K, 100 * mean(pe5$ref == 0)))

# 3D canvas: load if present, otherwise build it from the SAME source and
# UMAP config as your other embeddings (so it is comparable).
if (file.exists(EMB3D_PATH)) {
  E3 <- as.matrix(readRDS(EMB3D_PATH))
} else {
  message("umap_3d.rds not found -> building a 3D UMAP from source (one-off).")
  suppressPackageStartupMessages(library(umap))
  stopifnot(file.exists(SRC_DATA))
  Xsrc <- readRDS(SRC_DATA)
  Xsrc <- as.matrix(as.data.frame(lapply(Xsrc, function(c)
    if (is.factor(c)) as.numeric(as.character(c)) else as.numeric(c))))
  Xsrc <- scale(Xsrc)
  cfg <- umap::umap.defaults
  cfg$n_neighbors <- 15; cfg$min_dist <- 0.1; cfg$metric <- "euclidean"
  cfg$n_epochs <- 500; cfg$random_state <- SEED; cfg$n_components <- 3
  set.seed(SEED)
  E3 <- umap::umap(Xsrc, config = cfg)$layout
  saveRDS(E3, EMB3D_PATH)
}
stopifnot(nrow(E3) == nrow(E5), ncol(E3) >= 3)
colnames(E3) <- c("U1", "U2", "U3")

df3 <- data.frame(U1 = E3[, 1], U2 = E3[, 2], U3 = E3[, 3],
                  lam_death = pe5$lam_death, ref = pe5$ref)
df3$logld <- log10(df3$lam_death)
df3$hover <- sprintf("final %s | log10 lambda_death %.2f",
                     ifelse(df3$ref == 0, "noise", paste0("C", df3$ref)), df3$logld)

# continuous magma colourscale for plotly (base R, no extra package)
library(viridisLite)
mag <- viridis(16, option = "magma")
cs_magma <- lapply(seq_along(mag), function(i) list((i - 1) / (length(mag) - 1), mag[i]))

scene <- list(xaxis = list(title = "UMAP1"), yaxis = list(title = "UMAP2"),
              zaxis = list(title = "UMAP3"), aspectmode = "data")

# (a) 3D peeling — the improved version of the plot you liked
p3d_peel <- plot_ly(df3, x = ~U1, y = ~U2, z = ~U3, type = "scatter3d", mode = "markers",
                    text = ~hover, hoverinfo = "text",
                    marker = list(size = 2.6, opacity = 0.8, color = ~logld,
                                  colorscale = cs_magma, showscale = TRUE,
                                  colorbar = list(title = "log10(lambda_death)"))) %>%
  layout(title = "peeling in 3D (clustering on 5D UMAP, drawn on 3D UMAP)", scene = scene)

# (b) 3D final labels
K <- pe5$K
lab_pal <- c("0" = noise_col, setNames(grDevices::hcl.colors(K, "Dark 3"), as.character(seq_len(K))))
df3$ref_f <- factor(df3$ref, levels = c(0, seq_len(K)))
p3d_lab <- plot_ly(df3, x = ~U1, y = ~U2, z = ~U3, color = ~ref_f, colors = lab_pal,
                   type = "scatter3d", mode = "markers", text = ~hover, hoverinfo = "text",
                   marker = list(size = 2.6, opacity = 0.8)) %>%
  layout(title = "final HDBSCAN labels in 3D (0 = noise)", scene = scene,
         legend = list(title = list(text = "cluster")))

htmlwidgets::saveWidget(p3d_peel,
  file.path(DOCS_DIR, "peeling_3d_umap5_on_umap3_minpts25.html"),
  selfcontained = TRUE, title = "HDBSCAN peeling 3D")
htmlwidgets::saveWidget(p3d_lab,
  file.path(DOCS_DIR, "labels_3d_umap5_on_umap3_minpts25.html"),
  selfcontained = TRUE, title = "HDBSCAN labels 3D")
saveRDS(df3, file.path(OBJ_DIR, "peeling_df3d_umap5_minpts25.rds"))

message("Done. Figures and HTML in results/figures/06_robustness/hdbscan.")
