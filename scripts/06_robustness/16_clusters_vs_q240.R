################################################################################
# scripts/06_robustness/16_clusters_vs_q240.R
# Thesis  : not reported (clusters vs left-right self-placement Q240)
#
# For a given clustering method and a given k, produces the three-panel figure:
#   (1) PCA space, points coloured by CLUSTER label
#   (2) alluvial flow  CLUSTER  <->  Q240 group
#   (3) UMAP space, points coloured by Q240 group
# Discordant respondents (assigned to non-corresponding groups by the two
# partitions) are drawn in black in both scatterplots.
#
# Cluster labels are matched to Q240 groups by the Hungarian algorithm
# (clue::solve_LSAP) so that "discordant" is well defined: raw labels from two
# partitions are not comparable without an optimal one-to-one matching.
#
# NOTE: for k > 2 the matching of a data-driven partition onto an ORDINAL
#       grouping (left-right) is optimal in the agreement sense but has no
#       substantive ordering guarantee; state this in the caption.
################################################################################

#### 00. Setup #################################################################

rm(list = intersect(ls(), c("wvs", "clusters", "pca_scores", "umap_2d",
                            "fig_list")))

source(here::here("R", "00_setup.R"))

library(ggplot2)
library(ggalluvial)
library(patchwork)
library(clue)

obj_dir  <- DIR_OBJ
data_dir <- DIR_DATA
plot_dir <- ensure_dir(DIR_FIG, "06_robustness", "clusters_vs_q240")

#### 01. Config ################################################################

# Clustering columns to compare with Q240. Give the column name WITHOUT the k:
# the script appends "_k<k>" (and "_aligned" if that version exists).
METHODS <- c("KMeans_5pc", "KMeans_2pc", "PAM_Gower", "UMAP_KMeans")

K_VALUES <- 2:5

# Cut-points on the 1-10 left-right scale (upper bounds), one per k.
q240_schemes <- list(
  "2" = c(5, 10),
  "3" = c(3, 6, 10),
  "4" = c(3, 5, 7, 10),
  "5" = c(2, 4, 6, 8, 10)
)

#### 02. Data ##################################################################

wvs        <- readRDS(file.path(data_dir, "imputed_labeled.rds"))
clusters   <- readRDS(PATH_CLUSTERS)
pca_scores <- readRDS(file.path(obj_dir,  "pca_scores.rds"))
umap_2d    <- readRDS(file.path(obj_dir,  "UMAP", "umap_2d.rds"))   # adapt name

n <- nrow(wvs)

if ("id" %in% names(wvs) && "id" %in% names(clusters)) {
  clusters <- clusters[match(wvs$id, clusters$id), , drop = FALSE]
  stopifnot(identical(as.character(wvs$id), as.character(clusters$id)))
} else {
  stopifnot(nrow(clusters) == n)
}

pca_xy  <- as.data.frame(as.matrix(pca_scores)[, 1:2, drop = FALSE])
umap_xy <- as.data.frame(as.matrix(umap_2d)[, 1:2, drop = FALSE])
names(pca_xy)  <- c("PC1", "PC2")
names(umap_xy) <- c("UMAP1", "UMAP2")
stopifnot(nrow(pca_xy) == n, nrow(umap_xy) == n)

q240_num <- as.integer(as.character(wvs$Q240))
stopifnot(!anyNA(q240_num), all(q240_num %in% 1:10))

#### 03. Helpers ###############################################################

cut_q240 <- function(x, upper) {
  brk <- c(0, upper)
  lab <- paste0(head(brk, -1) + 1, "-", upper)
  cut(x, breaks = brk, labels = lab, include.lowest = TRUE, right = TRUE)
}

# Optimal one-to-one relabelling of g onto the levels of ref (max agreement).
align_to_ref <- function(g, ref) {
  g   <- droplevels(as.factor(g))
  ref <- droplevels(as.factor(ref))
  stopifnot(nlevels(g) == nlevels(ref))
  tab <- table(g, ref)
  map <- clue::solve_LSAP(as.matrix(tab), maximum = TRUE)
  factor(levels(ref)[as.integer(map)[as.integer(g)]], levels = levels(ref))
}

# Resolve a column name: prefers the *_aligned version when it exists.
resolve_col <- function(method, k) {
  cand <- c(paste0(method, "_k", k, "_aligned"), paste0(method, "_k", k))
  hit  <- cand[cand %in% names(clusters)]
  if (!length(hit)) return(NA_character_)
  hit[1]
}

pretty_map <- c("KMeans_2pc" = "KMeans-2PC", "KMeans_5pc" = "KMeans-5PC",
                "KMeans_15pc" = "KMeans-15PC", "KMeans_75pc" = "KMeans-75PC",
                "KMeans_198pc" = "KMeans-198PC", "PAM_Gower" = "PAM-Gower",
                "UMAP_KMeans" = "KM-UMAP",
                "km_pc2" = "KMeans-2PC", "km_pc5" = "KMeans-5PC",
                "km_umap" = "KM-UMAP")
pretty_lab <- function(x) if (x %in% names(pretty_map)) pretty_map[[x]] else gsub("_", "-", x)

# Palette shared by the two scatterplots so that colours mean the same thing.
grp_palette <- function(levs) {
  base <- c("#4E79A7", "#C1705C", "#8CA9C4", "#E0A490", "#3C5A73")
  setNames(base[seq_along(levs)], levs)
}

#### 04. Figure builder ########################################################

make_triptych <- function(method, k) {
  
  col <- resolve_col(method, k)
  if (is.na(col)) {
    message("skip: ", method, " k = ", k, " not available")
    return(NULL)
  }
  
  g_q   <- cut_q240(q240_num, q240_schemes[[as.character(k)]])
  g_clu <- align_to_ref(clusters[[col]], g_q)
  
  discord   <- g_clu != g_q
  pct_disc  <- 100 * mean(discord)
  meth_lab  <- pretty_lab(method)
  pal       <- grp_palette(levels(g_q))
  
  df <- data.frame(pca_xy, umap_xy,
                   cluster = g_clu, q240 = g_q, discord = discord)
  
  sub_txt <- sprintf("Black = discordant with Q240 groups (%.1f%% of obs)", pct_disc)
  
  base_theme <- theme_minimal(base_size = 9) +
    theme(plot.title    = element_text(size = 9, face = "bold"),
          plot.subtitle = element_text(size = 7, colour = "grey40"),
          panel.grid.minor = element_blank(),
          legend.position  = "right",
          legend.key.size  = unit(0.7, "lines"))
  
  # --- panel 1: PCA space, cluster labels -------------------------------------
  p1 <- ggplot(df, aes(PC1, PC2)) +
    geom_point(data = subset(df, !discord), aes(colour = cluster),
               size = 0.5, alpha = 0.55) +
    geom_point(data = subset(df, discord), colour = "black",
               size = 0.5, alpha = 0.8) +
    scale_colour_manual(values = pal, name = paste0(meth_lab, "\ncluster")) +
    labs(title = sprintf("PCA space - %s labels (k = %d)", meth_lab, k),
         subtitle = sub_txt, x = "PC 1", y = "PC 2") +
    base_theme +
    guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1)))
  
  # --- panel 2: alluvial cluster <-> Q240 -------------------------------------
  alluv <- as.data.frame(table(cluster = df$cluster, q240 = df$q240),
                         stringsAsFactors = TRUE)
  alluv <- subset(alluv, Freq > 0)
  
  p2 <- ggplot(alluv, aes(y = Freq, axis1 = cluster, axis2 = q240)) +
    geom_alluvium(aes(fill = cluster), width = 1/6, alpha = 0.55,
                  knot.pos = 0.3) +
    geom_stratum(aes(fill = after_stat(stratum)), width = 1/6,
                 colour = "white", linewidth = 0.3) +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)),
              size = 2.6, colour = "white", fontface = "bold") +
    scale_x_discrete(limits = c(meth_lab, "Q240"), expand = c(0.1, 0.1)) +
    scale_fill_manual(values = pal, guide = "none") +
    labs(title = sprintf("Cluster flow: %s <-> Q240 (k = %d)", meth_lab, k),
         y = "Observations", x = NULL) +
    base_theme
  
  # --- panel 3: UMAP space, Q240 groups ---------------------------------------
  p3 <- ggplot(df, aes(UMAP1, UMAP2)) +
    geom_point(data = subset(df, !discord), aes(colour = q240),
               size = 0.5, alpha = 0.55) +
    geom_point(data = subset(df, discord), colour = "black",
               size = 0.5, alpha = 0.8) +
    scale_colour_manual(values = pal, name = "Q240\ngroup") +
    labs(title = sprintf("UMAP space - Q240 groups (k = %d)", k),
         subtitle = sub_txt, x = "UMAP 1", y = "UMAP 2") +
    base_theme +
    guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1)))
  
  # --- assemble ---------------------------------------------------------------
  p <- (p1 | p2 | p3) +
    plot_annotation(
      title   = sprintf("%s vs Q240 (left-right self-placement) - k = %d",
                        meth_lab, k),
      caption = sprintf("Reference: Q240 groups | Hungarian alignment | %.1f%% discordant",
                        pct_disc),
      theme   = theme(plot.title   = element_text(face = "bold", size = 11),
                      plot.caption = element_text(size = 7, colour = "grey40"))
    )
  
  attr(p, "pct_disc") <- pct_disc
  p
}

#### 05. Build and save ########################################################

fig_list <- list()

for (m in METHODS) {
  for (k in K_VALUES) {
    p <- make_triptych(m, k)
    if (is.null(p)) next
    nm <- paste0(m, "_k", k)
    fig_list[[nm]] <- p
    ggsave(file.path(plot_dir, paste0("triptych_", nm, "_vs_Q240.pdf")),
           p, width = 11, height = 4, device = cairo_pdf)
    message("saved: ", nm, "  (", sprintf("%.1f", attr(p, "pct_disc")),
            "% discordant)")
  }
}

# quick look at one of them
print(fig_list[[1]])

