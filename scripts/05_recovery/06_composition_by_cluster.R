# =============================================================================
# scripts/05_recovery/06_composition_by_cluster.R
# Purpose : Sociodemographic composition of the three clusters (KM-PC5, k = 3):
#           (a) distribution of each category within each cluster;
#           (b) for each category, the share belonging to each cluster.
#           Variables: sex (Q260), occupational group (Q281), subjective class
#           (Q287), income decile (Q288), religious denomination (Q289),
#           ethnic group (Q290), age group (from Q262).
# Thesis  : Section 4.3.6; Figure 4.9
# Input   : results/objects/clusters/all_clusters_aligned.rds
#           data/processed/imputed_labeled.rds
# Output  : results/figures/05_recovery/composition/<cluster_col>/
#           {within_cluster,category_share}/<variable>.png
# Note    : the short cluster names refer to the labels aligned to KM-PC5.
#           U.S. reference statistics: https://datausa.io/profile/geo/united-states
# =============================================================================

source(here::here("R", "00_setup.R"))

# -- Settings --------------------------------------------------------------------
all_clusters_aligned <- readRDS(PATH_CLUSTERS)
imputed_labeled      <- readRDS(file.path(DIR_DATA, "imputed_labeled.rds"))

cluster_col <- "km_pc5_k3"   # Figure 4.9: KM-PC5, k = 3
imputed_labeled[[cluster_col]] <- all_clusters_aligned[[cluster_col]]

output_root <- ensure_dir(DIR_FIG, "05_recovery", "composition")

# -- Short cluster labels ----------------------------------------------------------
short_labels_k3 <- c(
  "1" = "Precarious and Disaffected",
  "2" = "Secular and Globally Oriented",
  "3" = "Religious and Order-Oriented"
)
short_labels_k4 <- c(
  "1" = "Vulnerable Traditionalists",
  "2" = "Faithful Modernists",
  "3" = "Secular Liberals",
  "4" = "Devout Conservatives"
)
short_labels_k5 <- c(
  "1" = "Vulnerable Traditionalists",
  "2" = "Religious Progressives",
  "3" = "Religious Conservatives",
  "4" = "Secular Nationalists",
  "5" = "Secular Liberals"
)

# -- Age groups (from Q262) --------------------------------------------------------
imputed_labeled <- imputed_labeled %>%
  mutate(
    Q262_num  = as.numeric(as.character(Q262)),
    age_group = cut(
      Q262_num,
      breaks         = c(0, 25, 35, 45, 55, 65,Inf),
      labels         = c("0-25", "25-35", "35-45", "45-55", "55-65", "65Plus."),
      right          = FALSE,
      include.lowest = TRUE
    )
  )

table(imputed_labeled$age_group)
round(table(imputed_labeled$Q290)/nrow(imputed_labeled),2)

# -- (a) One bar chart per cluster (facets), shares within the cluster -------------
plot_cluster_distribution <- function(data, var, cluster_col, title,
                                      cluster_labels = NULL) {
  df <- data %>%
    dplyr::filter(!is.na(.data[[var]]), !is.na(.data[[cluster_col]])) %>%
    dplyr::mutate(
      cluster = factor(.data[[cluster_col]]),
      cat     = factor(.data[[var]])
    ) %>%
    dplyr::count(cluster, cat) %>%
    dplyr::group_by(cluster) %>%
    dplyr::mutate(prop = n / sum(n)) %>%
    dplyr::ungroup()

  # labeller: use the label vector if given, otherwise label_both
  facet_labeller <- if (!is.null(cluster_labels)) {
    as_labeller(cluster_labels)
  } else {
    label_both
  }

  ggplot(df, aes(x = cat, y = prop, fill = cat)) +
    geom_col(color = "white", linewidth = 0.2) +
    facet_wrap(~ cluster, labeller = facet_labeller) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(title = title, x = var, y = "Proportion within cluster",
         fill = var) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title  = element_text(face = "bold"),
      strip.text  = element_text(face = "bold")
    )
}

# -- (b) For each category of a variable, the share in each cluster ---------------
plot_category_cluster_share <- function(data, var, cluster_col, title,
                                        cluster_labels = NULL,
                                        cluster_colors = NULL,
                                        position = "fill") {
  df <- data %>%
    dplyr::filter(!is.na(.data[[var]]), !is.na(.data[[cluster_col]])) %>%
    dplyr::mutate(
      cluster_id = factor(.data[[cluster_col]]),   # numeric: needed for the colours
      cat        = factor(.data[[var]])
    ) %>%
    dplyr::count(cat, cluster_id) %>%
    dplyr::group_by(cat) %>%
    dplyr::mutate(prop = n / sum(n)) %>%
    dplyr::ungroup()

  # levels present in the data, in order
  lev <- levels(df$cluster_id)

  # colours restricted to the clusters present (matched on the NUMBERS)
  fill_values <- if (!is.null(cluster_colors)) cluster_colors[lev] else NULL

  # legend labels: short labels if given, otherwise the number
  fill_labels <- if (!is.null(cluster_labels)) cluster_labels[lev] else lev

  ggplot(df, aes(x = cat, y = prop, fill = cluster_id)) +
    geom_col(color = "white", linewidth = 0.2) +
    scale_fill_manual(values = fill_values,
                      labels = fill_labels,
                      name   = "Cluster") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(title = title, x = var, y = "Share of category") +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title  = element_text(face = "bold")
    )
}

# -- Save helper -------------------------------------------------------------------
# Folder structure:
#   output_root/<cluster_col>/within_cluster/<var>.png
#   output_root/<cluster_col>/category_share/<var>.png
save_plot <- function(plot, var, cluster_col, type,
                      root = output_root,
                      width = 9, height = 6, dpi = 300) {
  dir <- ensure_dir(root, cluster_col, type)
  file <- file.path(dir, paste0(var, ".png"))
  ggsave(file, plot = plot, width = width, height = height, dpi = dpi)
  message("Saved: ", file)
  invisible(file)
}

vars_within <- c("Q260","Q281","Q287","Q288","Q289","Q290","age_group")
titles_within <- c(
  Q260      = "Respondent's sex",
  Q281      = "Respondent's occupational group",
  Q287      = "Perceived social class",
  Q288      = "Household income decile (1 = bottom 10%, 10 = top 10%)",
  Q289      = "Religious denomination",
  Q290      = "Ethnic group",
  age_group = "Respondent age"
)

# -- Within-cluster distribution ---------------------------------------------------
for (v in vars_within) {
  p <- plot_cluster_distribution(imputed_labeled, v, cluster_col,
                                 title = titles_within[[v]],
                                 cluster_labels = short_labels_k3)
  print(p)
  save_plot(p, v, cluster_col, type = "within_cluster")
}

# -- Share of each category by cluster (Figure 4.9) --------------------------------
for (v in vars_within) {
  p <- plot_category_cluster_share(imputed_labeled, v, cluster_col,
                                   title = titles_within[[v]],
                                   cluster_labels = short_labels_k3,
                                   cluster_colors = cluster_colors)
  print(p)
  save_plot(p, v, cluster_col, type = "category_share")
}
