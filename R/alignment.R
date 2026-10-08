# =============================================================================
# R/alignment.R
# Helpers to compare and align partitions across clustering configurations:
# Hungarian label alignment, global and per-cluster agreement, Cramér's V,
# and the pairwise concordance tables used in 02_clustering/04_concordance.R
# and 07_align_clusters.R (Table 4.2).
# =============================================================================

# Column names of the partitions with k clusters, one per configuration.
get_col_names <- function(k_val) {
  c(
    paste0("km_umap_",   k_val),
    paste0("pam_",       k_val),
    paste0("km_pc2_k",   k_val),
    paste0("km_pc5_k",   k_val),
    paste0("km_pc15_k",  k_val),
    paste0("km_pc75_k",  k_val),
    paste0("km_pc198_k", k_val)
  )
}

method_labels <- c(
  "km_umap"  = "KM-UMAP",
  "pam"      = "PAM",
  "km_pc2"   = "KM-PC2",
  "km_pc5"   = "KM-PC5",
  "km_pc15"  = "KM-PC15",
  "km_pc75"  = "KM-PC75",
  "km_pc198" = "KM-PC198"
)

# Drop the "_k" suffix: "km_pc5_k3" -> "km_pc5", "pam_3" -> "pam".
strip_k <- function(col) {
  col <- gsub("km_pc([0-9]+)_k[0-9]+", "km_pc\\1", col)
  col <- gsub("_k?[0-9]+$", "", col)
  col
}

cramers_v <- function(x, y) {
  tbl  <- table(x, y)
  chi2 <- suppressWarnings(chisq.test(tbl, correct = FALSE)$statistic)
  n    <- sum(tbl)
  k    <- min(nrow(tbl), ncol(tbl))
  as.numeric(sqrt(chi2 / (n * (k - 1))))
}

# Relabel `oth` onto `ref` with the Hungarian algorithm (maximum overlap).
hungarian_align <- function(ref, oth, k_val) {
  tab <- table(
    factor(ref, levels = as.character(seq_len(k_val))),
    factor(oth, levels = as.character(seq_len(k_val)))
  )
  asg <- clue::solve_LSAP(tab, maximum = TRUE)
  map <- setNames(as.character(seq_len(k_val)), as.character(asg))
  as.character(map[as.character(oth)])
}

# Share of respondents with the same label after Hungarian alignment.
global_agreement <- function(x, y, k_val) {
  y_al <- hungarian_align(x, y, k_val)
  mean(as.character(x) == as.character(y_al))
}

# For each cluster of x, share of its members that y (aligned) puts in the
# same cluster.
cluster_agreement <- function(x, y, col_i, col_j, k_val) {
  y_al  <- hungarian_align(x, y, k_val)
  x_chr <- as.character(x)
  y_chr <- as.character(y_al)
  clusts <- sort(unique(x_chr))
  do.call(rbind, lapply(clusts, function(cl) {
    idx <- x_chr == cl
    data.frame(
      method_i  = col_i,
      method_j  = col_j,
      cluster   = cl,
      n         = sum(idx),
      agreement = mean(y_chr[idx] == cl),
      stringsAsFactors = FALSE
    )
  }))
}

# Pairwise global agreement, Cramér's V and per-cluster agreement between all
# configurations, for every k. Returns list(global_long, clust_long).
compute_concordance <- function(clusters, k_values) {
  results_global <- list()
  results_clust  <- list()

  for (k_val in k_values) {
    cols <- get_col_names(k_val)
    cols <- cols[cols %in% colnames(clusters)]
    n_m  <- length(cols)
    if (n_m < 2) next

    mat_agr    <- matrix(NA, n_m, n_m, dimnames = list(cols, cols))
    mat_crv    <- matrix(NA, n_m, n_m, dimnames = list(cols, cols))
    clust_rows <- list()

    for (i in seq_len(n_m)) {
      for (j in seq_len(n_m)) {
        xi <- clusters[[cols[i]]]
        xj <- clusters[[cols[j]]]

        if (i == j) {
          mat_agr[i, j] <- 1
          mat_crv[i, j] <- 1
          next
        }

        mat_agr[i, j] <- global_agreement(xi, xj, k_val)
        mat_crv[i, j] <- cramers_v(xi, xj)

        if (i < j) {
          clust_rows[[paste(cols[i], cols[j], sep = "__")]] <-
            cluster_agreement(xi, xj, cols[i], cols[j], k_val)
        }
      }
    }

    df_agr <- as.data.frame(as.table(mat_agr)) %>%
      setNames(c("method_i", "method_j", "Agreement"))
    df_crv <- as.data.frame(as.table(mat_crv)) %>%
      setNames(c("method_i", "method_j", "CramerV"))

    results_global[[paste0("k", k_val)]] <- df_agr %>%
      left_join(df_crv, by = c("method_i", "method_j")) %>%
      mutate(
        k         = k_val,
        label_i   = method_labels[strip_k(as.character(method_i))],
        label_j   = method_labels[strip_k(as.character(method_j))],
        Agreement = round(Agreement, 3),
        CramerV   = round(CramerV,   3)
      )

    results_clust[[paste0("k", k_val)]] <- do.call(rbind, clust_rows) %>%
      mutate(
        k         = k_val,
        label_i   = method_labels[strip_k(method_i)],
        label_j   = method_labels[strip_k(method_j)],
        agreement = round(agreement, 3)
      )

    cat(sprintf("k = %d done\n", k_val))
  }

  list(
    global_long = do.call(rbind, results_global),
    clust_long  = do.call(rbind, results_clust)
  )
}
