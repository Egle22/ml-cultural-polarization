# =============================================================================
# scripts/06_robustness/05_information_vs_k.R
# Purpose : How much information does one more cluster really add?
#   (A) STRUCTURAL novelty : is C_k -> C_{k+1} a refinement or a reorganisation?
#   (B) SUBSTANTIVE novelty: does the top-10 of the importances recruit new
#                            variables / themes?
#   (C) decision rule      : the smallest k beyond which (A) and (B) are flat.
#   Functions only; the usage example at the end is wrapped in if (FALSE).
# Thesis  : Section 4.2.1 (choice of k = 3, supporting analysis)
# Input   : results/objects/clusters/all_clusters_aligned.rds (k = 2..6)
# Output  : console
# =============================================================================

source(here::here("R", "00_setup.R"))

## -----------------------------------------------------------------------------
## 0. Helpers: entropy, conditional entropy, Variation of Information
## -----------------------------------------------------------------------------
H <- function(x) {
  p <- as.numeric(table(x)) / length(x)
  p <- p[p > 0]
  -sum(p * log(p))
}

## H(a | b): residual uncertainty about a once b is known
H_cond <- function(a, b) {
  tb <- table(a, b) / length(a)
  pb <- colSums(tb)
  s <- 0
  for (j in seq_along(pb)) {
    if (pb[j] <= 0) next
    p <- tb[, j] / pb[j]
    p <- p[p > 0]
    s <- s + pb[j] * (-sum(p * log(p)))
  }
  s
}

## Meila (2007): VI = H(a|b) + H(b|a). A true metric on the space of partitions.
VI  <- function(a, b) H_cond(a, b) + H_cond(b, a)
NVI <- function(a, b) VI(a, b) / log(length(a))   # normalised to [0,1]

## Refinement index: each NEW cluster is mapped to its dominant OLD cluster;
## returns the share of units consistent with that mapping.
## = 1  -> pure split (no reorganisation, no new information)
## < 1  -> units are reshuffled across old clusters
refinement_index <- function(old, new) {
  tb <- table(new, old)
  sum(apply(tb, 1, max)) / sum(tb)
}

## -----------------------------------------------------------------------------
## 1. (A) STRUCTURAL NOVELTY along the sequence of k
## -----------------------------------------------------------------------------
## cl_list: list of label vectors (same length n), names = k
structural_seq <- function(cl_list) {
  ks <- as.integer(names(cl_list))
  ord <- order(ks); cl_list <- cl_list[ord]; ks <- ks[ord]

  do.call(rbind, lapply(seq_len(length(ks) - 1L), function(i) {
    a <- cl_list[[i]]; b <- cl_list[[i + 1L]]
    data.frame(
      k_from            = ks[i],
      k_to              = ks[i + 1L],
      VI                = VI(a, b),
      NVI               = NVI(a, b),
      H_old_given_new   = H_cond(a, b),   # ~0  => pure refinement
      H_new_given_old   = H_cond(b, a),   # always grows: it is the split
      refinement        = refinement_index(a, b),
      stringsAsFactors  = FALSE
    )
  }))
}

## -----------------------------------------------------------------------------
## 2. Bootstrap stability (Hennig): a k should be chosen only if clusters hold
## -----------------------------------------------------------------------------
## cluster_fun(idx) must re-cluster the subsample idx and return the labels.
## Example for PAM/Gower at the end of the file.
boot_stability <- function(labels, cluster_fun, B = 50, seed = 42) {
  set.seed(seed)
  n <- length(labels); K <- length(unique(labels))
  jac <- matrix(NA_real_, nrow = B, ncol = K)

  for (b in seq_len(B)) {
    idx <- sample.int(n, n, replace = TRUE)
    idx <- unique(idx)
    new <- cluster_fun(idx)
    old <- labels[idx]
    for (k in seq_len(K)) {
      A <- which(old == sort(unique(labels))[k])
      if (!length(A)) next
      jac[b, k] <- max(vapply(sort(unique(new)), function(g) {
        B_ <- which(new == g)
        length(intersect(A, B_)) / length(union(A, B_))
      }, numeric(1)))
    }
  }
  data.frame(cluster = seq_len(K), jaccard_mean = colMeans(jac, na.rm = TRUE))
}

## -----------------------------------------------------------------------------
## 3. (B) SUBSTANTIVE NOVELTY: top-N of the importances
## -----------------------------------------------------------------------------
## imp_list: list of NAMED numeric vectors (names = Q182, Q164, ...), names = k

top_vars <- function(imp, n = 10) {
  n <- min(n, length(imp))
  names(sort(imp, decreasing = TRUE))[seq_len(n)]
}

jaccard <- function(a, b) length(intersect(a, b)) / length(union(a, b))

## Truncated Rank-Biased Overlap: weights the top positions more and avoids the
## arbitrary top-10 cut. p ~ 0.9 => the first ~10 positions weigh 86%.
rbo <- function(a, b, p = 0.9) {
  d <- min(length(a), length(b))
  if (d == 0) return(NA_real_)
  s <- 0
  for (i in seq_len(d)) {
    s <- s + (length(intersect(a[seq_len(i)], b[seq_len(i)])) / i) * p^(i - 1)
  }
  ((1 - p) * s) / (1 - p^d)   # normalised: 1 = identical lists
}

## Saturation curve: how many DISTINCT variables have been seen up to k
saturation_curve <- function(top_list) {
  cum <- Reduce(union, top_list, accumulate = TRUE)
  data.frame(k = as.integer(names(top_list)),
             n_vars_cumulative = vapply(cum, length, integer(1)))
}

## -----------------------------------------------------------------------------
## 4. Map of WVS themes (fixed BEFORE looking at the results)
## -----------------------------------------------------------------------------
theme_map <- c(
  ## religiosity
  Q6 = "religion", Q15 = "religion", Q64 = "religion", Q94 = "religion",
  Q164 = "religion", Q167 = "religion", Q168 = "religion", Q169 = "religion",
  Q170 = "religion", Q171 = "religion", Q172 = "religion",
  ## science vs faith
  Q160 = "science", Q161 = "science",
  ## sexual morality / family
  Q36 = "sexual_morality", Q182 = "sexual_morality", Q184 = "sexual_morality",
  Q186 = "sexual_morality", Q188 = "sexual_morality", Q193 = "sexual_morality",
  ## civic morality
  Q179 = "civic_morality", Q181 = "civic_morality", Q189 = "civic_morality",
  Q192 = "civic_morality", Q195 = "civic_morality",
  ## immigration
  Q121 = "immigration", Q122 = "immigration", Q126 = "immigration",
  ## confidence in national institutions
  Q66 = "inst_national", Q69 = "inst_national", Q71 = "inst_national",
  Q76 = "inst_national", Q79 = "inst_national", Q80 = "inst_national",
  ## confidence in international institutions
  Q82 = "inst_intl", Q84 = "inst_intl", Q85 = "inst_intl", Q86 = "inst_intl",
  Q87 = "inst_intl", Q88 = "inst_intl", Q89 = "inst_intl",
  ## security / neighbourhood
  Q134 = "security", Q136 = "security", Q137 = "security",
  ## economy / role of the state
  Q54 = "economy", Q108 = "economy", Q111 = "economy"
)

theme_of <- function(v) {
  out <- unname(theme_map[v])
  out[is.na(out)] <- "other"
  out
}

## Thematic entropy of the top-N: high = the clustering draws on different domains
theme_entropy <- function(vars) H(theme_of(vars))

substantive_seq <- function(imp_list, n = 10) {
  ks <- as.integer(names(imp_list))
  ord <- order(ks); imp_list <- imp_list[ord]; ks <- ks[ord]

  tops <- lapply(imp_list, top_vars, n = n)
  names(tops) <- ks
  full <- lapply(imp_list, function(x) names(sort(x, decreasing = TRUE)))

  per_k <- data.frame(
    k             = ks,
    theme_entropy = vapply(tops, theme_entropy, numeric(1)),
    n_themes      = vapply(tops, function(v) length(unique(theme_of(v))), integer(1)),
    stringsAsFactors = FALSE
  )
  per_k <- merge(per_k, saturation_curve(tops), by = "k")

  trans <- do.call(rbind, lapply(seq_len(length(ks) - 1L), function(i) {
    data.frame(k_from = ks[i], k_to = ks[i + 1L],
               jaccard_top = jaccard(tops[[i]], tops[[i + 1L]]),
               rbo         = rbo(full[[i]], full[[i + 1L]]),
               new_vars    = length(setdiff(tops[[i + 1L]], tops[[i]])),
               new_themes  = length(setdiff(theme_of(tops[[i + 1L]]),
                                            theme_of(tops[[i]]))),
               stringsAsFactors = FALSE)
  }))

  list(per_k = per_k, transitions = trans)
}

## -----------------------------------------------------------------------------
## 5. (C) Decision rule
## -----------------------------------------------------------------------------
## The chosen k = the smallest k such that ALL following transitions are
##   (i)  almost pure refinements    (refinement >= ref_thr)
##   (ii) without thematic recruitment (new_themes == 0)
choose_k <- function(struct, subst, ref_thr = 0.95) {
  m <- merge(struct[, c("k_from", "k_to", "refinement")],
             subst$transitions[, c("k_from", "k_to", "new_themes", "new_vars")],
             by = c("k_from", "k_to"))
  m <- m[order(m$k_from), ]
  m$adds_nothing <- m$refinement >= ref_thr & m$new_themes == 0

  sel <- NA_integer_
  for (i in seq_along(m$k_from)) {
    if (all(m$adds_nothing[i:nrow(m)])) { sel <- m$k_from[i]; break }
  }
  list(table = m, k_selected = sel)
}

## -----------------------------------------------------------------------------
## 6. USAGE — adapt to your own objects
## -----------------------------------------------------------------------------
if (FALSE) {

  all_cl <- readRDS(PATH_CLUSTERS)

  ## labels for one method (the canonical file keeps k = 2..6)
  ks      <- 2:6
  cl_list <- setNames(lapply(ks, function(k) all_cl[[paste0("pam_", k)]]), ks)

  ## importances: list of named vectors, one per k
  ## imp_list <- readRDS(...)

  struct <- structural_seq(cl_list)
  subst  <- substantive_seq(imp_list, n = 10)
  dec    <- choose_k(struct, subst)

  print(struct); print(subst$per_k); print(dec$k_selected)

  ## bootstrap stability for PAM/Gower at the chosen k
  ## gd  <- readRDS(file.path(DIR_CLUST, "gower_dist.rds"))
  ## fun <- function(idx) cluster::pam(as.dist(as.matrix(gd)[idx, idx]),
  ##                                   k = dec$k_selected, diss = TRUE)$clustering
  ## boot_stability(cl_list[[as.character(dec$k_selected)]], fun, B = 50)
}
