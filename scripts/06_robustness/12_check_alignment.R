## =============================================================================
## scripts/06_robustness/12_check_alignment.R
## Thesis  : quality check of Table 4.2 (not reported)
## Input   : results/objects/clusters/all_clusters_aligned.rds
## Output  : console
## =============================================================================
## Verify that the labels in all_clusters_aligned.rds are genuinely aligned
## across methods, using the Hungarian algorithm (clue::solve_LSAP).
##
## Logic: if labels are already optimally aligned, then re-running the
## Hungarian algorithm on any pair of methods must return the IDENTITY
## permutation. Any departure from the identity is a misalignment.
##
## ARI is computed alongside because it is invariant to relabelling: it
## separates "same partition, wrong labels" (high ARI, low diagonal
## agreement -> fixable) from "genuinely different partitions" (low ARI ->
## no relabelling can help).

rm(list = intersect(ls(), c("clusters", "labs", "res_pairs", "res_ref",
                            "M", "K", "obj_path")))

source(here::here("R", "00_setup.R"))   # clue, mclust (adjustedRandIndex)

#### ---- 0. Paths ---------------------------------------------------------

obj_path <- PATH_CLUSTERS
stopifnot(file.exists(obj_path))

clusters <- readRDS(obj_path)

#### ---- 1. Inspect structure (run once, then adjust get_labels below) ----

cat("\n--- structure of all_clusters_aligned ---\n")
str(clusters, max.level = 2)

#### ---- 2. Label extractor ----------------------------------------------
## ADJUST THIS FUNCTION to match the actual layout printed above.
## Three common layouts are handled; the stopifnot at the end fails loudly
## if none applies.

get_labels <- function(obj, method, k) {
  out <- NULL
  
  # layout A: nested list  obj[[method]][[as.character(k)]]
  if (is.list(obj) && !is.null(obj[[method]]) &&
      !is.null(obj[[method]][[as.character(k)]])) {
    out <- obj[[method]][[as.character(k)]]
  }
  
  # layout B: flat list with names like "KM-PC2_k5" or "KM-PC2_3"
  if (is.null(out)) {
    nm <- grep(paste0("^", method, "[_-]k?", k, "$"), names(obj), value = TRUE)
    if (length(nm) == 1L) out <- obj[[nm]]
  }
  
  # layout C: data.frame, one column per method x k
  if (is.null(out) && is.data.frame(obj)) {
    nm <- grep(paste0(method, ".*", k), names(obj), value = TRUE)
    if (length(nm) == 1L) out <- obj[[nm]]
  }
  
  stopifnot(!is.null(out))
  # normalise to a factor with levels "1".."k"
  out <- factor(as.character(as.integer(as.character(out))),
                levels = as.character(seq_len(k)))
  stopifnot(!anyNA(out))
  out
}

#### ---- 3. Configuration -------------------------------------------------

M <- c("KM-PC2", "KM-PC5", "KM-UMAP", "PAM")   # method keys as stored
K <- c(3, 5)                                    # k values to check

#### ---- 4. Core check: is the identity permutation optimal? --------------

check_pair <- function(a, b, k) {
  ## contingency table, methods a (rows) x b (cols)
  tab <- table(a, b)
  stopifnot(dim(tab) == c(k, k))
  
  ## Hungarian: maximise the diagonal (maximum = TRUE)
  perm <- as.integer(clue::solve_LSAP(tab, maximum = TRUE))
  
  n        <- length(a)
  agree_now <- sum(diag(tab)) / n                      # as stored
  agree_opt <- sum(tab[cbind(seq_len(k), perm)]) / n   # after optimal remap
  is_ident  <- identical(perm, seq_len(k))
  
  list(perm = perm,
       identity = is_ident,
       agree_now = agree_now,
       agree_opt = agree_opt,
       gain = agree_opt - agree_now,
       ari = mclust::adjustedRandIndex(a, b),
       tab = tab)
}

#### ---- 5. Pairwise sweep ------------------------------------------------

res_pairs <- do.call(rbind, lapply(K, function(k) {
  labs <- lapply(M, function(m) get_labels(clusters, m, k))
  names(labs) <- M
  stopifnot(length(unique(vapply(labs, length, 1L))) == 1L)
  
  cmb <- combn(M, 2, simplify = FALSE)
  do.call(rbind, lapply(cmb, function(p) {
    r <- check_pair(labs[[p[1]]], labs[[p[2]]], k)
    data.frame(k = k, a = p[1], b = p[2],
               identity  = r$identity,
               perm      = paste(r$perm, collapse = ""),
               agree_now = round(r$agree_now, 3),
               agree_opt = round(r$agree_opt, 3),
               gain      = round(r$gain, 3),
               ari       = round(r$ari, 3),
               stringsAsFactors = FALSE)
  }))
}))

cat("\n--- pairwise alignment check ---\n")
print(res_pairs, row.names = FALSE)

#### ---- 6. Reference-based check ----------------------------------------
## Alignment is usually done against ONE reference method. In that case
## reference-vs-others must be the identity; other pairs need not be.
## Set the reference actually used in the alignment script.

REF <- "PAM"

res_ref <- subset(res_pairs, a == REF | b == REF)
cat("\n--- alignment against reference:", REF, "---\n")
print(res_ref, row.names = FALSE)

if (all(res_ref$identity)) {
  cat("\nOK: labels are optimally aligned to", REF, "\n")
} else {
  cat("\nWARNING: alignment to", REF, "is NOT optimal for:\n")
  print(subset(res_ref, !identity), row.names = FALSE)
}

#### ---- 7. Confusion matrices for any failing pair ----------------------

bad <- subset(res_pairs, !identity)
if (nrow(bad) > 0) {
  for (i in seq_len(nrow(bad))) {
    k <- bad$k[i]
    a <- get_labels(clusters, bad$a[i], k)
    b <- get_labels(clusters, bad$b[i], k)
    r <- check_pair(a, b, k)
    cat("\n=== k =", k, ":", bad$a[i], "(rows) vs", bad$b[i], "(cols) ===\n")
    print(r$tab)
    cat("optimal permutation:", paste(r$perm, collapse = " "), "\n")
    cat("i.e. remap", bad$b[i], "labels:",
        paste(sprintf("%d->%d", r$perm, seq_len(k)), collapse = ", "), "\n")
    cat("agreement now:", round(r$agree_now, 3),
        "-> optimal:", round(r$agree_opt, 3),
        "| ARI:", round(r$ari, 3), "\n")
  }
} else {
  cat("\nAll pairs are already at the identity permutation.\n")
}

#### ---- 8. Interpretation guide -----------------------------------------

cat("\n--- how to read the output ---\n")
cat("identity = TRUE  and gain = 0        -> alignment is correct\n")
cat("identity = FALSE and gain > 0.02     -> real misalignment, remap needed\n")
cat("identity = FALSE and gain < 0.01     -> ties; partitions barely overlap\n")
cat("ARI < 0.30                           -> partitions genuinely differ;\n")
cat("                                        NO relabelling can fix this\n")