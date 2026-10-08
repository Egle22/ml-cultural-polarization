# =============================================================================
# scripts/05_recovery/01_recovery_k3.R
# Purpose : Out-of-sample recovery of cluster membership at k = 3 only, from
#           the held-out blocks. Three sets of runs, each over the three
#           predictor blocks (pol, dem, pol_dem; geography excluded):
#
#             FULL     KM-PC2, KM-PC5, KM-UMAP, PAM on all n = 2,579
#             HDB      HDBSCAN on the assigned respondents only (noise dropped)
#             MATCHED  the same four methods, restricted to exactly the
#                      respondents HDBSCAN assigns
#
#           MATCHED exists because HDBSCAN metrics computed on its own
#           assigned subset are NOT comparable with FULL: the subset is the
#           densest part of the embedding, where any classifier does better.
#           Comparing HDB against FULL would be a selection-bias artefact.
#           Compare HDB against MATCHED instead.
# Thesis  : Section 4.3; Tables 4.4, 4.6 (inputs), 4.7, 4.8; Figure 4.8 (via 02_recall_k3.R)
# Input   : results/objects/clusters/all_clusters_aligned.rds
#           data/processed/imputed_labeled.rds
# Output  : results/objects/recovery_k3/<block>/{metrics_k3,recall_k3,
#           rf_importance_k3,lasso_coefs_k3,term_map,run_meta}.rds
#           results/objects/recovery_k3/metrics_k3_all.rds
# Models  : multinomial logit (nnet), multinomial LASSO (cv.glmnet,
#           lambda.1se, type.measure = "class"), random forest (ranger,
#           300 trees, impurity importance); 5-fold stratified CV; SEED = 42.
# Note    : runs in parallel with future (multisession, future.seed = TRUE):
#           documented exception to the "sequential lapply" rule, kept so that
#           the thesis numbers are reproduced exactly (see CLAUDE.md).
#           Factors enter as unordered (ordered = FALSE); model-matrix columns
#           are mapped to variables with the "assign" attribute.
# =============================================================================

source(here::here("R", "00_setup.R"))

suppressMessages({
  library(dplyr);  library(tidyr)
  library(nnet);   library(glmnet); library(ranger)
  library(future); library(future.apply); library(progressr)
})

# ── 0. Parameters ────────────────────────────────────────────────────####

K         <- 3L
BLOCKS    <- c("pol", "dem", "pol_dem")
METHODS   <- c("KM-PC2", "KM-PC5", "KM-UMAP", "PAM")
NOISE_LAB <- c("0", "-1")          # HDBSCAN noise codes

NFOLDS    <- 5L
GLMNET_NF <- 5L
LAMBDA    <- "lambda.1se"
TYPE_MEAS <- "class"
N_TREES   <- 300L
TOP_N     <- 20L
SEED      <- 42L
RF_IMPORTANCE <- "impurity"     # "impurity" is ~8x faster
N_WORKERS <- max(1L, min(6L, parallel::detectCores(logical = FALSE) - 1L))

AGE_COL <- "Q262"; INCOME_COL <- "Q288"

out_root <- file.path(DIR_OBJ, "recovery_k3")
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)

# ── 1. Data ──────────────────────────────────────────────────────────####
all_cl   <- readRDS(PATH_CLUSTERS)
data_all <- readRDS(file.path(DIR_DATA, "imputed_labeled.rds"))
N <- nrow(all_cl); stopifnot(nrow(data_all) == N)
cat("N =", N, "\n")

col_for <- function(m) switch(m,
                              "KM-PC2" = "km_pc2_k3", "KM-PC5" = "km_pc5_k3",
                              "KM-UMAP" = "km_umap_3", "PAM" = "pam_3",
                              stop("unknown method: ", m))
stopifnot(all(vapply(METHODS, col_for, "") %in% names(all_cl)))

# HDBSCAN column: taken from the data, not hard-coded
hdb_col <- grep("^hdb", names(all_cl), value = TRUE)
if (length(hdb_col) != 1L)
  stop("expected exactly one hdb* column, found: ", paste(hdb_col, collapse = ", "))
cat("HDBSCAN column:", hdb_col, "\n")

hdb_raw    <- as.character(all_cl[[hdb_col]])
is_noise   <- hdb_raw %in% NOISE_LAB
idx_assign <- which(!is_noise)
n_assign   <- length(idx_assign)
cat(sprintf("HDBSCAN: %d assigned (%.1f%%), %d noise (%.1f%%), %d clusters\n",
            n_assign, 100 * n_assign / N, sum(is_noise), 100 * sum(is_noise) / N,
            length(unique(hdb_raw[!is_noise]))))

# ── 2. Design matrix for one predictor block ─────────────────────────####
build_X <- function(block) {
  cols <- switch(block,
                 "pol"     = Political,
                 "dem"     = Demographics,
                 "pol_dem" = c(Political, Demographics),
                 stop("unknown block: ", block))
  keep <- intersect(cols, colnames(data_all))
  d <- dplyr::select(data_all, dplyr::all_of(keep))
  
  bin_age <- function(x) cut(suppressWarnings(as.numeric(as.character(x))),
                             breaks = c(0,15,25,35,45,60,Inf), right = FALSE,
                             labels = c("0_15","16_25","26_34","35_44","45_59","60plus"))
  bin_inc <- function(x) cut(suppressWarnings(as.numeric(as.character(x))),
                             breaks = c(0,3,7,10), right = TRUE, labels = c("Low","Middle","High"))
  looks_num <- function(x) is.numeric(x) ||
    mean(is.na(suppressWarnings(as.numeric(as.character(x))))) < 0.5
  if (AGE_COL    %in% names(d) && looks_num(d[[AGE_COL]]))    d[[AGE_COL]]    <- bin_age(d[[AGE_COL]])
  if (INCOME_COL %in% names(d) && looks_num(d[[INCOME_COL]])) d[[INCOME_COL]] <- bin_inc(d[[INCOME_COL]])
  
  # ordered = FALSE: ordered factors would give polynomial contrasts and
  # silently invalidate the LASSO and logit coefficients.
  d <- d %>% mutate(across(everything(), ~ droplevels(factor(., ordered = FALSE))))
  stopifnot(!any(vapply(d, is.ordered, logical(1))), !anyNA(d))
  drop1 <- names(which(sapply(d, nlevels) < 2))
  if (length(drop1)) d <- dplyr::select(d, -all_of(drop1))
  
  mm <- model.matrix(~ ., data = d)
  stopifnot(nrow(mm) == N)
  X <- mm[, -1, drop = FALSE]
  tm <- data.frame(term_raw = colnames(mm)[-1],
                   variable = names(d)[attr(mm, "assign")[-1]],
                   stringsAsFactors = FALSE)
  colnames(X) <- make.names(colnames(X), unique = TRUE)
  tm$term <- colnames(X)
  tm <- left_join(tm, section_lookup, by = "variable") %>%
    mutate(section_label = ifelse(is.na(section_label), "Other", section_label))
  list(X = X, term_map = tm, n_vars = ncol(d))
}

# ── 3. Model helpers ─────────────────────────────────────────────────####
make_folds <- function(y, nfolds, seed) {
  set.seed(seed); fold <- integer(length(y))
  for (cl in unique(y)) {
    i <- which(y == cl); fold[i] <- sample(rep_len(seq_len(nfolds), length(i)))
  }
  fold
}

cv_all <- function(X, y, cfg) {
  folds <- make_folds(y, cfg$nfolds, cfg$seed)
  Xdf <- as.data.frame(X)
  oof <- list(Logit = character(length(y)), LASSO = character(length(y)),
              RF = character(length(y)))
  imp_rf <- list(); coef_la <- list()
  
  for (f in sort(unique(folds))) {
    tr <- folds != f; te <- folds == f
    
    fit <- nnet::multinom(.y ~ ., data = data.frame(.y = y[tr], Xdf[tr, , drop = FALSE]),
                          trace = FALSE, maxit = 500L, MaxNWts = 1e5)
    oof$Logit[te] <- as.character(predict(fit, newdata = Xdf[te, , drop = FALSE],
                                          type = "class"))
    
    set.seed(cfg$seed)
    cvf <- glmnet::cv.glmnet(X[tr, , drop = FALSE], y[tr], family = "multinomial",
                             alpha = 1, nfolds = cfg$glmnet_nf,
                             type.measure = cfg$type_measure)
    oof$LASSO[te] <- as.character(as.vector(
      predict(cvf, newx = X[te, , drop = FALSE], s = cfg$lambda, type = "class")))
    cl <- coef(cvf, s = cfg$lambda)
    coef_la[[as.character(f)]] <- bind_rows(lapply(names(cl), function(g) {
      m <- as.matrix(cl[[g]])
      data.frame(class = g, term = rownames(m), coef = as.numeric(m[, 1]),
                 stringsAsFactors = FALSE)
    }))
    
    set.seed(cfg$seed)
    rf <- ranger::ranger(x = Xdf[tr, , drop = FALSE], y = y[tr],
                         num.trees = cfg$n_trees, importance = cfg$rf_importance,
                         num.threads = 1L, classification = TRUE)
    oof$RF[te] <- as.character(predict(rf, data = Xdf[te, , drop = FALSE],
                                       num.threads = 1L)$predictions)
    imp_rf[[as.character(f)]] <- data.frame(
      term = names(rf$variable.importance),
      imp  = as.numeric(rf$variable.importance), stringsAsFactors = FALSE)
  }
  
  imp <- bind_rows(imp_rf) %>% group_by(term) %>%
    summarise(imp = mean(imp), .groups = "drop") %>% arrange(desc(imp))
  cf <- bind_rows(coef_la) %>% filter(term != "(Intercept)") %>%
    group_by(class, term) %>%
    summarise(mean_coef = mean(coef), nonzero_folds = sum(coef != 0),
              .groups = "drop") %>% filter(nonzero_folds > 0)
  
  list(oof = lapply(oof, factor, levels = levels(y)), imp_rf = imp, coef_lasso = cf)
}

eval_pred <- function(pred, truth) {
  lev <- sort(unique(c(as.character(truth), as.character(pred))))
  cm <- table(pred = factor(as.character(pred), levels = lev),
              truth = factor(as.character(truth), levels = lev))
  n <- sum(cm); k_eff <- sum(colSums(cm) > 0)
  p_o <- sum(diag(cm)) / n
  p_e <- sum(rowSums(cm) * colSums(cm)) / n^2
  data.frame(n = n, acc = p_o, k_eff = k_eff,
             kappa_emp  = (p_o - p_e) / (1 - p_e),
             kappa_unif = (p_o - 1/k_eff) / (1 - 1/k_eff))
}

per_class <- function(pred, truth) {
  lev <- sort(unique(as.character(truth)))
  cm <- table(pred = factor(as.character(pred), levels = lev),
              truth = factor(as.character(truth), levels = lev))
  tp <- diag(cm); supp <- colSums(cm); pn <- rowSums(cm)
  rec <- ifelse(supp > 0, tp/supp, NA_real_); prec <- ifelse(pn > 0, tp/pn, NA_real_)
  data.frame(cluster = names(tp), support = as.integer(supp),
             recall = as.numeric(rec), precision = as.numeric(prec),
             f1 = as.numeric(ifelse(prec + rec > 0, 2*prec*rec/(prec+rec), 0)),
             stringsAsFactors = FALSE)
}

CFG <- list(nfolds = NFOLDS, seed = SEED, lambda = LAMBDA,
            glmnet_nf = GLMNET_NF, type_measure = TYPE_MEAS,
            n_trees = N_TREES, rf_importance = RF_IMPORTANCE)

# ── 4. Task grid ─────────────────────────────────────────────────────####
# run_type: FULL / HDB / MATCHED. rows = which respondents enter the fit.
tasks <- bind_rows(
  data.frame(run_type = "FULL",    method = METHODS,   stringsAsFactors = FALSE),
  data.frame(run_type = "HDB",     method = "HDBSCAN", stringsAsFactors = FALSE),
  data.frame(run_type = "MATCHED", method = METHODS,   stringsAsFactors = FALSE)
)
tasks$tag <- paste(tasks$run_type, tasks$method, sep = "_")

label_for <- function(run_type, method) {
  if (method == "HDBSCAN") {
    y <- factor(hdb_raw[idx_assign]); rows <- idx_assign
  } else {
    lab <- as.character(all_cl[[col_for(method)]])
    rows <- if (run_type == "MATCHED") idx_assign else seq_len(N)
    y <- factor(lab[rows])
  }
  list(y = droplevels(y), rows = rows)
}

# ── 5. Run ───────────────────────────────────────────────────────────####
options(future.globals.maxSize = 2 * 1024^3)
handlers("txtprogressbar")
t_all <- Sys.time()

for (blk in BLOCKS) {
  cat("\n\n══════  BLOCK:", blk, " ══════\n")
  B <- build_X(blk)
  cat(sprintf("%d variables -> %d dummies\n", B$n_vars, ncol(B$X)))
  X_full <- B$X
  
  out_dir <- file.path(out_root, blk)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  plan(multisession, workers = N_WORKERS)
  res <- with_progress({
    p <- progressor(steps = nrow(tasks))
    future_lapply(seq_len(nrow(tasks)), function(i) {
      lb <- label_for(tasks$run_type[i], tasks$method[i])
      r  <- cv_all(X_full[lb$rows, , drop = FALSE], lb$y, CFG)
      p(tasks$tag[i])
      list(tag = tasks$tag[i], run_type = tasks$run_type[i],
           method = tasks$method[i], truth = lb$y, n = length(lb$rows),
           oof = r$oof, imp_rf = r$imp_rf, coef_lasso = r$coef_lasso)
    }, future.seed = TRUE,
    future.globals = list(
      tasks = tasks, X_full = X_full, CFG = CFG, all_cl = all_cl, N = N,
      hdb_raw = hdb_raw, idx_assign = idx_assign, col_for = col_for,
      label_for = label_for, cv_all = cv_all, make_folds = make_folds, p = p),
    future.packages = c("dplyr", "nnet", "glmnet", "ranger"))
  })
  plan(sequential)
  names(res) <- tasks$tag
  
  metrics <- bind_rows(lapply(res, function(r)
    bind_rows(lapply(names(r$oof), function(cl)
      data.frame(classifier = cl, eval_pred(r$oof[[cl]], r$truth)))) %>%
      mutate(block = blk, run_type = r$run_type, method = r$method, .before = 1)))
  
  recall <- bind_rows(lapply(res, function(r)
    bind_rows(lapply(names(r$oof), function(cl)
      data.frame(classifier = cl, per_class(r$oof[[cl]], r$truth)))) %>%
      mutate(block = blk, run_type = r$run_type, method = r$method, .before = 1)))
  
  imp <- bind_rows(lapply(res, function(r)
    r$imp_rf %>% mutate(block = blk, run_type = r$run_type, method = r$method))) %>%
    left_join(B$term_map[, c("term", "variable", "section_label")], by = "term")
  
  cfl <- bind_rows(lapply(res, function(r)
    if (nrow(r$coef_lasso) == 0) NULL else
      r$coef_lasso %>% mutate(block = blk, run_type = r$run_type,
                              method = r$method))) %>%
    { if (nrow(.) == 0) . else left_join(., B$term_map[, c("term","variable","section_label")],
                                         by = "term") }
  
  cat("\n──── k = 3 |", blk, "────\n")
  print(as.data.frame(metrics %>% filter(classifier == "LASSO") %>%
                        select(run_type, method, n, acc, kappa_emp)), row.names = FALSE, digits = 3)
  
  saveRDS(metrics, file.path(out_dir, "metrics_k3.rds"))
  saveRDS(recall,  file.path(out_dir, "recall_k3.rds"))
  saveRDS(imp,     file.path(out_dir, "rf_importance_k3.rds"))
  saveRDS(cfl,     file.path(out_dir, "lasso_coefs_k3.rds"))
  saveRDS(B$term_map, file.path(out_dir, "term_map.rds"))
  saveRDS(list(block = blk, k = K, n_full = N, n_assigned = n_assign,
               hdb_col = hdb_col, n_vars = B$n_vars, n_dummies = ncol(X_full),
               seed = SEED, nfolds = NFOLDS, lambda = LAMBDA,
               rf_importance = RF_IMPORTANCE, run_time = Sys.time()),
          file.path(out_dir, "run_meta.rds"))
}

cat(sprintf("\n\nTotal: %.1f min. Output in %s\n",
            as.numeric(difftime(Sys.time(), t_all, units = "mins")), out_root))

# ── 6. Cross-block summary ───────────────────────────────────────────####
allm <- bind_rows(lapply(BLOCKS, function(b)
  readRDS(file.path(out_root, b, "metrics_k3.rds"))))
saveRDS(allm, file.path(out_root, "metrics_k3_all.rds"))

cat("\n──── LASSO, kappa_emp at k = 3 ────\n")
print(allm %>% filter(classifier == "LASSO") %>%
        select(run_type, method, block, kappa_emp) %>%
        pivot_wider(names_from = block, values_from = kappa_emp) %>%
        arrange(run_type, method), n = 40, digits = 3)

# HDBSCAN vs the matched runs: the only fair comparison, same respondents.
cat("\n──── HDBSCAN vs the same four methods on the SAME respondents ────\n")
print(allm %>% filter(classifier == "LASSO", run_type %in% c("HDB", "MATCHED")) %>%
        select(method, block, n, kappa_emp) %>%
        pivot_wider(names_from = block, values_from = kappa_emp) %>%
        arrange(method), n = 40, digits = 3)

cat("\n──── Selection effect: MATCHED minus FULL, same partition ────\n")
print(allm %>% filter(classifier == "LASSO", run_type %in% c("FULL", "MATCHED")) %>%
        select(run_type, method, block, kappa_emp) %>%
        pivot_wider(names_from = run_type, values_from = kappa_emp) %>%
        mutate(selection_effect = MATCHED - FULL) %>%
        arrange(block, method), n = 40, digits = 3)

