# =============================================================================
# scripts/05_recovery/04_recovery_all_k.R
# Purpose : Out-of-sample recovery of cluster membership from one predictor
#           block, for the four configurations (KM-PC2, KM-PC5, KM-UMAP, PAM)
#           and k = 2..6. Parallel over (method x k) tasks; deterministic
#           given SEED. Normally run once per block by 05_run_all_blocks.R.
# Thesis  : Appendix A.6; Tables A.11 (pol), A.12 (dem), A.13 (pol_dem), A.14
# Input   : results/objects/clusters/all_clusters_aligned.rds
#           data/processed/imputed_labeled.rds
#           (geo blocks) results/objects/Geo/geo_usa_labeled.rds, data/processed/geo_raw.rds
# Output  : results/objects/recovery_all_k/class_<block>/{metrics_long,
#           metrics_wide,results_list,rf_importance,rf_importance_by_variable,
#           lasso_coefs,logit_coefs,rf_top20,lasso_top20,lasso_null_models,
#           term_map,run_meta}.rds and single_cluster/*.rds
# Models  : multinomial logit (nnet), multinomial LASSO (cv.glmnet, inner
#           5-fold, lambda.1se, type.measure = "class"), random forest
#           (ranger, 300 trees, permutation importance); 5-fold stratified CV;
#           SEED = 42. Metrics: accuracy, NIR (p_maj), kappa_emp, kappa_unif.
# Note    : runs in parallel with future (multisession, future.seed = TRUE):
#           documented exception to the "sequential lapply" rule, kept so that
#           the thesis numbers are reproduced exactly (see CLAUDE.md).
#           Factors enter as unordered (ordered = FALSE); model-matrix columns
#           are mapped to variables with the "assign" attribute.
#           05_run_all_blocks.R refuses to start if more than one copy of this
#           script exists in scripts/.
# =============================================================================

BENCHMARK_VERSION <- "2026-08-27"   # 05_run_all_blocks.R checks this

source(here::here("R", "00_setup.R"))

suppressMessages({
  library(dplyr);  library(tidyr)
  library(nnet);   library(glmnet)
  library(future); library(future.apply); library(progressr)
})

# ── 0. Parameters ────────────────────────────────────────────────────####

# 05_run_all_blocks.R sets BATCH_MODE <- TRUE and defines BLOCK / USE_GEO
# BEFORE source()-ing this script. Those values must not be overwritten,
# or every batch iteration silently reruns the same block.
if (!exists("BATCH_MODE", inherits = FALSE) || !isTRUE(BATCH_MODE)) {
  BLOCK   <- "pol_dem"    # "pol" | "dem" | "pol_dem" | "geo"
  USE_GEO <- TRUE         # add geographic block to pol/dem/pol_dem
}
stopifnot(BLOCK %in% c("pol", "dem", "pol_dem", "geo"),
          is.logical(USE_GEO), length(USE_GEO) == 1L)
cat("\n>>> Requested: BLOCK =", BLOCK, "| USE_GEO =", USE_GEO,
    "| script version", BENCHMARK_VERSION, "\n")

GEO_VARS  <- c("geo_region", "geo_urban", "geo_town")
MIN_GEO_N <- 30L          # geo levels rarer than this are lumped into "Other"

K_VALUES  <- 2:6
METHODS   <- c("KM-PC2", "KM-PC5", "KM-UMAP", "PAM")
NFOLDS    <- 5L
GLMNET_NF <- 5L
LAMBDA    <- "lambda.1se"
TYPE_MEAS <- "class"      # "class" | "deviance"
N_TREES   <- 300L
TOP_N     <- 20L
SEED      <- 42L

RF_ENGINE <- "ranger"
RF_IMPORTANCE <- "permutation"   # "permutation" (slow, ~80 min/block)
# "impurity"    (fast, ~10 min/block)
REF_MODE  <- "first"
N_WORKERS <- max(1L, min(6L, parallel::detectCores(logical = FALSE) - 1L))

AGE_COL    <- "Q262"
INCOME_COL <- "Q288"

if (RF_ENGINE == "ranger") library(ranger) else library(randomForest)

# Sequential source() calls share the global environment: drop objects from a
# previous block so a run that fails midway cannot leave stale results behind.
rm(list = intersect(ls(), c(
  "data_pd", "data_geo", "geo_assoc", "X_full", "term_map", "results_list",
  "metrics_long", "metrics_wide", "confusion_list", "per_class_long",
  "bal_acc", "imp_rf_all", "coef_lasso_all", "coef_logit_all", "rf_top20",
  "lasso_top20", "rf_imp_by_var", "y_list", "tasks", "errs", "BLOCK_KEY")))

# ── 1. Output directories ────────────────────────────────────────────####
DIR_FOR_BLOCK <- c(
  pol         = "class_only_pol",
  dem         = "class_only_dem",
  pol_dem     = "class_pol_dem",
  geo         = "class_only_geo",
  pol_geo     = "class_pol_geo",
  dem_geo     = "class_dem_geo",
  pol_dem_geo = "class_pol_dem_geo"
)

BLOCK_KEY <- if (BLOCK == "geo") "geo" else
  if (isTRUE(USE_GEO)) paste0(BLOCK, "_geo") else BLOCK
stopifnot(BLOCK_KEY %in% names(DIR_FOR_BLOCK))

out_dir_obj <- file.path(DIR_OBJ, "recovery_all_k")
out_dir     <- file.path(out_dir_obj, DIR_FOR_BLOCK[[BLOCK_KEY]])
out_dir_cl  <- file.path(out_dir, "single_cluster")
out_dir_tex <- file.path(out_dir, "tex")
for (d in c(out_dir, out_dir_cl, out_dir_tex))
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
cat("Block:", BLOCK_KEY, "| output:", out_dir, "\n")

# ── 2. Data ──────────────────────────────────────────────────────────####
all_clusters_aligned <- readRDS(PATH_CLUSTERS)
data_all <- readRDS(file.path(DIR_DATA, "imputed_labeled.rds"))

N <- nrow(all_clusters_aligned)
stopifnot(nrow(data_all) == N)
cat("N =", N, "\n")

use_geo <- (BLOCK == "geo") || isTRUE(USE_GEO)

base_cols <- switch(BLOCK,
                    "pol"     = Political,
                    "dem"     = Demographics,
                    "pol_dem" = c(Political, Demographics),
                    "geo"     = character(0),
                    stop("unknown block: ", BLOCK))

pd_cols <- intersect(base_cols, colnames(data_all))
lost <- setdiff(base_cols, colnames(data_all))
if (length(lost))
  warning("Block variables not found in data_all (", length(lost), "): ",
          paste(lost, collapse = ", "))
cat("Variables from data_all:", length(pd_cols), "of", length(base_cols), "\n")

data_pd <- if (length(pd_cols)) {
  dplyr::select(data_all, dplyr::all_of(pd_cols))
} else {
  data_all[, 0, drop = FALSE]
}

# ── 2bis. Geographic block ───────────────────────────────────────────####
na_level <- function(x, lab = "Unknown") {
  x <- as.character(x); x[is.na(x)] <- lab; factor(x)
}
lump <- function(x, min_n, other = "Other") {
  x <- as.character(x); tab <- table(x)
  x[x %in% names(tab)[tab < min_n]] <- other
  factor(x)
}
cramers_v <- function(a, b) {
  tb <- table(a, b)
  if (any(dim(tb) < 2)) return(NA_real_)
  chi <- suppressWarnings(chisq.test(tb, correct = FALSE)$statistic)
  sqrt(as.numeric(chi) / (sum(tb) * (min(dim(tb)) - 1)))
}

geo_assoc <- NULL
if (use_geo) {
  geo     <- readRDS(file.path(DIR_OBJ, "Geo", "geo_usa_labeled.rds"))
  geo_raw <- readRDS(file.path(DIR_DATA, "geo_raw.rds"))
  idx <- match(as.integer(geo_raw$D_INTERVIEW), geo$D_INTERVIEW)
  stopifnot(!anyNA(idx), !anyDuplicated(idx), nrow(geo_raw) == N)
  
  missing_geo <- setdiff(GEO_VARS, names(geo))
  if (length(missing_geo))
    stop("GEO_VARS not found in geo_usa_labeled.rds: ",
         paste(missing_geo, collapse = ", "))
  
  data_geo <- as.data.frame(lapply(GEO_VARS, function(v)
    lump(na_level(geo[[v]][idx]), MIN_GEO_N)))
  names(data_geo)    <- GEO_VARS
  rownames(data_geo) <- NULL
  stopifnot(nrow(data_geo) == N,
            !any(names(data_geo) %in% names(data_pd)))
  
  cat("\n──────  Geographic block  ──────\n")
  print(sapply(data_geo, nlevels))
  
  if (length(GEO_VARS) > 1) {
    pairs <- utils::combn(GEO_VARS, 2, simplify = FALSE)
    geo_assoc <- bind_rows(lapply(pairs, function(p)
      data.frame(var1 = p[1], var2 = p[2],
                 cramers_v = cramers_v(data_geo[[p[1]]], data_geo[[p[2]]]))))
    cat("\nPairwise association within geo block (Cramer's V):\n")
    print(geo_assoc, row.names = FALSE, digits = 3)
  }
  
  data_pd <- dplyr::bind_cols(data_pd, data_geo)
  rm(geo, geo_raw, data_geo, idx)
}

stopifnot(ncol(data_pd) > 0, nrow(data_pd) == N)

# ── 3. Age and income binning ────────────────────────────────────────####
looks_numeric <- function(x) {
  if (is.numeric(x)) return(TRUE)
  mean(is.na(suppressWarnings(as.numeric(as.character(x))))) < 0.5
}
bin_age <- function(x) cut(suppressWarnings(as.numeric(as.character(x))),
                           breaks = c(0, 15, 25, 35, 45, 60, Inf), right = FALSE,
                           labels = c("0_15","16_25","26_34","35_44","45_59","60plus"))
bin_income <- function(x) cut(suppressWarnings(as.numeric(as.character(x))),
                              breaks = c(0, 3, 7, 10), right = TRUE,
                              labels = c("Low", "Middle", "High"))

if (AGE_COL %in% names(data_pd) && looks_numeric(data_pd[[AGE_COL]])) {
  data_pd[[AGE_COL]] <- bin_age(data_pd[[AGE_COL]]); cat("\nAge binned.\n")
}
if (INCOME_COL %in% names(data_pd) && looks_numeric(data_pd[[INCOME_COL]])) {
  data_pd[[INCOME_COL]] <- bin_income(data_pd[[INCOME_COL]]); cat("Income binned.\n")
}

# ── 4. Factors + reference level ─────────────────────────────────────####
# ordered = FALSE is essential: ordered factors make model.matrix() emit
# polynomial contrasts, silently invalidating LASSO and logit.
data_pd <- data_pd %>%
  mutate(across(everything(), ~ droplevels(factor(., ordered = FALSE))))
stopifnot(!any(vapply(data_pd, is.ordered, logical(1))))

single_lvl <- names(which(sapply(data_pd, nlevels) < 2))
if (length(single_lvl)) {
  message("Dropped single-level variables: ", paste(single_lvl, collapse = ", "))
  data_pd <- dplyr::select(data_pd, -all_of(single_lvl))
}

if (REF_MODE == "modal") {
  data_pd <- data_pd %>%
    mutate(across(everything(),
                  ~ stats::relevel(., ref = names(which.max(table(.))))))
}
cat("Reference level:", REF_MODE, "\n")

# ── 5. Design matrix ─────────────────────────────────────────────────####
stopifnot(!anyNA(data_pd))
mm <- model.matrix(~ ., data = data_pd)
stopifnot(nrow(mm) == N)

X_full    <- mm[, -1, drop = FALSE]
raw_terms <- colnames(mm)[-1]
assign_id <- attr(mm, "assign")[-1]
var_of    <- names(data_pd)[assign_id]
level_of  <- mapply(function(t, v) sub(paste0("^", v), "", t), raw_terms, var_of)
colnames(X_full) <- make.names(colnames(X_full), unique = TRUE)

term_map <- data.frame(term = colnames(X_full), term_raw = raw_terms,
                       variable = var_of, level = as.character(level_of),
                       stringsAsFactors = FALSE) %>%
  left_join(section_lookup, by = "variable") %>%
  mutate(section_label = ifelse(
    is.na(section_label),
    ifelse(variable %in% GEO_VARS, "Geography", "Other"), section_label))

cat(sprintf("Predictors: %d variables -> %d dummies\n",
            ncol(data_pd), ncol(X_full)))
rm(mm, data_all)

# ── 6. Helpers ───────────────────────────────────────────────────────####
make_folds <- function(y, nfolds, seed) {
  set.seed(seed)
  fold <- integer(length(y))
  for (cl in unique(y)) {
    i <- which(y == cl)
    fold[i] <- sample(rep_len(seq_len(nfolds), length(i)))
  }
  fold
}

cv_multinom <- function(X, y, folds, cfg) {
  Xdf <- as.data.frame(X); oof <- character(length(y))
  models <- list(); conv <- logical(0)
  for (f in sort(unique(folds))) {
    tr <- folds != f; te <- folds == f
    fit <- nnet::multinom(.y ~ ., data = data.frame(.y = y[tr], Xdf[tr, , drop = FALSE]),
                          trace = FALSE, maxit = cfg$maxit, MaxNWts = 1e5)
    oof[te] <- as.character(predict(fit, newdata = Xdf[te, , drop = FALSE],
                                    type = "class"))
    models[[as.character(f)]] <- fit
    conv <- c(conv, isTRUE(fit$convergence == 0))
  }
  list(oof_pred = factor(oof, levels = levels(y)), models = models,
       accuracy = mean(oof == as.character(y)), converged = all(conv))
}

cv_lasso_multi <- function(X, y, folds, cfg) {
  oof <- character(length(y)); models <- list()
  for (f in sort(unique(folds))) {
    tr <- folds != f; te <- folds == f
    set.seed(cfg$seed)
    cvf <- glmnet::cv.glmnet(X[tr, , drop = FALSE], y[tr], family = "multinomial",
                             alpha = 1, nfolds = cfg$glmnet_nf,
                             type.measure = cfg$type_measure)
    oof[te] <- as.character(as.vector(
      predict(cvf, newx = X[te, , drop = FALSE], s = cfg$lambda, type = "class")))
    models[[as.character(f)]] <- cvf
  }
  list(oof_pred = factor(oof, levels = levels(y)), models = models,
       accuracy = mean(oof == as.character(y)))
}

fit_rf <- function(Xtr, ytr, cfg) {
  if (cfg$rf_engine == "ranger")
    ranger::ranger(x = Xtr, y = ytr, num.trees = cfg$n_trees,
                   importance = cfg$rf_importance, num.threads = 1L,
                   classification = TRUE)
  else
    randomForest::randomForest(x = Xtr, y = ytr, ntree = cfg$n_trees,
                               importance = TRUE)
}
predict_rf <- function(fit, Xte, cfg) {
  if (cfg$rf_engine == "ranger")
    as.character(predict(fit, data = Xte, num.threads = 1L)$predictions)
  else as.character(predict(fit, newdata = Xte))
}

cv_rf_multi <- function(X, y, folds, cfg) {
  Xdf <- as.data.frame(X); oof <- character(length(y)); models <- list()
  for (f in sort(unique(folds))) {
    tr <- folds != f; te <- folds == f
    set.seed(cfg$seed)
    rf <- fit_rf(Xdf[tr, , drop = FALSE], y[tr], cfg)
    oof[te] <- predict_rf(rf, Xdf[te, , drop = FALSE], cfg)
    models[[as.character(f)]] <- rf
  }
  list(oof_pred = factor(oof, levels = levels(y)), models = models,
       accuracy = mean(oof == as.character(y)))
}

# ── 7. Importance extraction ─────────────────────────────────────────####
rf_importance <- function(res, cfg) {
  do.call(rbind, lapply(res$models, function(rf) {
    if (cfg$rf_engine == "ranger")
      data.frame(term = names(rf$variable.importance),
                 MeanDecreaseAccuracy = as.numeric(rf$variable.importance),
                 MeanDecreaseGini = NA_real_, stringsAsFactors = FALSE)
    else {
      imp <- randomForest::importance(rf)
      data.frame(term = rownames(imp),
                 MeanDecreaseAccuracy = imp[, "MeanDecreaseAccuracy"],
                 MeanDecreaseGini     = imp[, "MeanDecreaseGini"],
                 stringsAsFactors = FALSE)
    }
  })) %>%
    group_by(term) %>%
    summarise(MeanDecreaseAccuracy = mean(MeanDecreaseAccuracy),
              MeanDecreaseGini = mean(MeanDecreaseGini), .groups = "drop")
}

lasso_importance <- function(res, lambda) {
  bind_rows(lapply(res$models, function(cvf) {
    cl <- coef(cvf, s = lambda)
    do.call(rbind, lapply(names(cl), function(k) {
      m <- as.matrix(cl[[k]])
      data.frame(class = k, term = rownames(m), coef = as.numeric(m[, 1]),
                 stringsAsFactors = FALSE)
    }))
  })) %>%
    filter(term != "(Intercept)") %>%
    group_by(class, term) %>%
    summarise(mean_coef = mean(coef), nonzero_folds = sum(coef != 0),
              .groups = "drop") %>%
    filter(nonzero_folds > 0)
}

logit_coefs <- function(res, y) {
  base_lvl <- levels(y)[1]
  bind_rows(lapply(res$models, function(fit) {
    co <- coef(fit)
    if (is.null(dim(co)))
      data.frame(class = levels(y)[2], term = names(co),
                 coef = as.numeric(co), stringsAsFactors = FALSE)
    else
      do.call(rbind, lapply(rownames(co), function(k)
        data.frame(class = k, term = colnames(co), coef = as.numeric(co[k, ]),
                   stringsAsFactors = FALSE)))
  })) %>%
    filter(term != "(Intercept)") %>%
    group_by(class, term) %>%
    summarise(mean_coef = mean(coef), .groups = "drop") %>%
    mutate(odds_ratio = exp(mean_coef), baseline_class = base_lvl)
}

# ── 8. Single task ───────────────────────────────────────────────────####
run_task <- function(tag, method, k, X, y, cfg) {
  folds <- make_folds(y, cfg$nfolds, cfg$seed)
  
  r_lo <- tryCatch(cv_multinom(X, y, folds, cfg),    error = function(e) e)
  r_la <- tryCatch(cv_lasso_multi(X, y, folds, cfg), error = function(e) e)
  r_rf <- tryCatch(cv_rf_multi(X, y, folds, cfg),    error = function(e) e)
  bad <- function(r) inherits(r, "error")
  
  list(tag = tag, method = method, k = k, truth = y,
       acc_logit = if (bad(r_lo)) NA_real_ else r_lo$accuracy,
       acc_lasso = if (bad(r_la)) NA_real_ else r_la$accuracy,
       acc_rf    = if (bad(r_rf)) NA_real_ else r_rf$accuracy,
       oof_logit = if (bad(r_lo)) NULL else r_lo$oof_pred,
       oof_lasso = if (bad(r_la)) NULL else r_la$oof_pred,
       oof_rf    = if (bad(r_rf)) NULL else r_rf$oof_pred,
       logit_converged = if (bad(r_lo)) NA else r_lo$converged,
       imp_rf     = if (bad(r_rf)) NULL else rf_importance(r_rf, cfg),
       coef_lasso = if (bad(r_la)) NULL else lasso_importance(r_la, cfg$lambda),
       coef_logit = if (bad(r_lo)) NULL else logit_coefs(r_lo, y),
       errors = list(logit = if (bad(r_lo)) conditionMessage(r_lo) else NULL,
                     lasso = if (bad(r_la)) conditionMessage(r_la) else NULL,
                     rf    = if (bad(r_rf)) conditionMessage(r_rf) else NULL))
}

# ── 9. Task grid ─────────────────────────────────────────────────────####
col_for <- function(method, k)
  switch(method,
         "KM-PC2"  = paste0("km_pc2_k", k),
         "KM-PC5"  = paste0("km_pc5_k", k),
         "KM-UMAP" = paste0("km_umap_", k),
         "PAM"     = paste0("pam_",     k),
         stop("unknown method: ", method))

tasks <- expand.grid(method = METHODS, k = K_VALUES, stringsAsFactors = FALSE)
tasks <- tasks[mapply(function(m, k) col_for(m, k) %in% names(all_clusters_aligned),
                      tasks$method, tasks$k), ]
tasks <- tasks[order(-tasks$k, tasks$method), ]
tasks$tag <- sprintf("%s_k%d", tasks$method, tasks$k)
rownames(tasks) <- NULL

y_list <- lapply(seq_len(nrow(tasks)), function(i)
  droplevels(factor(all_clusters_aligned[[col_for(tasks$method[i], tasks$k[i])]])))
names(y_list) <- tasks$tag
stopifnot(all(vapply(y_list, function(y) !anyNA(y) && length(y) == N, logical(1))),
          all(mapply(function(y, k) nlevels(y) == k, y_list, tasks$k)))

CFG <- list(nfolds = NFOLDS, seed = SEED, lambda = LAMBDA,
            glmnet_nf = GLMNET_NF, type_measure = TYPE_MEAS,
            n_trees = N_TREES, rf_engine = RF_ENGINE,
            rf_importance = RF_IMPORTANCE, maxit = 500L)

cat(sprintf("\nTasks: %d | workers: %d | RF importance: %s\n",
            nrow(tasks), N_WORKERS, RF_IMPORTANCE))

# ── 10. Parallel execution ───────────────────────────────────────────####
options(future.globals.maxSize = 2 * 1024^3)
plan(multisession, workers = N_WORKERS)
handlers("txtprogressbar")

t0 <- Sys.time()
results_list <- with_progress({
  prog <- progressor(steps = nrow(tasks))
  future_lapply(
    seq_len(nrow(tasks)),
    function(i) {
      r <- run_task(tasks$tag[i], tasks$method[i], tasks$k[i],
                    X_full, y_list[[i]], CFG)
      prog(sprintf("%s", tasks$tag[i]))
      r
    },
    future.seed = TRUE,
    future.globals = list(
      tasks = tasks, X_full = X_full, y_list = y_list, CFG = CFG, prog = prog,
      run_task = run_task, make_folds = make_folds,
      cv_multinom = cv_multinom, cv_lasso_multi = cv_lasso_multi,
      cv_rf_multi = cv_rf_multi, fit_rf = fit_rf, predict_rf = predict_rf,
      rf_importance = rf_importance, lasso_importance = lasso_importance,
      logit_coefs = logit_coefs),
    future.packages = c("dplyr", "nnet", "glmnet",
                        if (RF_ENGINE == "ranger") "ranger" else "randomForest"))
})
plan(sequential)
names(results_list) <- tasks$tag
cat(sprintf("\nElapsed: %.1f min\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# ── 11. Diagnostics ──────────────────────────────────────────────────####
errs <- unlist(lapply(results_list, function(r) unlist(r$errors)))
if (length(errs)) { warning("Failed fits:"); print(errs) }

nonconv <- names(which(!vapply(results_list,
                               function(r) isTRUE(r$logit_converged), logical(1))))
if (length(nonconv))
  warning("multinom hit maxit (not converged): ", paste(nonconv, collapse = ", "))

# ── 12. Accuracy + chance-corrected metrics ──────────────────────────####
eval_pred <- function(pred, truth) {
  if (is.null(pred))
    return(data.frame(acc = NA_real_, p_e_emp = NA_real_, p_maj = NA_real_,
                      k_eff = NA_integer_, kappa_emp = NA_real_,
                      kappa_unif = NA_real_, kappa_maj = NA_real_))
  lev <- unique(c(as.character(truth), as.character(pred)))
  lev <- lev[order(suppressWarnings(as.numeric(lev)), lev)]
  cm  <- table(pred  = factor(as.character(pred),  levels = lev),
               truth = factor(as.character(truth), levels = lev))
  Nn <- sum(cm); k_eff <- sum(colSums(cm) > 0)
  p_o <- sum(diag(cm)) / Nn
  p_emp <- sum(rowSums(cm) * colSums(cm)) / Nn^2
  p_uni <- 1 / k_eff; p_maj <- max(colSums(cm)) / Nn
  data.frame(acc = p_o, p_e_emp = p_emp, p_maj = p_maj, k_eff = k_eff,
             kappa_emp  = (p_o - p_emp) / (1 - p_emp),
             kappa_unif = (p_o - p_uni) / (1 - p_uni),
             kappa_maj  = (p_o - p_maj) / (1 - p_maj))
}

metrics_long <- bind_rows(lapply(results_list, function(r)
  bind_rows(data.frame(classifier = "Logit", eval_pred(r$oof_logit, r$truth)),
            data.frame(classifier = "LASSO", eval_pred(r$oof_lasso, r$truth)),
            data.frame(classifier = "RF",    eval_pred(r$oof_rf,    r$truth))) %>%
    mutate(block = BLOCK_KEY, method = r$method, k = r$k, .before = 1))) %>%
  arrange(k, method, classifier)

metrics_wide <- metrics_long %>%
  select(k, method, classifier, acc, kappa_emp, kappa_unif) %>%
  pivot_wider(names_from = classifier,
              values_from = c(acc, kappa_emp, kappa_unif),
              names_glue = "{classifier}_{.value}") %>%
  mutate(block = BLOCK_KEY, baseline = 1 / k, .after = method) %>%
  arrange(k, method)

cat("\n────────  ACCURACY / KAPPA (5-fold CV)  ────────\n")
print(as.data.frame(metrics_wide), row.names = FALSE, digits = 3)

# ── 13. Per-cluster metrics + confusion matrices ─────────────────────####
conf_mat <- function(pred, truth) {
  if (is.null(pred)) return(NULL)
  lev <- unique(c(as.character(truth), as.character(pred)))
  lev <- lev[order(suppressWarnings(as.numeric(lev)), lev)]
  table(pred  = factor(as.character(pred),  levels = lev),
        truth = factor(as.character(truth), levels = lev))
}

per_class <- function(cm) {
  tp <- diag(cm); supp <- colSums(cm); pn <- rowSums(cm)
  rec  <- ifelse(supp > 0, tp / supp, NA_real_)
  prec <- ifelse(pn   > 0, tp / pn,   NA_real_)
  data.frame(cluster = names(tp), support = as.integer(supp),
             n_pred = as.integer(pn), recall = as.numeric(rec),
             precision = as.numeric(prec),
             f1 = as.numeric(ifelse(prec + rec > 0,
                                    2 * prec * rec / (prec + rec), 0)),
             stringsAsFactors = FALSE)
}

# unname() is required: results_list is named, so do.call(c, .) would prefix
# every inner name with the outer tag ("KM-PC2_k2.KM-PC2_k2_Logit") and the
# regex below would parse the method as "KM-PC2_k2.KM-PC2".
confusion_list <- do.call(c, unname(lapply(results_list, function(r) {
  preds <- list(Logit = r$oof_logit, LASSO = r$oof_lasso, RF = r$oof_rf)
  out <- lapply(preds, conf_mat, truth = r$truth)
  names(out) <- sprintf("%s_k%d_%s", r$method, r$k, names(preds))
  out[!vapply(out, is.null, logical(1))]
})))
stopifnot(!any(grepl(".", names(confusion_list), fixed = TRUE)))

per_class_long <- bind_rows(lapply(names(confusion_list), function(tag) {
  meta <- regmatches(tag, regexec("^(.*)_k(\\d+)_(\\w+)$", tag))[[1]]
  data.frame(block = BLOCK_KEY, method = meta[2], k = as.integer(meta[3]),
             classifier = meta[4], per_class(confusion_list[[tag]]),
             stringsAsFactors = FALSE)
})) %>% arrange(k, method, classifier, as.integer(cluster))
stopifnot(all(per_class_long$method %in% METHODS))

bal_acc <- per_class_long %>%
  group_by(block, method, k, classifier) %>%
  summarise(balanced_acc  = mean(recall, na.rm = TRUE),
            min_recall    = min(recall,  na.rm = TRUE),
            worst_cluster = cluster[which.min(recall)], .groups = "drop") %>%
  arrange(k, method, classifier)

cat("\n────────  Balanced accuracy / worst cluster  ────────\n")
print(as.data.frame(bal_acc), row.names = FALSE, digits = 3)

# ── 14. Aggregated importance ────────────────────────────────────────####
# At lambda.1se LASSO can shrink EVERY coefficient to zero when the predictor
# block is weak (this happens with geo alone): lasso_importance() then returns
# a 0-row frame, and cbind() on 0 rows errors out. Guard on nrow, not on NULL.
tag_rows <- function(df, method, k) {
  if (is.null(df) || nrow(df) == 0L) return(NULL)
  cbind(method = method, k = k, df, stringsAsFactors = FALSE)
}
empty_like <- function(...) data.frame(..., stringsAsFactors = FALSE)

imp_rf_all <- bind_rows(lapply(results_list, function(r)
  tag_rows(r$imp_rf, r$method, r$k)))
if (nrow(imp_rf_all) == 0L) {
  warning("No RF importance returned for block ", BLOCK_KEY)
  imp_rf_all <- empty_like(method = character(), k = integer(), term = character(),
                           MeanDecreaseAccuracy = numeric(),
                           MeanDecreaseGini = numeric())
}
imp_rf_all <- left_join(imp_rf_all, term_map, by = "term")

coef_lasso_all <- bind_rows(lapply(results_list, function(r)
  tag_rows(r$coef_lasso, r$method, r$k)))
if (nrow(coef_lasso_all) == 0L) {
  warning("LASSO selected nothing anywhere in block ", BLOCK_KEY,
          ": all coefficients shrunk to zero at ", LAMBDA,
          ". This is a substantive result, not a failure — report it.")
  coef_lasso_all <- empty_like(method = character(), k = integer(),
                               class = character(), term = character(),
                               mean_coef = numeric(), nonzero_folds = integer())
}

coef_logit_all <- bind_rows(lapply(results_list, function(r)
  tag_rows(r$coef_logit, r$method, r$k)))
if (nrow(coef_logit_all) == 0L)
  coef_logit_all <- empty_like(method = character(), k = integer(),
                               class = character(), term = character(),
                               mean_coef = numeric(), odds_ratio = numeric(),
                               baseline_class = character())
coef_logit_all <- left_join(coef_logit_all, term_map, by = "term")

lasso_null_models <- bind_rows(lapply(results_list, function(r)
  data.frame(method = r$method, k = r$k,
             lasso_empty = is.null(r$coef_lasso) || nrow(r$coef_lasso) == 0L,
             stringsAsFactors = FALSE))) %>% arrange(k, method)
if (any(lasso_null_models$lasso_empty)) {
  cat("\n────────  LASSO selected no variables in:  ────────\n")
  print(dplyr::filter(lasso_null_models, lasso_empty), row.names = FALSE)
}

rf_top20 <- imp_rf_all %>% group_by(method, k) %>%
  slice_max(MeanDecreaseAccuracy, n = TOP_N, with_ties = FALSE) %>%
  mutate(rank = row_number()) %>% ungroup() %>% arrange(method, k, rank)

lasso_top20 <- coef_lasso_all %>%
  group_by(method, k, term) %>%
  summarise(max_abs = max(abs(mean_coef)), .groups = "drop") %>%
  left_join(term_map, by = "term") %>%
  group_by(method, k) %>%
  slice_max(max_abs, n = TOP_N, with_ties = FALSE) %>%
  mutate(rank = row_number()) %>% ungroup() %>% arrange(method, k, rank)

rf_imp_by_var <- imp_rf_all %>%
  group_by(method, k, variable, section_label) %>%
  summarise(imp_total = sum(MeanDecreaseAccuracy, na.rm = TRUE),
            n_dummies = n(), .groups = "drop") %>%
  arrange(method, k, desc(imp_total))

# ── 15. Save ─────────────────────────────────────────────────────────####
saveRDS(results_list,      file.path(out_dir, "results_list.rds"))
saveRDS(metrics_long,      file.path(out_dir, "metrics_long.rds"))
saveRDS(metrics_wide,      file.path(out_dir, "metrics_wide.rds"))
saveRDS(imp_rf_all,        file.path(out_dir, "rf_importance.rds"))
saveRDS(rf_imp_by_var,     file.path(out_dir, "rf_importance_by_variable.rds"))
saveRDS(coef_lasso_all,    file.path(out_dir, "lasso_coefs.rds"))
saveRDS(coef_logit_all,    file.path(out_dir, "logit_coefs.rds"))
saveRDS(rf_top20,          file.path(out_dir, "rf_top20.rds"))
saveRDS(lasso_top20,       file.path(out_dir, "lasso_top20.rds"))
saveRDS(lasso_null_models, file.path(out_dir, "lasso_null_models.rds"))
saveRDS(term_map,          file.path(out_dir, "term_map.rds"))
saveRDS(per_class_long,    file.path(out_dir_cl, "per_class_metrics.rds"))
saveRDS(confusion_list,    file.path(out_dir_cl, "confusion_matrices.rds"))
saveRDS(bal_acc,           file.path(out_dir_cl, "balanced_accuracy.rds"))
if (!is.null(geo_assoc))
  saveRDS(geo_assoc,       file.path(out_dir, "geo_association.rds"))

saveRDS(list(block = BLOCK_KEY, base_block = BLOCK, use_geo = use_geo,
             geo_vars = if (use_geo) GEO_VARS else NA,
             min_geo_n = MIN_GEO_N, n_obs = N, n_vars = ncol(data_pd),
             n_dummies = ncol(X_full), ref_mode = REF_MODE,
             rf_engine = RF_ENGINE, rf_importance = RF_IMPORTANCE,
             seed = SEED, nfolds = NFOLDS, n_trees = N_TREES,
             lambda = LAMBDA, type_measure = TYPE_MEAS,
             script_version = BENCHMARK_VERSION, run_time = Sys.time()),
        file.path(out_dir, "run_meta.rds"))

cat("\nDone. Output in:", out_dir, "\n")