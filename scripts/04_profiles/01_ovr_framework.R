# =============================================================================
# scripts/04_profiles/01_ovr_framework.R
# Purpose : Framework for the one-vs-rest logits that characterise the
#           clusters (sourced by 02_ovr_lasso.R, 03_method_overlap.R and
#           04_latex_tables.R). Predictors keep their measurement level: the
#           ordinal scales are NOT turned into dummies.
# Thesis  : Sections 3.6, 4.2; Table A.7
# Input   : results/objects/clusters/all_clusters_aligned.rds
#           data/processed/imputed_labeled.rds (value block only)
# Output  : objects in the session (dat_model, X_full, col_var, LABELS,
#           ORD_TABLE, REF_TABLE, helpers); OUT_LOG = results/objects/logit_ovr_ord
#
# -- CHOICES -------------------------------------------------------------------
# [0] Each predictor enters the model according to its type:
#       ORDINAL : ONE numeric column = position of the level on the scale
#                 (1..L), centred on the neutral category if
#                 CENTER_ORDINAL = TRUE. One coefficient = effect of ONE STEP
#                 on the scale (linearity in the ordinal score is assumed).
#       NOMINAL : dummies (L-1 columns) with a neutral reference category.
#       BINARY  : one 0/1 column (special case of nominal: the two
#                 treatments coincide).
#
# [1] The ordinal / nominal classification is EXPLICIT and checkable:
#       - factors already `ordered` -> ordinal;
#       - variables in NOMINAL_DEFAULT (read from the WVS-7 questionnaire)
#         -> nominal;
#       - 2-level factors -> binary (1 column);
#       - everything else -> ordinal, because WVS scales are ordered in the
#         order of their codes.
#     Manual corrections with FORCE_NOMINAL / FORCE_ORDINAL.
#     No ordered factor reaches glm()/glmnet(): ordinal items become numeric
#     scores and nominal items unordered factors (see CLAUDE.md).
#
# [2] Ordinal items use the POSITION of the level, not the WVS code: some
#     scales are reversed (e.g. Q121: 5 = "very good"; Q1: 1 = "very
#     important"). The direction of each coefficient is read in ORD_TABLE,
#     which prints the first/last level of every scale.
#     WARNING: this assumes the factor levels are in questionnaire order and
#     NOT in alphabetical order.
#
# [3] Practical consequence: the design matrix goes from ~700 dummies to
#     just over 200 columns. The LASSO selects QUESTIONS, not single
#     categories, and the unpenalised logits are much more stable (less
#     separation).
#
# [4] One k at a time: K_CHOSEN. All methods, all clusters of that k.
#
# [5] Model: glm(family = binomial), maximum likelihood, no penalty.
#     Separation / non-convergence are caught with tryCatch and flagged.
#
# Columns of the model matrix are mapped back to variables with the "assign"
# attribute (never with regular expressions on the column names).
# =============================================================================

source(here::here("R", "00_setup.R"))   # section_lookup, section_colors, Political, Demographics

# -- Shared parameters ---------------------------------------------------------
K_CHOSEN       <- 3L                                    # <<<<<< the only k analysed (thesis: k = 3)
METHODS        <- c("KM-PC2", "KM-PC5", "KM-UMAP", "PAM")
SEED           <- 42
N_TOP          <- 20
CENTER_ORDINAL <- TRUE   # ordinal score centred on the neutral category
                         # (changes only the intercept, not the coefficients:
                         #  the intercept becomes the "neutral on everything" case)

# -- Paths ---------------------------------------------------------------------
OUT_OBJ <- DIR_OBJ
OUT_LOG <- ensure_dir(OUT_OBJ, "logit_ovr_ord")

# -- Guard against stale objects -----------------------------------------------
rm(list = intersect(ls(),
                    c("predictors_raw", "dat_model", "X_full", "col_var",
                      "LABELS", "VAR_TYPE", "ORD_TABLE", "REF_TABLE")))

# -- Base data -----------------------------------------------------------------
all_clusters_aligned <- readRDS(PATH_CLUSTERS)

predictors_raw <- readRDS(file.path(DIR_DATA, "imputed_labeled.rds")) |>
  dplyr::select(-dplyr::all_of(c(Political, Demographics))) |>
  dplyr::mutate(dplyr::across(dplyr::everything(), as.factor))
# NB: as.factor() leaves factors that are already `ordered` untouched.

n <- nrow(predictors_raw)
stopifnot(nrow(all_clusters_aligned) == n, ncol(predictors_raw) > 0)

# =============================================================================
#  [1] Ordinal / nominal classification
# =============================================================================
# "True" nominal items of the WVS-7 value block (from the questionnaire):
#   Q42       three attitudes towards society (revolution / reform / defence)
#   Q56       standard of living vs parents (better / worse / same -> not monotone)
#   Q91-Q93   knowledge questions (answer A / B / C)
#   Q111      environment vs economic growth (+ "other answer")
#   Q119      women less corrupt (contains "hard to say" = 0)
#   Q149      freedom vs equality ;  Q150  freedom vs security
#   Q152-Q157 postmaterialist index (4 unordered goals)
#   Q174-Q175 meaning of religion (two alternatives)
#   Q221-Q222 local / national vote (contains "not allowed to vote")
#   Q254      national pride (contains "I am not [nationality]")
NOMINAL_DEFAULT <- c("Q42", "Q56", "Q91", "Q92", "Q93", "Q111", "Q119",
                     "Q149", "Q150", "Q152", "Q153", "Q154", "Q155",
                     "Q156", "Q157", "Q174", "Q175", "Q221", "Q222", "Q254")

# manual overrides, e.g. FORCE_NOMINAL <- c("Q173"); FORCE_ORDINAL <- c("Q42")
FORCE_NOMINAL <- character(0)
FORCE_ORDINAL <- character(0)

# pattern for the "neutral / middle" category
NEUTRAL_PAT <- "neither|don'?t mind|hard to say|about the same|might do|inactive|^some$|^fair$|neutral"
# levels that BREAK the ordering if present in a scale treated as ordinal
SUSPECT_PAT <- "hard to say|not applicable|not allowed|don'?t know|no answer|^other|i am not"

classify_var <- function(x, nm) {
  if (nlevels(x) < 2)                                  return("drop")
  if (nm %in% FORCE_ORDINAL)                           return("ordinal")
  if (nm %in% c(FORCE_NOMINAL, NOMINAL_DEFAULT))       return("nominal")
  if (is.ordered(x))                                   return("ordinal")
  if (nlevels(x) == 2)                                 return("binary")
  "ordinal"
}

VAR_TYPE <- vapply(names(predictors_raw),
                   function(nm) classify_var(predictors_raw[[nm]], nm),
                   character(1))

# -- Neutral category: index (ordinal) and level (nominal) ---------------------
neutral_index <- function(x) {
  lv <- levels(x); low <- tolower(trimws(lv))
  hit <- which(grepl(NEUTRAL_PAT, low))
  if (length(hit) >= 1) return(hit[1])
  ceiling(length(lv) / 2)          # odd scales: exact midpoint
}                                  # even scales: lower central level

pick_reference <- function(x) {
  lv <- levels(x); nl <- length(lv); low <- tolower(trimws(lv))
  hit <- which(grepl(NEUTRAL_PAT, low))
  if (length(hit) >= 1) return(lv[hit[1]])
  if (nl == 2) {
    neg <- which(grepl("not mentioned|^no$|^not |don'?t belong|would never|^never", low))
    if (length(neg) >= 1) return(lv[neg[1]])
    return(lv[1])
  }
  lv[ceiling(nl / 2)]
}

# -- Modelling data frame ------------------------------------------------------
dat_model <- predictors_raw
ord_rows  <- list()
ref_rows  <- list()

for (nm in names(dat_model)) {
  x <- droplevels(dat_model[[nm]])
  if (nlevels(x) < 2) { VAR_TYPE[[nm]] <- "drop"; next }   # constant: dropped

  if (VAR_TYPE[[nm]] == "ordinal") {
    idx   <- neutral_index(x)
    score <- as.numeric(as.integer(x))                     # position of the level
    if (CENTER_ORDINAL) score <- score - idx
    dat_model[[nm]] <- score
    ord_rows[[nm]] <- data.frame(
      variable      = nm,
      n_levels      = nlevels(x),
      first_level   = levels(x)[1],
      last_level    = levels(x)[nlevels(x)],
      neutral_level = levels(x)[idx],
      centred       = CENTER_ORDINAL,
      suspect       = any(grepl(SUSPECT_PAT, tolower(levels(x)))),
      stringsAsFactors = FALSE)
  } else {
    ref <- pick_reference(x)
    dat_model[[nm]] <- relevel(factor(x, levels = levels(x), ordered = FALSE),
                               ref = ref)
    ref_rows[[nm]] <- data.frame(
      variable = nm, type = VAR_TYPE[[nm]], n_levels = nlevels(x),
      reference = ref, stringsAsFactors = FALSE)
  }
}

dat_model <- dat_model[, VAR_TYPE[names(dat_model)] != "drop", drop = FALSE]
ORD_TABLE <- if (length(ord_rows)) do.call(rbind, ord_rows) else NULL
REF_TABLE <- if (length(ref_rows)) do.call(rbind, ref_rows) else NULL
rownames(ORD_TABLE) <- NULL; rownames(REF_TABLE) <- NULL

# -- Design matrix: 1 column per ordinal, L-1 dummies per nominal --------------
build_design <- function(df) {
  mm   <- stats::model.matrix(~ ., data = df)
  asg  <- attr(mm, "assign")
  tl   <- attr(stats::terms(~ ., data = df), "term.labels")
  keep <- asg != 0                                   # drop the intercept
  list(mm = mm[, keep, drop = FALSE], col_var = tl[asg[keep]])
}

dm      <- build_design(dat_model)
X_full  <- dm$mm            # colnames = term names (= variable if ordinal)
col_var <- dm$col_var       # original variable of each column
stopifnot(nrow(X_full) == n, length(col_var) == ncol(X_full))

# map term -> original variable / type of term
term_to_var <- function(tr) {
  i <- match(tr, colnames(X_full))
  out <- rep(NA_character_, length(tr))
  out[!is.na(i)] <- col_var[i[!is.na(i)]]
  out
}
term_type <- function(tr) {
  v <- term_to_var(tr)
  ifelse(is.na(v), "intercept", ifelse(tr == v, "ordinal", "dummy"))
}

# -- Report of the parametrisation ---------------------------------------------
cat(sprintf("\nPredictors: %d usable / %d (%d constant dropped)\n",
            ncol(dat_model), ncol(predictors_raw),
            sum(VAR_TYPE == "drop")))
cat(sprintf("  ordinal  (1 column)  : %d\n", sum(VAR_TYPE == "ordinal")))
cat(sprintf("  binary   (1 column)  : %d\n", sum(VAR_TYPE == "binary")))
cat(sprintf("  nominal  (L-1 dummy) : %d\n", sum(VAR_TYPE == "nominal")))
cat(sprintf("  -> total columns of the design matrix: %d\n\n", ncol(X_full)))

cat("Ordinal scales (coefficient direction: from first_level to last_level) — first 12:\n")
print(utils::head(ORD_TABLE, 12), row.names = FALSE)

susp <- ORD_TABLE$variable[ORD_TABLE$suspect]
if (length(susp))
  cat("\nWARNING — ordinal items with a level that breaks the ordering",
      "(consider FORCE_NOMINAL):\n  ", paste(susp, collapse = ", "), "\n", sep = "")

if (!is.null(REF_TABLE)) {
  cat("\nReference categories of the nominal/binary items — first 12:\n")
  print(utils::head(REF_TABLE, 12), row.names = FALSE)
}

# =============================================================================
#  Cluster labels for the chosen k only
# =============================================================================
col_for <- function(method, k) {
  switch(method,
         "KM-PC2"  = paste0("km_pc2_k", k),
         "KM-PC5"  = paste0("km_pc5_k", k),
         "KM-UMAP" = paste0("km_umap_", k),
         "PAM"     = paste0("pam_",     k),
         stop("unknown method: ", method))
}

LABELS <- setNames(lapply(METHODS, function(m) {
  cn <- col_for(m, K_CHOSEN)
  if (!cn %in% names(all_clusters_aligned))
    stop("Column missing in all_clusters_aligned: ", cn,
         "  (method ", m, ", k = ", K_CHOSEN, ")")
  droplevels(factor(all_clusters_aligned[[cn]]))
}), METHODS)

stopifnot(all(vapply(LABELS, length, integer(1)) == n))
cat(sprintf("\nk = %d | clusters per method: %s\n", K_CHOSEN,
            paste(sprintf("%s=%d", METHODS,
                          vapply(LABELS, nlevels, integer(1))), collapse = ", ")))

# -- One-vs-rest logit helpers -------------------------------------------------
ovr_y <- function(labels, cl)
  factor(ifelse(as.character(labels) == as.character(cl), "in", "out"),
         levels = c("out", "in"))

# glm on the chosen VARIABLES, each in its natural form
# (ordinal = numeric; nominal = factor with neutral reference)
fit_logit <- function(y, vars, method, k, cluster, selector) {
  vars <- intersect(vars, names(dat_model))
  if (length(vars) == 0)
    return(list(coefs = NULL, status = "no_vars", n_vars = 0L))
  d  <- cbind(.y = y, dat_model[vars])
  fm <- stats::reformulate(vars, response = ".y")
  m  <- tryCatch(suppressWarnings(
          glm(fm, data = d, family = binomial(),
              control = glm.control(maxit = 100))),
        error = function(e) e)
  if (inherits(m, "error"))
    return(list(coefs = NULL, status = paste0("error: ", conditionMessage(m)),
                n_vars = length(vars)))
  s      <- summary(m)$coefficients
  status <- if (!isTRUE(m$converged)) "not_converged"
            else if (any(abs(stats::coef(m)[-1]) > 15, na.rm = TRUE)) "possible_separation"
            else "ok"
  coefs <- data.frame(
    method = method, k = k, cluster = as.character(cluster), selector = selector,
    term = rownames(s), variable = term_to_var(rownames(s)),
    type = term_type(rownames(s)),
    estimate = s[, 1], std_error = s[, 2], z = s[, 3], p_value = s[, 4],
    odds_ratio = exp(s[, 1]), status = status, row.names = NULL)
  list(coefs = coefs, status = status, n_vars = length(vars))
}

# -- Grouped data (useful only if few variables remain) ------------------------
# NB: with ordinal scores almost all patterns are singletons; grouping is
#     informative only for very short selections.
.pattern_key <- function(X) do.call(paste, c(lapply(X, as.character), sep = "\r"))

grouping_feasibility <- function(vars) {
  vars <- intersect(vars, names(dat_model))
  if (length(vars) == 0) return(NULL)
  key <- .pattern_key(dat_model[vars]); tab <- table(key)
  data.frame(n_vars = length(vars), n_obs = length(key),
             n_patterns = length(tab),
             obs_per_pat = round(length(key) / length(tab), 2),
             pct_singleton = round(mean(tab == 1) * 100, 1),
             n_pat_ge5 = sum(tab >= 5),
             cov_ge5_pct = round(sum(tab[tab >= 5]) / length(key) * 100, 1))
}

fit_logit_grouped <- function(y, vars) {
  vars <- intersect(vars, names(dat_model)); if (length(vars) == 0) return(NULL)
  yb <- as.integer(y == "in"); X <- dat_model[vars]
  key <- .pattern_key(X)
  succ <- tapply(yb, key, sum); trials <- tapply(yb, key, length)
  rep_idx <- match(names(succ), key)
  dg <- X[rep_idx, , drop = FALSE]
  dg$.succ <- as.integer(succ); dg$.fail <- as.integer(trials - succ)
  fm <- stats::reformulate(vars, response = "cbind(.succ, .fail)")
  suppressWarnings(glm(fm, data = dg, family = binomial(),
                       control = glm.control(maxit = 100)))
}

residual_check <- function(mod, which = c(1, 2, 3, 5), main = "") {
  op <- par(mfrow = c(2, 2), oma = c(0, 0, 2, 0)); on.exit(par(op))
  plot(mod, which = which)
  if (nzchar(main)) mtext(main, outer = TRUE, cex = 1.1)
  dev <- deviance(mod); dfr <- df.residual(mod)
  cat(sprintf("Residual deviance = %.1f on %d df  (GOF p = %.3f; valid only if trials are adequate)\n",
              dev, dfr, pchisq(dev, dfr, lower.tail = FALSE)))
}

cat("\nOrdinal framework loaded. Logit output in:", OUT_LOG, "\n")
