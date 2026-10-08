# =============================================================================
# scripts/04_profiles/02_ovr_lasso.R
# Purpose : LASSO + two one-vs-rest logits for k = K_CHOSEN only (set in
#           01_ovr_framework.R; thesis: k = 3), on all methods and all
#           clusters of that k.
#   (A) Logistic LASSO (cv.glmnet) on the mixed design matrix: one column
#       for each ordinal scale, dummies only for the nominal items.
#   (B) Logit with the TOP_N most important QUESTIONS for the LASSO
#       (ordinal items as scores, nominal items as full factors).
#   (C) Logit with the N_REG strongest REGRESSORS: for an ordinal item the
#       regressor coincides with the question, for a nominal item it is the
#       single category. (B) and (C) therefore differ ONLY when nominal
#       items are present.
# Thesis  : Sections 3.6, 4.2 (one-vs-rest LASSO, unpenalised refit)
# Input   : objects from 01_ovr_framework.R
# Output  : results/objects/logit_ovr_ord/ (suffix k<K_CHOSEN>):
#   logit_lasso_selected_k*.rds          variables selected by the LASSO
#   logit_lasso_coefs_glmnet_k*.rds      non-zero LASSO coefficients (+ OR)
#   logit_lasso_top20_selected_k*.rds    (B) the TOP_N questions
#   logit_lasso_top20_coefs_k*.rds       (B) glm coefficients
#   logit_lasso_top20reg_selected_k*.rds (C) the N_REG regressors
#   logit_lasso_top20reg_coefs_k*.rds    (C) glm coefficients
#   ordinal_scale_table_k*.rds           audit of the ordinal scales
#   neutral_reference_table_k*.rds       audit of the reference categories
#
# CHOICES:
#   - lambda = lambda.1se; nfolds = 10; set.seed(SEED = 42) before each fit
#   - standardize = TRUE  <-- REQUIRED here: the columns are no longer all
#     0/1 (a 10-point scale has sd ~ 2.5, a dummy ~ 0.5). Without
#     standardisation the L1 penalty would mechanically favour long scales.
#     glmnet still returns the coefficients on the original scale.
#   - the coefficients of the ordinal items are the effect of ONE STEP on the
#     scale; or_full_range = effect of the whole range (L-1 steps).
# =============================================================================

source(here::here("scripts", "04_profiles", "01_ovr_framework.R"))

LAMBDA          <- "lambda.1se"
NFOLDS          <- 10
TOP_N           <- 20      # number of QUESTIONS   (B)
N_REG           <- 20      # number of REGRESSORS  (C)
RECOMPUTE_LASSO <- FALSE   # TRUE to rerun cv.glmnet ignoring the saved files

tag_k  <- sprintf("k%d", K_CHOSEN)
f_sel  <- file.path(OUT_LOG, sprintf("logit_lasso_selected_%s.rds",     tag_k))
f_coef <- file.path(OUT_LOG, sprintf("logit_lasso_coefs_glmnet_%s.rds", tag_k))

rm(list = intersect(ls(), c("lasso_selected", "lasso_coefs_raw",
                            "top_lasso_vars", "top_lasso_terms",
                            "lasso_top_coefs", "lasso_topreg_coefs")))

# ════════════════════════════════════════════════════════════════════
#  (A) LASSO
# ════════════════════════════════════════════════════════════════════
if (!RECOMPUTE_LASSO && file.exists(f_sel) && file.exists(f_coef)) {
  cat("Found LASSO output for k =", K_CHOSEN, "-> reusing it.\n")
  lasso_selected  <- readRDS(f_sel)
  lasso_coefs_raw <- readRDS(f_coef)
} else {
  cat("Running cv.glmnet for k =", K_CHOSEN, "...\n")
  sel_list <- list(); lasso_coef_list <- list()
  for (m in METHODS) {
    labels <- LABELS[[m]]
    for (cl in levels(labels)) {
      tag <- sprintf("%s_k%d_c%s", m, K_CHOSEN, cl)
      yb  <- as.integer(ovr_y(labels, cl) == "in")
      set.seed(SEED)
      cv <- tryCatch(cv.glmnet(X_full, yb, family = "binomial", alpha = 1,
                               nfolds = NFOLDS, standardize = TRUE),
                     error = function(e) e)
      if (inherits(cv, "error")) { cat(sprintf("LASSO %s: error\n", tag)); next }

      cf <- as.matrix(coef(cv, s = LAMBDA)); nz <- cf[, 1][cf[, 1] != 0]
      if (length(nz) > 0)
        lasso_coef_list[[tag]] <- data.frame(
          method = m, k = K_CHOSEN, cluster = as.character(cl),
          term = names(nz), variable = term_to_var(names(nz)),
          type = term_type(names(nz)),
          coef_lasso = as.numeric(nz), odds_ratio = exp(as.numeric(nz)),
          lambda = cv[[LAMBDA]], stringsAsFactors = FALSE)

      co <- cf[-1, 1]
      sel_vars <- unique(col_var[co != 0]); sel_vars <- sel_vars[!is.na(sel_vars)]
      sel_list[[tag]] <- data.frame(method = m, k = K_CHOSEN,
                                    cluster = as.character(cl),
                                    variable = sel_vars, stringsAsFactors = FALSE)
      cat(sprintf("LASSO %-22s  vars=%3d\n", tag, length(sel_vars)))
    }
  }
  lasso_selected <- dplyr::bind_rows(sel_list) |>
    dplyr::left_join(section_lookup, by = "variable") |>
    dplyr::mutate(section_label = ifelse(is.na(section_label), "Other", section_label)) |>
    dplyr::arrange(method, cluster, variable)
  lasso_coefs_raw <- dplyr::bind_rows(lasso_coef_list) |>
    dplyr::left_join(section_lookup, by = "variable") |>
    dplyr::mutate(section_label = ifelse(is.na(section_label), "Other", section_label))
}

# odds ratio over the whole range of the scale (ordinal items only)
if (!is.null(ORD_TABLE)) {
  span <- setNames(ORD_TABLE$n_levels - 1L, ORD_TABLE$variable)
  lasso_coefs_raw$or_full_range <- ifelse(
    lasso_coefs_raw$type == "ordinal",
    exp(lasso_coefs_raw$coef_lasso * span[lasso_coefs_raw$variable]),
    NA_real_)
}

# ════════════════════════════════════════════════════════════════════
#  (B) LOGIT with the TOP_N QUESTIONS
# ════════════════════════════════════════════════════════════════════
top_lasso_vars <- lasso_coefs_raw |>
  dplyr::filter(term != "(Intercept)", !is.na(variable)) |>
  dplyr::group_by(method, k, cluster, variable) |>
  dplyr::summarise(max_abs = max(abs(coef_lasso)), .groups = "drop") |>
  dplyr::left_join(section_lookup, by = "variable") |>
  dplyr::mutate(section_label = ifelse(is.na(section_label), "Other", section_label)) |>
  dplyr::group_by(method, k, cluster) |>
  dplyr::slice_max(max_abs, n = TOP_N, with_ties = FALSE) |>
  dplyr::mutate(rank = dplyr::row_number()) |>
  dplyr::arrange(method, cluster, rank) |>
  dplyr::ungroup()

# NB: ranking by max|coef| compares different scales. For a ranking
# on a common scale, order by abs(coef) * sd of the column:
#   sd_col <- apply(X_full, 2, sd)
vars_for <- function(m, cl)
  top_lasso_vars$variable[top_lasso_vars$method == m &
                            top_lasso_vars$cluster == as.character(cl)]

var_coef_list <- list()
for (m in METHODS) {
  labels <- LABELS[[m]]
  for (cl in levels(labels)) {
    vars <- vars_for(m, cl)
    if (length(vars) == 0) next
    fit <- fit_logit(ovr_y(labels, cl), vars, m, K_CHOSEN, cl,
                     selector = sprintf("LASSO-top%dq", TOP_N))
    if (!is.null(fit$coefs))
      var_coef_list[[sprintf("%s_c%s", m, cl)]] <- fit$coefs
    cat(sprintf("%d-QUESTIONS %-18s  vars=%2d  status=%s\n",
                TOP_N, sprintf("%s_c%s", m, cl), length(vars), fit$status))
  }
}
lasso_top_coefs <- dplyr::bind_rows(var_coef_list)

# question-level glm for a single case, for manual inspection
get_glm_var <- function(m, cl) {
  vars <- intersect(vars_for(m, cl), names(dat_model))
  if (length(vars) == 0) stop("no variables for ", m, " c", cl)
  y <- ovr_y(LABELS[[m]], cl)
  glm(stats::reformulate(vars, ".y"),
      data = cbind(.y = y, dat_model[vars]),
      family = binomial(), control = glm.control(maxit = 100))
}

# ════════════════════════════════════════════════════════════════════
#  (C) LOGIT with the N_REG strongest REGRESSORS
# ════════════════════════════════════════════════════════════════════
top_lasso_terms <- lasso_coefs_raw |>
  dplyr::filter(term != "(Intercept)") |>
  dplyr::group_by(method, k, cluster) |>
  dplyr::slice_max(abs(coef_lasso), n = N_REG, with_ties = FALSE) |>
  dplyr::mutate(rank = dplyr::row_number()) |>
  dplyr::arrange(method, cluster, rank) |>
  dplyr::ungroup()

get_glm_reg <- function(m, cl, nreg = N_REG) {
  d <- subset(lasso_coefs_raw, method == m & k == K_CHOSEN &
                cluster == as.character(cl) & term != "(Intercept)")
  terms_sel <- head(d$term[order(-abs(d$coef_lasso))], nreg)
  terms_sel <- terms_sel[terms_sel %in% colnames(X_full)]
  if (length(terms_sel) == 0) stop("no regressors for ", m, " c", cl)
  Xs   <- X_full[, terms_sel, drop = FALSE]
  safe <- make.names(terms_sel, unique = TRUE)
  colnames(Xs) <- safe
  y   <- as.integer(ovr_y(LABELS[[m]], cl) == "in")
  mod <- glm(.y ~ ., data = data.frame(.y = y, Xs, check.names = FALSE),
             family = binomial(), control = glm.control(maxit = 100))
  attr(mod, "term_map") <- setNames(terms_sel, safe)
  mod
}

# tidy table (summary + OR) with the original term names
show_glm <- function(mod, digits = 3) {
  s  <- summary(mod)$coefficients
  tm <- attr(mod, "term_map")
  orig <- rownames(s)
  if (!is.null(tm))
    orig <- ifelse(rownames(s) == "(Intercept)", "(Intercept)", tm[rownames(s)])
  or <- exp(cbind(OR = coef(mod), confint.default(mod)))
  status <- if (!isTRUE(mod$converged)) "not_converged"
            else if (any(abs(coef(mod)[-1]) > 15, na.rm = TRUE)) "possible_separation"
            else "ok"
  data.frame(term = orig, type = term_type(orig),
             estimate = round(s[, 1], digits), std_error = round(s[, 2], digits),
             z = round(s[, 3], digits), p_value = signif(s[, 4], 3),
             OR = round(or[, 1], digits), lwr = round(or[, 2], digits),
             upr = round(or[, 3], digits), status = status, row.names = NULL)
}

reg_coef_list <- list()
for (m in METHODS) {
  for (cl in levels(LABELS[[m]])) {
    mod <- tryCatch(get_glm_reg(m, cl), error = function(e) e)
    if (inherits(mod, "error")) next
    tab <- cbind(method = m, k = K_CHOSEN, cluster = as.character(cl),
                 selector = sprintf("LASSO-top%dreg", N_REG), show_glm(mod))
    reg_coef_list[[sprintf("%s_c%s", m, cl)]] <- tab
    cat(sprintf("%d-REGR    %-18s  n=%2d  status=%s\n",
                N_REG, sprintf("%s_c%s", m, cl), nrow(tab) - 1, tab$status[1]))
  }
}
lasso_topreg_coefs <- dplyr::bind_rows(reg_coef_list)

# -- Print: all methods, all clusters of K_CHOSEN ------------------------------
for (m in METHODS) {
  for (cl in levels(LABELS[[m]])) {
    cat(sprintf("\n=====  %s — k=%d  c=%s  (%d LASSO regressors)  =====\n",
                m, K_CHOSEN, cl, N_REG))
    mod <- tryCatch(get_glm_reg(m, cl), error = function(e) e)
    if (inherits(mod, "error")) { cat("  (cluster not available)\n"); next }
    cat(sprintf("AIC=%.1f | dev.resid=%.1f on %d df | iter=%d\n",
                AIC(mod), deviance(mod), df.residual(mod), mod$iter))
    print(show_glm(mod), row.names = FALSE)
  }
}
# question-level version for a single case:
#   print(show_glm(get_glm_var("KM-PC2", 2)), row.names = FALSE)
#   residual_check(get_glm_var("KM-PC2", 2), main = "KM-PC2 c2")

# -- Save (all at the end) ------------------------------------------------------
saveRDS(lasso_selected,     file.path(OUT_LOG, sprintf("logit_lasso_selected_%s.rds",        tag_k)))
saveRDS(lasso_coefs_raw,    file.path(OUT_LOG, sprintf("logit_lasso_coefs_glmnet_%s.rds",    tag_k)))
saveRDS(top_lasso_vars,     file.path(OUT_LOG, sprintf("logit_lasso_top20_selected_%s.rds",  tag_k)))
saveRDS(lasso_top_coefs,    file.path(OUT_LOG, sprintf("logit_lasso_top20_coefs_%s.rds",     tag_k)))
saveRDS(top_lasso_terms,    file.path(OUT_LOG, sprintf("logit_lasso_top20reg_selected_%s.rds", tag_k)))
saveRDS(lasso_topreg_coefs, file.path(OUT_LOG, sprintf("logit_lasso_top20reg_coefs_%s.rds",  tag_k)))
saveRDS(ORD_TABLE,          file.path(OUT_LOG, sprintf("ordinal_scale_table_%s.rds",         tag_k)))
saveRDS(REF_TABLE,          file.path(OUT_LOG, sprintf("neutral_reference_table_%s.rds",     tag_k)))

cat("\nFatto (k =", K_CHOSEN, "). Output in:", OUT_LOG, "\n")

