# =============================================================================
# scripts/04_profiles/04_latex_tables.R
# Purpose : LaTeX tables of the one-vs-rest logits on the N_REG strongest
#           LASSO regressors, for k = K_CHOSEN only (run after
#           02_ovr_lasso.R): one longtable per METHOD, with the clusters as
#           inner blocks, plus the coding table of the ordinal scales.
# Thesis  : Section 4.2, Appendix A.3.3
# Input   : results/objects/logit_ovr_ord/logit_lasso_top20reg_coefs_k<K>.rds
# Output  : results/tables/profiles/logit_top20reg_<METHOD>_k<K>.tex,
#           logit_top20reg_all_k<K>.tex, ordinal_scales_k<K>.tex
#
# -- CHOICES -------------------------------------------------------------------
# [0] The tables are generated from the SAVED OBJECTS, not from the text
#     printed to the console: the console output wraps the columns, and
#     copying it by hand is the fastest way to a typo in the thesis.
#
# [1] Columns: LASSO rank, variable, type, contrast, estimate, standard
#     error, z, p, odds ratio, 95% CI. "Contrast" is the category against the
#     reference for dummies and "per step" for ordinal items (one step on the
#     scale). The `status` column appears only if some model is not "ok";
#     otherwise convergence is stated in the note.
#
# [2] 95% CI in a single cell [lwr, upr]: both bounds are there, but the
#     table stays readable. SPLIT_CI = TRUE keeps them in two columns.
#
# [3] Rows with p >= 0.05 shaded in grey (\rowcolor): the unpenalised logit
#     is a post-LASSO refit, so the p-values are optimistic and must be read
#     as descriptive. See the note.
#
# [4] Block header for each cluster with n, events, AIC and residual
#     deviance, recomputed with get_glm_reg() (available when 02_ovr_lasso.R
#     has been run in the same session).
#
# Required LaTeX preamble:
#   \usepackage{booktabs, longtable, array, colortbl, xcolor}
#   \usepackage{pdflscape}      % only if LANDSCAPE = TRUE
# =============================================================================

# -- Loading -------------------------------------------------------------------
# Always read the coefficients for the current K_CHOSEN from disk (never a
# stale object left in the session by a run with another k).
if (!exists("K_CHOSEN") || !exists("dat_model"))
  source(here::here("scripts", "04_profiles", "01_ovr_framework.R"))
lasso_topreg_coefs <- readRDS(file.path(
  OUT_LOG, sprintf("logit_lasso_top20reg_coefs_k%d.rds", K_CHOSEN)))

SPLIT_CI    <- FALSE   # TRUE -> lwr and upr in two columns
LANDSCAPE   <- TRUE    # wraps the longtables in \begin{landscape}
GREY_NONSIG <- TRUE    # shades the rows with p >= 0.05
TAB_DIR     <- ensure_dir(DIR_TAB, "profiles")

tag_k <- sprintf("k%d", K_CHOSEN)

# -- Formatting helpers --------------------------------------------------------
esc <- function(x) {
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("([&%$#_{}])", "\\\\\\1", x)
  gsub("~", "\\\\textasciitilde{}", x)
}

fmt_n <- function(x, d = 3) ifelse(is.na(x), "--", sprintf(paste0("%.", d, "f"), x))

fmt_p <- function(p) {
  out <- character(length(p))
  ex  <- ifelse(p > 0, floor(log10(p)), 0)
  sci <- !is.na(p) & p < 1e-3
  out[!sci] <- sprintf("%.3f", p[!sci])
  out[sci]  <- sprintf("$%.1f{\\times}10^{%d}$", (p / 10^ex)[sci], ex[sci])
  out[is.na(p)] <- "--"
  out
}

# -- Data preparation ----------------------------------------------------------
d <- lasso_topreg_coefs
d$cluster  <- as.character(d$cluster)
d$variable <- term_to_var(d$term)
d$contrast <- mapply(function(tr, v, ty) {
  if (ty == "intercept") return("--")
  if (ty == "ordinal")   return("per step")
  substring(tr, nchar(v) + 1)
}, d$term, d$variable, d$type, USE.NAMES = FALSE)
d$variable[d$type == "intercept"] <- "(Intercept)"

# LASSO rank: the order of the rows is already that of |coef_lasso|
d <- do.call(rbind, lapply(split(d, list(d$method, d$cluster), drop = TRUE),
                           function(x) { x$rank <- seq_len(nrow(x)) - 1L; x }))
d$rank[d$type == "intercept"] <- NA_integer_
rownames(d) <- NULL

SHOW_STATUS <- any(d$status != "ok")

# fit statistics for each block
fit_stats <- function(m, cl) {
  if (!exists("get_glm_reg")) return(NULL)
  mod <- tryCatch(get_glm_reg(m, cl), error = function(e) NULL)
  if (is.null(mod)) return(NULL)
  y <- mod$model$.y
  sprintf("$n = %d$ ($n_{\\text{in}} = %d$), AIC $= %.1f$, residual deviance $= %.1f$ on %d df",
          length(y), sum(y == 1), AIC(mod), deviance(mod), df.residual(mod))
}

# -- One longtable per method --------------------------------------------------
ncols_ci  <- if (SPLIT_CI) 2L else 1L
head_ci   <- if (SPLIT_CI) "Lower & Upper" else "95\\% CI"
align_ci  <- if (SPLIT_CI) "rr" else "c"
head_stat <- if (SHOW_STATUS) " & Status" else ""
alig_stat <- if (SHOW_STATUS) "l" else ""
NCOL      <- 8L + ncols_ci + as.integer(SHOW_STATUS)

header <- paste0("\\# & Variable & Type & Contrast & Estimate & SE & $z$ & $p$ & OR & ",
                 head_ci, head_stat, " \\\\")

make_table <- function(m) {
  dm  <- d[d$method == m, ]
  cls <- sort(unique(dm$cluster))
  lab <- sprintf("tab:logit-%s-%s", gsub("[^A-Za-z0-9]", "", m), tag_k)
  cap <- sprintf(paste("One-vs-rest logistic regression on the %d strongest LASSO",
                       "regressors, %s clustering, $k = %d$."),
                 max(dm$rank, na.rm = TRUE), esc(m), K_CHOSEN)

  L <- c(
    if (LANDSCAPE) "\\begin{landscape}",
    "\\begingroup\\footnotesize",
    sprintf("\\begin{longtable}{r l l l r r r r r %s%s}", align_ci, alig_stat),
    sprintf("\\caption{%s}\\label{%s}\\\\", cap, lab),
    "\\toprule", header, "\\midrule", "\\endfirsthead",
    sprintf("\\multicolumn{%d}{l}{\\emph{\\tablename\\ \\thetable{} -- continued}}\\\\", NCOL),
    "\\toprule", header, "\\midrule", "\\endhead",
    sprintf("\\midrule \\multicolumn{%d}{r}{\\emph{continued on next page}}\\\\", NCOL),
    "\\endfoot", "\\bottomrule", "\\endlastfoot")

  for (cl in cls) {
    dc <- dm[dm$cluster == cl, ]
    fs <- fit_stats(m, cl)
    L <- c(L,
      sprintf("\\multicolumn{%d}{l}{\\textbf{Cluster %s}%s}\\\\[2pt]", NCOL, esc(cl),
              if (is.null(fs)) "" else paste0(" \\quad \\small ", fs)))
    for (i in seq_len(nrow(dc))) {
      r  <- dc[i, ]
      ci <- if (SPLIT_CI) paste(fmt_n(r$lwr), "&", fmt_n(r$upr))
            else sprintf("[%s,\\,%s]", fmt_n(r$lwr), fmt_n(r$upr))
      row <- paste(
        ifelse(is.na(r$rank), "", r$rank),
        esc(r$variable), r$type, esc(r$contrast),
        fmt_n(r$estimate), fmt_n(r$std_error), fmt_n(r$z, 2),
        fmt_p(r$p_value), fmt_n(r$OR), ci,
        if (SHOW_STATUS) esc(r$status) else NULL,
        sep = " & ")
      grey <- GREY_NONSIG && !is.na(r$p_value) && r$p_value >= 0.05
      L <- c(L, paste0(if (grey) "\\rowcolor{gray!12} ", row, " \\\\"))
    }
    if (cl != cls[length(cls)]) L <- c(L, "\\addlinespace[4pt]")
  }

  L <- c(L,
    sprintf(paste0("\\multicolumn{%d}{p{0.95\\linewidth}}{\\footnotesize \\emph{Notes:} ",
                   "regressors selected by one-vs-rest LASSO ($\\lambda_{1\\text{se}}$) and ",
                   "refitted without penalty; $p$-values are therefore descriptive, not ",
                   "inferential. Ordinal variables enter as a single score, so the estimate ",
                   "is the effect of one step along the scale; the direction of each scale ",
                   "is given in Table~\\ref{tab:ordinal-scales}. Dummy contrasts are read ",
                   "against the neutral reference category. %s%s}\\\\"),
            NCOL,
            if (GREY_NONSIG) "Shaded rows have $p \\geq 0.05$. " else "",
            if (SHOW_STATUS) "" else "All models converged."),
    "\\end{longtable}", "\\endgroup",
    if (LANDSCAPE) "\\end{landscape}")

  L
}

# -- Write the files -----------------------------------------------------------
files <- character(0)
for (m in METHODS) {
  fn <- file.path(TAB_DIR, sprintf("logit_top20reg_%s_%s.tex",
                                   gsub("[^A-Za-z0-9]", "", m), tag_k))
  writeLines(make_table(m), fn)
  files <- c(files, fn)
  cat("written:", fn, "\n")
}

master <- file.path(TAB_DIR, sprintf("logit_top20reg_all_%s.tex", tag_k))
writeLines(c(
  sprintf("%% generated by scripts/04_profiles/04_latex_tables.R -- k = %d", K_CHOSEN),
  "% preamble: \\usepackage{booktabs, longtable, array, colortbl, xcolor, pdflscape}",
  sprintf("\\input{%s}", sub("\\.tex$", "", basename(files)))), master)
cat("master:", master, "\n")

# -- Coding table of the ordinal scales (referenced in the notes) --------------
if (!is.null(ORD_TABLE)) {
  ot <- ORD_TABLE[order(ORD_TABLE$variable), ]
  L <- c(
    "\\begingroup\\footnotesize",
    "\\begin{longtable}{l r l l l}",
    paste("\\caption{Coding of the ordinal scales: the score increases from the first",
          "to the last level.}\\label{tab:ordinal-scales}\\\\"),
    "\\toprule",
    "Variable & Levels & Score 1 & Score max & Neutral (centred at) \\\\",
    "\\midrule \\endfirsthead",
    "\\toprule",
    "Variable & Levels & Score 1 & Score max & Neutral (centred at) \\\\",
    "\\midrule \\endhead",
    "\\midrule \\multicolumn{5}{r}{\\emph{continued on next page}}\\\\ \\endfoot",
    "\\bottomrule \\endlastfoot",
    sprintf("%s & %d & %s & %s & %s \\\\", esc(ot$variable), ot$n_levels,
            esc(ot$first_level), esc(ot$last_level), esc(ot$neutral_level)),
    "\\end{longtable}", "\\endgroup")
  fn <- file.path(TAB_DIR, sprintf("ordinal_scales_%s.tex", tag_k))
  writeLines(L, fn)
  cat("written:", fn, "\n")
}

cat("\nDone. Tables in:", TAB_DIR, "\n")
