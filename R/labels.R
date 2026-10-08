# =============================================================================
# R/labels.R
# Factor-labelling helpers used to turn WVS numeric codes into labelled
# factors (scripts/01_data/04_preprocess_impute.R and 05_label_cleaned_data.R).
# Restored verbatim from the version that produced the saved data
# (R/functions_labels.R at commit 99bd27f).
#
# NOTE: fct_ord() returns ORDERED factors. They are only used for the labelled
# dataset (imputed_labeled.rds). Every modelling script converts them back to
# unordered factors or numeric scores before LASSO/logit: never pass an ordered
# factor to glmnet/glm (polynomial contrasts).
# =============================================================================

fct <- function(x, v, l) {
  factor(x, levels = c(miss_lvl, v), labels = c(miss_lbl, l))
}

fct_ord <- function(x, v, l) {
  factor(x, levels = c(miss_lvl, v), labels = c(miss_lbl, l), ordered = TRUE)
}

fct_10 <- function(x) {
  fct_ord(x, 1:10, as.character(1:10))
}
