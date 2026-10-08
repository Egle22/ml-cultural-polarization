# =============================================================================
# scripts/05_recovery/05_run_all_blocks.R
# Purpose : Run 04_recovery_all_k.R once per predictor block (dem, dem_geo,
#           pol, pol_geo, pol_dem, pol_dem_geo), then compare the blocks.
#           A batch of 2026-08-26 wasted 9.5 hours rerunning pol_geo six times
#           because a stale copy of the benchmark script was sourced. The
#           pre-flight checks below make that failure mode impossible: the
#           batch refuses to start if more than one copy exists or if the copy
#           lacks the fixes, and aborts on the first block that lands in the
#           wrong directory.
# Thesis  : Appendix A.6; Tables A.11-A.14
# Input   : scripts/05_recovery/04_recovery_all_k.R (sourced per block)
# Output  : results/objects/recovery_all_k/class_<block>/ (via the benchmark)
#           results/objects/recovery_all_k/{cross_block_metrics,
#           geo_incremental_gain,block_overlap}.rds
# Note    : in the archived run the pol_dem_geo block never finished and
#           cross_block_metrics.rds predates the final blocks (see README).
# =============================================================================

source(here::here("R", "00_setup.R"))

out_dir_obj <- file.path(DIR_OBJ, "recovery_all_k")
EXPECTED_VERSION <- "2026-08-27"

# -- Pre-flight 1: exactly one copy of the script -------------------------------
copies <- list.files(here("scripts"), pattern = "^04_recovery_all_k\\.R$",
                     recursive = TRUE, full.names = TRUE)
if (length(copies) == 0L)
  stop("04_recovery_all_k.R not found anywhere under ", here("scripts"))
if (length(copies) > 1L)
  stop("Found ", length(copies), " copies of 04_recovery_all_k.R:\n  ",
       paste(copies, collapse = "\n  "),
       "\nDelete all but one before running the batch.")
SCRIPT <- copies[1]
cat("Script:", SCRIPT, "\n")

# -- Pre-flight 2: the copy carries the fixes -----------------------------------
src <- readLines(SCRIPT, warn = FALSE)
need <- c(
  "BATCH_MODE"        = "BLOCK/USE_GEO would be overwritten on every iteration",
  "tag_rows"          = "Section 14 crashes when LASSO selects nothing",
  "unname(lapply"     = "the method column in per_class_metrics comes out corrupted",
  "BENCHMARK_VERSION" = "no version sentinel"
)
missing_fix <- names(need)[!vapply(names(need),
                                   function(p) any(grepl(p, src, fixed = TRUE)),
                                   logical(1))]
if (length(missing_fix))
  stop("This copy of the script is stale. Missing: ",
       paste(sprintf("%s (%s)", missing_fix, need[missing_fix]),
             collapse = "; "))

ver <- sub('.*"(.*)".*', "\\1",
           grep("^BENCHMARK_VERSION", src, value = TRUE)[1])
if (!identical(ver, EXPECTED_VERSION))
  stop("Script version ", ver, ", expected ", EXPECTED_VERSION)
cat("Version check passed:", ver, "\n")

# ── Runs (cheapest first, so a problem surfaces in ~20 min) ──────────####
runs <- list(
  list(BLOCK = "dem",     USE_GEO = FALSE, key = "dem"),          # dem alone
  list(BLOCK = "dem",     USE_GEO = TRUE,  key = "dem_geo"),      # dem + geo
  list(BLOCK = "pol",     USE_GEO = FALSE, key = "pol"),          # pol alone
  list(BLOCK = "pol",     USE_GEO = TRUE,  key = "pol_geo"),      # pol + geo
  list(BLOCK = "pol_dem", USE_GEO = FALSE, key = "pol_dem"),      # pol + dem
  list(BLOCK = "pol_dem", USE_GEO = TRUE,  key = "pol_dem_geo")   # pol+dem+geo
)

blocks <- c(pol         = "class_only_pol",
            dem         = "class_only_dem",
            pol_geo     = "class_pol_geo",
            dem_geo     = "class_dem_geo",
            pol_dem     = "class_pol_dem",
            pol_dem_geo = "class_pol_dem_geo")
stopifnot(setequal(vapply(runs, `[[`, character(1), "key"), names(blocks)))

# Set SKIP_DONE to TRUE to resume a batch without redoing finished blocks.
SKIP_DONE <- FALSE

# ── Execution ────────────────────────────────────────────────────────####
BATCH_MODE <- TRUE
failed <- character(0)
t_all  <- Sys.time()

for (r in runs) {
  
  done_f <- file.path(out_dir_obj, blocks[[r$key]], "run_meta.rds")
  if (SKIP_DONE && file.exists(done_f) &&
      identical(readRDS(done_f)$block, r$key)) {
    cat("\n-- skipping", r$key, "(already done)\n"); next
  }
  
  cat("\n\n══════  BLOCK:", r$BLOCK, "| geo:", r$USE_GEO,
      "| expected key:", r$key, " ══════\n")
  BLOCK     <- r$BLOCK
  USE_GEO   <- r$USE_GEO
  BLOCK_KEY <- NA_character_          # so a stale value cannot pass the check
  t_run     <- Sys.time()
  
  ok <- tryCatch({
    source(SCRIPT, local = FALSE, echo = FALSE)
    TRUE
  }, error = function(e) {
    message("!! BLOCK ", r$key, " FAILED: ", conditionMessage(e))
    FALSE
  })
  
  if (!ok) { failed <- c(failed, r$key); next }
  
  # Abort the whole batch on a mismatch: continuing would burn hours
  # recomputing the same block, which is exactly what happened before.
  if (!identical(BLOCK_KEY, r$key))
    stop("BLOCK_KEY mismatch: expected '", r$key, "', got '", BLOCK_KEY,
         "'. The sourced script is not honouring BATCH_MODE. Batch aborted.")
  
  m <- readRDS(file.path(out_dir_obj, blocks[[r$key]], "run_meta.rds"))
  if (!identical(m$block, r$key))
    stop("run_meta.rds in ", blocks[[r$key]], " reports block '", m$block,
         "'. Batch aborted.")
  
  cat(sprintf("--> %s OK: %d vars, %d dummies, %.1f min\n",
              BLOCK_KEY, m$n_vars, m$n_dummies,
              as.numeric(difftime(Sys.time(), t_run, units = "mins"))))
}

BATCH_MODE <- FALSE
cat(sprintf("\n\nTotal elapsed: %.1f min\n",
            as.numeric(difftime(Sys.time(), t_all, units = "mins"))))
if (length(failed)) cat("FAILED blocks:", paste(failed, collapse = ", "), "\n")

# ── Provenance ───────────────────────────────────────────────────────####
meta <- bind_rows(lapply(names(blocks), function(b) {
  f <- file.path(out_dir_obj, blocks[[b]], "run_meta.rds")
  if (!file.exists(f))
    return(data.frame(expected = b, found = NA_character_,
                      stringsAsFactors = FALSE))
  m <- readRDS(f)
  data.frame(expected = b, found = m$block, base_block = m$base_block,
             use_geo = m$use_geo, n_vars = m$n_vars, n_dummies = m$n_dummies,
             version = m$script_version %||% NA_character_,
             run_time = format(m$run_time, "%Y-%m-%d %H:%M"),
             stringsAsFactors = FALSE)
}))
cat("\n────────  Run provenance  ────────\n")
print(meta, row.names = FALSE)

stale <- meta$expected[is.na(meta$found) | meta$expected != meta$found]
if (length(stale))
  warning("Blocks missing or written by the wrong run: ",
          paste(stale, collapse = ", "))

# ── Load metrics ─────────────────────────────────────────────────────####
cmp <- bind_rows(lapply(names(blocks), function(b) {
  f <- file.path(out_dir_obj, blocks[[b]], "metrics_long.rds")
  if (file.exists(f)) readRDS(f) else NULL
}))
stopifnot(nrow(cmp) > 0)

BLOCK_ORDER <- names(blocks)
cmp <- cmp %>% mutate(block = factor(block, levels = BLOCK_ORDER))

wide_by <- function(metric, clf = "LASSO") {
  cmp %>%
    filter(classifier == clf) %>%
    select(k, method, block, dplyr::all_of(metric)) %>%
    pivot_wider(names_from = block, values_from = dplyr::all_of(metric)) %>%
    arrange(k, method)
}

# kappa_emp corrects against the empirical marginals: with unbalanced clusters
# it is the honest measure. kappa_unif rewards simply predicting the majority.
cat("\n────────  LASSO — kappa_emp by block  ────────\n")
print(wide_by("kappa_emp"), n = 40, digits = 3)

cat("\n────────  LASSO — accuracy by block  ────────\n")
print(wide_by("acc"), n = 40, digits = 3)

cat("\n────────  RF — kappa_emp by block  ────────\n")
print(wide_by("kappa_emp", "RF"), n = 40, digits = 3)

# ── Incremental value of geography ───────────────────────────────────####
geo_gain <- cmp %>%
  mutate(base = sub("_geo$", "", as.character(block)),
         has_geo = grepl("_geo$", as.character(block))) %>%
  filter(base %in% c("pol", "dem", "pol_dem")) %>%
  select(k, method, classifier, base, has_geo, kappa_emp, acc) %>%
  pivot_wider(names_from = has_geo, values_from = c(kappa_emp, acc),
              names_glue = "{.value}_{ifelse(has_geo, 'geo', 'nogeo')}") %>%
  mutate(d_kappa = kappa_emp_geo - kappa_emp_nogeo,
         d_acc   = acc_geo - acc_nogeo) %>%
  arrange(base, k, method, classifier)

cat("\n────────  Gain from adding geography (LASSO)  ────────\n")
geo_gain %>% filter(classifier == "LASSO") %>%
  select(base, k, method, kappa_emp_nogeo, kappa_emp_geo, d_kappa) %>%
  print(n = 60, digits = 3)

cat("\n────────  Mean gain from geography, by base block  ────────\n")
geo_gain %>%
  group_by(base, classifier) %>%
  summarise(mean_d_kappa = mean(d_kappa, na.rm = TRUE),
            mean_d_acc   = mean(d_acc,   na.rm = TRUE),
            n = dplyr::n(), .groups = "drop") %>%
  arrange(base, classifier) %>%
  print(digits = 3)

# ── Do the two blocks overlap? pol vs dem vs pol_dem ─────────────────####
overlap <- cmp %>%
  filter(classifier == "LASSO", block %in% c("pol", "dem", "pol_dem")) %>%
  select(k, method, block, kappa_emp) %>%
  pivot_wider(names_from = block, values_from = kappa_emp) %>%
  mutate(best_single = pmax(pol, dem, na.rm = TRUE),
         gain_from_combining = pol_dem - best_single) %>%
  arrange(k, method)

cat("\n────────  Combined block vs. better single block  ────────\n")
print(overlap, n = 40, digits = 3)
cat("\nMean gain from combining:",
    round(mean(overlap$gain_from_combining, na.rm = TRUE), 4), "\n")

saveRDS(cmp,      file.path(out_dir_obj, "cross_block_metrics.rds"))
saveRDS(geo_gain, file.path(out_dir_obj, "geo_incremental_gain.rds"))
saveRDS(overlap,  file.path(out_dir_obj, "block_overlap.rds"))
cat("\nComparison objects saved in:", out_dir_obj, "\n")