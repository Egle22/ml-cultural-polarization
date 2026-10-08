# =============================================================================
# scripts/01_data/02_missing_analysis.R
# Purpose : Count missing values by questionnaire section and by type of
#           non-response (don't know, refused, not applicable, ...).
# Thesis  : Appendix A.2, Table A.2
# Input   : data/processed/data_raw_questions.rds
# Output  : data/processed/all_missing.rds
#           data/processed/table_miss_per_sec.rds
#           data/processed/miss_per_question.rds
# =============================================================================

source(here::here("R", "00_setup.R"))

data <- readRDS(file.path(DIR_DATA, "data_raw_questions.rds"))

# Building sections
df_section <- build_sections(
  data = data,
  questions = questions,
  section_labels = section_labels
)

# Missing analysis per section
all_missing <- mapply(
  analyze_missing_section,
  df = df_section,
  section_name = names(df_section),
  SIMPLIFY = FALSE
) |> do.call(what = rbind)

rownames(all_missing) <- NULL

# Summary table
total_missing <- sum(all_missing$n_miss_tot)

table_miss_per_sec2 <- all_missing %>%
  group_by(sec) %>%
  summarise(
    n_var = n(),
    miss_per_sec = paste0(
      round(sum(n_miss_tot) / (sum(n_miss_tot) + sum(n_tot)) * 100, 1),
      "%"
    ),
    miss_tot = paste0(
      round(sum(n_miss_tot) / total_missing * 100, 1),
      "%"
    ),
    across(starts_with("n_") | starts_with("X."), sum),
    .groups = "drop"
  ) %>%
  rename(
    dont_know      = X.1,
    refused        = X.2,
    not_applicable = X.3,
    not_asked      = X.4,
    missing_other  = X.5
  ) %>%
  arrange(desc(n_miss_tot))

# Missing values per variable
miss_per_question <- all_missing %>%
  filter(pct_miss > 1) %>%
  arrange(sec, desc(pct_miss)) %>%
  rename(
    dont_know      = X.1,
    refused        = X.2,
    not_applicable = X.3,
    not_asked      = X.4,
    missing_other  = X.5
  )

print(as.data.frame(miss_per_question))

# Save results
saveRDS(all_missing,         file.path(DIR_DATA, "all_missing.rds"))
saveRDS(table_miss_per_sec2, file.path(DIR_DATA, "table_miss_per_sec.rds"))
saveRDS(miss_per_question,   file.path(DIR_DATA, "miss_per_question.rds"))
