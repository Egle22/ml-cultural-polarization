# =============================================================================
# scripts/01_data/03_cleaning.R
# Purpose : Light cleaning used for the codebook: drop respondents with
#           systematic non-response and demographic items with structural
#           missingness. (The analysis sample is built in 04_preprocess_impute.R.)
# Thesis  : Appendix A.1 (codebook), A.2
# Input   : data/processed/data_raw_questions.rds, data/processed/df_section.rds
# Output  : data/processed/missing_per_subject.rds
#           data/processed/miss_per_question_dem.rds
#           data/processed/data_clean.rds
# =============================================================================

source(here::here("R", "00_setup.R"))

data       <- readRDS(file.path(DIR_DATA, "data_raw_questions.rds"))
df_section <- readRDS(file.path(DIR_DATA, "df_section.rds"))

# -- Identify subjects with unusually many missing values ---------------------
missing_per_subj <- data %>%
  mutate(
    n_m1 = rowSums(. == -1, na.rm = TRUE),
    n_m2 = rowSums(. == -2, na.rm = TRUE),
    n_m3 = rowSums(. == -3, na.rm = TRUE),
    n_m5 = rowSums(. == -5, na.rm = TRUE),
    n_miss_tot = n_m1 + n_m2 + n_m3 + n_m5
  ) %>%
  select(n_m1, n_m2, n_m3, n_m5, n_miss_tot) %>%
  mutate(row_id = row_number())

subj_72  <- missing_per_subj %>% filter(n_m1 == 72)
subj_272 <- missing_per_subj %>% filter(n_m2 == 272)
subj_275 <- missing_per_subj %>% filter(n_miss_tot == 275)

d_rows <- unique(c(subj_72$row_id, subj_272$row_id, subj_275$row_id))

data_clean <- data %>%
  mutate(row_id = row_number()) %>%
  filter(!row_id %in% d_rows) %>%
  select(-row_id)

# -- Analyse missing values in demographics -----------------------------------
Demography_q <- df_section$Demographics

demo_missing <- single_section_missing(Demography_q)

miss_per_question_dem <- demo_missing %>%
  arrange(desc(pct_miss)) %>%
  select(-X.3, -X.4, -X.5) %>%
  rename(
    dont_know = X.1,
    refused   = X.2
  )

# -- Remove demographic variables with structural missingness -----------------
data_clean <- data_clean %>%
  select(-Q282, -Q276, -Q280, -Q281)

# -- Save cleaned data --------------------------------------------------------
saveRDS(missing_per_subj,      file.path(DIR_DATA, "missing_per_subject.rds"))
saveRDS(miss_per_question_dem, file.path(DIR_DATA, "miss_per_question_dem.rds"))
saveRDS(data_clean,            file.path(DIR_DATA, "data_clean.rds"))
