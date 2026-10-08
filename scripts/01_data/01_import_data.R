# =============================================================================
# scripts/01_data/01_import_data.R
# Purpose : Import the raw WVS Wave 7 U.S. file and keep the questionnaire
#           items listed in R/constants.R (`questions`).
# Thesis  : Section 2.1
# Input   : data/raw/WVS_Wave_7_United_States_Csv_v5.1.csv (not redistributed)
# Output  : data/processed/data_raw_questions.rds
#           data/processed/df_section.rds  (list of data frames, one per section)
# =============================================================================

source(here::here("R", "00_setup.R"))

data <- read_delim(
  file.path(DIR_RAW, "WVS_Wave_7_United_States_Csv_v5.1.csv"),
  delim = ";",
  escape_double = FALSE,
  trim_ws = TRUE
)

Q_var <- unlist(questions, use.names = FALSE)
data <- data %>% dplyr::select(all_of(Q_var))

df_section <- build_sections(data, questions, section_labels)

# Political and sociodemographic items (held out from the clustering stage)
Political    <- c(names(df_section$`Political Interest`), names(df_section$`Political Culture`))
Demographics <- names(df_section$Demographics)

str(data)
saveRDS(data,       file.path(ensure_dir(DIR_DATA), "data_raw_questions.rds"))
saveRDS(df_section, file.path(DIR_DATA, "df_section.rds"))
