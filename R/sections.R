# ------------------------------------------------------------------------------
# build_sections()
#
# Create a list of datasets corresponding to each thematic section of the survey.
#
# The function uses the 'questions' list, where each element contains the
# variable names belonging to a specific section of the questionnaire.
#
# For each section:
#   - the corresponding variables are selected from the dataset
#   - a dataframe containing only those variables is created
#
# The output is a named list of dataframes, where:
#   - each element represents a survey section
#   - names are taken from 'section_labels' for readability
#
# Example output structure:
#   sections
#     ├─ Social Values
#     ├─ Happiness & Wellbeing
#     ├─ Social Capital & Trust
#     ├─ Economic Values
#     └─ ...
#
# This function simplifies the creation of section-specific datasets and avoids
# manually writing multiple select() statements.

build_sections <- function(data, questions, section_labels) {
  setNames(
    lapply(names(questions), function(s) dplyr::select(data, dplyr::any_of(questions[[s]]))),
    section_labels
  )
}
