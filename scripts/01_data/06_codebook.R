# =============================================================================
# scripts/01_data/06_codebook.R
# Purpose : Build the variable codebook and print it as an HTML table. The
#           printed table is the GitHub Pages codebook (docs/index.html): save
#           it from the RStudio Viewer ("Export > Save as Web Page") to
#           docs/index.html. It contains only per-variable aggregate shares.
# Thesis  : Appendix A.1
# Input   : data/processed/data_ord_fct.rds
# Output  : data/processed/codebook.rds, data/processed/codebook.csv
# =============================================================================

source(here::here("R", "00_setup.R"))

data_ord_fct <- readRDS(file.path(DIR_DATA, "data_ord_fct.rds"))

codebook <- build_codebook(
  data = data_ord_fct,
  questions = questions,
  section_labels = section_labels
)

saveRDS(codebook, file.path(DIR_DATA, "codebook.rds"))
write.csv(codebook, file.path(DIR_DATA, "codebook.csv"), row.names = FALSE)

print_codebook_table(
  codebook = codebook,
  section_labels = section_labels
)
