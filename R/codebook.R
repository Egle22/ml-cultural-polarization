# ------------------------------------------------------------------------------
# build_codebook()
#
# Create a variable-level codebook from a labelled dataset.
#
# This function scans all variables in the dataset and builds a summary table
# describing their properties. For each variable it:
#
#  • identifies the questionnaire section using the 'questions' list
#  • determines the variable type (numeric, ordered factor, nominal factor)
#  • counts the number of missing values
#  • computes the percentage of missing values
#  • summarises the distribution of valid responses
#
# For factor variables:
#   the distribution is reported as category percentages.
#   categories with ≥50% frequency are highlighted in bold.
#
# For numeric variables:
#   the distribution is summarised using min / mean / max.
#
# INPUT
#   data            labelled dataset (e.g. data_ord_fct)
#   questions       list mapping questionnaire sections to variable names
#   section_labels  human-readable section names (same order as questions)
#
# OUTPUT
#   dataframe with one row per variable containing:
#     Section | Variable | Type | N_missing | Pct_missing | Distribution
#
# USAGE
#   codebook <- build_codebook(data_ord_fct, questions, section_labels)
#
# ------------------------------------------------------------------------------

build_codebook <- function(data, questions, section_labels,
                           miss_lbls = c("Missing/NA", "Not applicable",
                                         "No answer/Refused", "Don't know")) {
  
  # Map each variable to its questionnaire section
  var_section <- data.frame(
    var     = unlist(questions, use.names = FALSE),
    section = rep(section_labels, times = sapply(questions, length)),
    stringsAsFactors = FALSE
  )
  
  # Build one row of the codebook for each variable
  rows <- lapply(names(data), function(col) {
    
    x   <- data[[col]]
    sec <- var_section$section[var_section$var == col]
    sec <- if (length(sec) == 0) "Unknown" else sec[1]
    
    n_total <- length(x)
    n_na    <- sum(is.na(x))
    
    if (is.factor(x)) {
      n_miss <- sum(x %in% miss_lbls, na.rm = TRUE) + n_na
      type   <- if (is.ordered(x)) "factor (ordered)" else "factor (nominal)"
      
      valid_x <- x[!x %in% miss_lbls & !is.na(x)]
      tbl     <- sort(table(droplevels(valid_x)), decreasing = TRUE)
      n_valid <- sum(tbl)
      
      distrib <- if (n_valid > 0) {
        paste(
          sapply(seq_along(tbl), function(i) {
            pct <- round(as.numeric(tbl[i]) / n_valid * 100, 1)
            lbl <- paste0(names(tbl)[i], " = ", pct, "%")
            if (pct >= 50) paste0("<b>", lbl, "</b>") else lbl
          }),
          collapse = " | "
        )
      } else {
        "—"
      }
      
    } else {
      x_num   <- suppressWarnings(as.numeric(as.character(x)))
      x_valid <- x_num[!is.na(x_num) & x_num >= 0]
      n_miss  <- sum(is.na(x_num) | x_num < 0, na.rm = TRUE)
      type    <- "numeric"
      
      distrib <- if (length(x_valid) > 0) {
        paste0(
          "min=", min(x_valid),
          " | mean=", round(mean(x_valid), 1),
          " | max=", max(x_valid)
        )
      } else {
        "—"
      }
    }
    
    data.frame(
      Section      = sec,
      Variable     = col,
      Type         = type,
      N_missing    = n_miss,
      Pct_missing  = round(n_miss / n_total * 100, 1),
      Distribution = distrib,
      stringsAsFactors = FALSE
    )
  })
  
  codebook <- do.call(rbind, rows)
  
  # Order sections according to the predefined order
  codebook$Section <- factor(codebook$Section, levels = section_labels)
  codebook <- codebook[order(codebook$Section), ]
  
  rownames(codebook) <- NULL
  
  return(codebook)
}


# ------------------------------------------------------------------------------
# print_codebook_table()
#
# Print a formatted HTML table of the codebook using kableExtra.
#
# The function:
#   • formats the codebook dataframe
#   • groups variables by questionnaire section
#   • preserves the section order defined in 'section_labels'
#
# It is used only for visualization / reporting purposes.
#
# INPUT
#   codebook        dataframe produced by build_codebook()
#   section_labels  section order used for grouping rows
#
# OUTPUT
#   styled HTML table.
#
# USAGE
#   print_codebook_table(codebook, section_labels)
#
# ------------------------------------------------------------------------------ 

print_codebook_table <- function(codebook, section_labels,
                                 caption = "WVS Wave 7 USA – Variable Codebook") {
  
  section_counts <- as.integer(table(codebook$Section)[section_labels])
  names(section_counts) <- section_labels
  section_counts <- section_counts[section_counts > 0]
  
  codebook %>%
    kbl(
      caption   = caption,
      col.names = c("Section", "Variable", "Type", "N missing", "% missing", "Distribution"),
      align     = c("l", "l", "l", "r", "r", "l"),
      booktabs  = TRUE,
      escape    = FALSE
    ) %>%
    kable_styling(
      bootstrap_options = c("striped", "hover", "condensed"),
      full_width        = FALSE,
      font_size         = 11
    ) %>%
    column_spec(6, width = "35em") %>%
    pack_rows(index = section_counts)
}

