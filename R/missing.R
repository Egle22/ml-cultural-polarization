# getting the mode 
get_mode <- function(x) {
  ux <- na.omit(unique(x))
  ux[which.max(tabulate(match(x, ux)))]
}



# Missing values per columns and type 
analyze_missing_section <- function(df, section_name) {
  n_rows <- nrow(df)
  lapply(names(df), function(col) {
    x      <- df[[col]]
    counts <- sapply(missing_codes, function(code) sum(x == code, na.rm = TRUE))
    names(counts) <- names(missing_labels)
    data.frame(
      sec    = section_name,
      var  = col,
      n_tot   = sum(!is.na(x) & !(x %in% missing_codes)),
      n_miss_tot = sum(counts) + sum(is.na(x)),
      pct_miss   = round((sum(counts) + sum(is.na(x))) / n_rows * 100, 1),
      t(counts),  
      row.names = NULL
    )
  }) |> do.call(what = rbind)
}




# Missing value per subject 
missing_per_row <- function(df) {
  counts <- sapply(missing_codes, function(code) rowSums(df == code, na.rm = TRUE))
  colnames(counts) <- paste0("n_", names(missing_codes))
  as.data.frame(counts) %>% mutate(n_miss_tot = rowSums(.), row_id = row_number())
}


#single section missing 
single_section_missing <- function(df) {
  n_rows <- nrow(df)
  lapply(names(df), function(col) {
    x      <- df[[col]]
    counts <- sapply(missing_codes, function(code) sum(x == code, na.rm = TRUE))
    names(counts) <- names(missing_labels)
    data.frame(
      var  = col,
      n_tot   = sum(!is.na(x) & !(x %in% missing_codes)),
      n_miss_tot = sum(counts) + sum(is.na(x)),
      pct_miss   = round((sum(counts) + sum(is.na(x))) / n_rows * 100, 1),
      t(counts),  
      row.names = NULL
    )
  }) |> do.call(what = rbind)
  
}

# Refused or don't know 
ref_dont <- function(df) {
  as.data.frame(lapply(df, function(col) {
    factor(case_when(
      col == -1 ~ "dont_know",
      col == -2 ~ "refused",
      TRUE      ~ "0"
    ), levels = c("0", "dont_know", "refused"))
  }))
}
