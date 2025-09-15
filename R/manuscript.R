remove_dups <- function(x) {
  x <- as.character(x)
  x[x == lag(x)] = ""
  x
}

bold_best <- function(df, col1, col2, target = 0) {
  col1 <- as.character(ensym(col1))
  col2 <- as.character(ensym(col2))

  if(is.null(df[[col1]])) stop(glue::glue("Column {col1} does not exist"))
  if(is.null(df[[col2]])) stop(glue::glue("Column {col2} does not exist"))

  n1 <- parse_number(df[[col1]])
  n2 <- parse_number(df[[col2]])
  
  d1 <- abs(n1 - target)
  d2 <- abs(n2 - target)

  df[, col1] <- cell_spec(df[[col1]], bold = ifelse(is.na(d1 < d2), FALSE, d1 < d2), format = "latex") 
  df[, col2] <- cell_spec(df[[col2]], bold = ifelse(is.na(d2 < d1), FALSE, d2 < d1), format = "latex") 
  
  df
}
