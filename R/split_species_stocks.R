#' @title Split Species-Stock columns into two columns
#' @description
#' Separates species and stocks by deliminator and creates two columns within the provided dataset
#'
#' @param df a data.frame with the column to be split
#' @param species_stock_col the column name containing both the species and stock name if available. The function handles species without stocks.
#' @param delim the deliminator used to seperate species and stocks. Should be unique and not used elsewhere in the string
#'
#' @return the provided data frame with the columns 'Species' and 'Stocks' added
#'
#'@export
#'

split_species_stocks <- function(df, species_stock_col, delim = "-\\s*") {
  # Force column into a 2-column matrix
  # "-\\s*" means: match a hyphen, followed by zero or more spaces
  split_mat <- stringr::str_split_fixed(
    df[, species_stock_col],
    pattern = delim,
    n = 2
  )

  # Convert empty strings to NA for cleaner data handling
  split_mat[split_mat == ""] <- NA

  # Extract as separate vectors if needed
  species <- split_mat[, 1]
  stock <- split_mat[, 2]

  df <- cbind(species, stock, df)
  colnames(df)[1:2] <- c('Species', 'Stock')

  return(df)
}
