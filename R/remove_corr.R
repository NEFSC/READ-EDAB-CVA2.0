#' @title Remove Correlated Environmental Covariates from Species/Environmental Data Frame
#' @description
#' This function removed correlated covariates to prepare data for modeling
#'
#' @param se data.frame of fisheries and environmental data
#' @param pa_col column name for presence/absence column
#' @param xy_col a vector with a length of 2 indicating the longitude and latitude column names
#' @param month_col,year_col column names for month and year columns respectively
#'
#' @return a data frame with correlated covariates removed

remove_corr <- function(se, pa_col, xy_col, month_col, year_col) {
  ind <- which(
    colnames(se) == pa_col |
      colnames(se) == xy_col[1] |
      colnames(se) == xy_col[2] |
      colnames(se) == month_col |
      colnames(se) == year_col
  )
  corInd <- caret::findCorrelation(cor(se[, -ind]), names = T) #find correlated variables
  if (length(corInd) != 0) {
    se <- se[, -which(colnames(se) == corInd)]
  }
  return(se)
}
