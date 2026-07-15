#' @title Remove Correlated Environmental Covariates from Species/Environmental Data Frame
#' @description
#' This function removed correlated covariates to prepare data for modeling
#'
#' @param df data.frame of fisheries and environmental data
#' @param var_names a vector containing the names of the environmental variables to check for correlation
#'
#' @return a data frame with correlated covariates removed
#'
#'@export

remove_corr <- function(df, var_names) {
  ind <- names(df) %in% var_names
  corInd <- caret::findCorrelation(stats::cor(df[, ind]), names = T) #find correlated variables
  if (length(corInd) != 0) {
    df <- df[, -which(colnames(df) == corInd)]
  }
  return(df)
}
