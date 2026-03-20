#' @title Calculate Model Confidence Score for each Species
#' @description Calculates mean and standard deviation of the model confidence score

#' @param species a data frame containing all of the expert scores for one species in a column called \code{Score}. See example for how to generate this list from the combined model confidence score spreadsheets.

#' @return a data frame including mean and standard deviation of data quality scores for each attribute with the length of the number of attributes.
#'
#' @examples
#' \dontrun{
#' species.data.list <- split(data, data$Species)
#' #data is the combined data frame of
#' #all model confidence spreadsheets from scorers
#'
#' species.conf <- lapply(species.data.list, model.confidence)
#'
#' speciesMC <- do.call(rbind, species.conf)
#' #this combines the list of vectors into a data.frame
#' }

model.confidence <- function(species){
  #calculate mean data quality
  meanConf <- mean(species$Score, na.rm = T)
  sdConf <- sd(species$Score, na.rm = T)
  #format data to fit into eventual csv
  cn <- as.data.frame(Species = species$Species[1], meanConfidence = meanConf, sdConfidence = sdConf)
  return(cn)
}
