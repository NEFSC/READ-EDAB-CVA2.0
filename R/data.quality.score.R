#' @title Calculate Data Quality Scores for each Attribute & Species
#' @description Calculates mean data quality score

#' @param species.attributes a data frame containing all of the expert scores for one species. See example for how to generate this list from the FCVA output.

#' @return a vector of mean data quality scores for each attribute with the length of the number of attributes.
#'
#' @examples
#' \dontrun{
#' data <-read.csv('expert_scores.csv') #sensitivity data from FSCVA Portal
#' species.data.list <- split(data, data$Stock.Name)
#'
#' species.dqs <- lapply(species.data.list, data.quality.score)
#'
#' speciesDQ <- do.call(rbind, species.dqs)
#' #this combines the list of vectors into a data.frame
#' }

data.quality.score <- function(species.attributes){
  #calculate mean data quality
  dq <- aggregate(Data.Quality ~ Attribute.Name, data = species.attributes, FUN = mean, na.rm = T)
  #format data to fit into eventual csv
  dq <- as.data.frame(t(dq[,2]))
  colnames(dq) <- unique(species.attributes$Attribute.Name)
  rownames(dq) <- unique(species.attributes$Stock.Name)
  return(dq)
}
