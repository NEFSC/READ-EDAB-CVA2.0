#' @title Calculate Sensitivity
#' @description Wrapper function that calculate species sensitivity using \code{attribute.score} and \code{logic.rule} for individual attributes and apply logic rule to calculate total sensitivity. The option to perform bootstrapping is provided. Based on code from the Southeast Fisheries Science Center (SEFSC)

#' @param species.attributes a list of data frames, with each data frame corresponding to all of the expert scores for one species. See example for how to generate this list from the FCVA output. 
#' @param bootstrap binary TRUE/FALSE to turn on/off bootstrapping of sensitivity scores. Will impact outputs. 
#' @param samples number of samples to run for bootstrapping. Default is 10,000.

#' @return a data.frame where each row is a species, and columns represent attribute scores, as well as the Total Sensitivity score column returned by the logic rule

#' @examples
#' data <-read.csv('expert_scores.csv') #must contain columns called 'Stock.Name' and 'Attribute.Name' and should only contain sensitivity data
#' library(plyr)
#' species.data.list <- dlply(data, .(Stock.Name))
#' species.sensitivity <- lapply(species.data.list, calculate.sensitivity, bootstrap = F)
#' sensitivity.bootstrap <- lapply(species.data.list, calculate.sensitivity, bootstrap = T)
#' sensitivity.certainty <- lapply(sensitivity.bootstrap, bootstrap.certainty, match = T, sensitivity.dataframe = species.sensitivity)

calculate.sensitivity <- function(species.attributes, bootstrap = TRUE, samples = 10000){
  #species.attributes is a list of dataframes, with one data frame per species
  #bootstrap is boolean TRUE/FALSE to trigger bootstrapping 
  #samples is the number of samples to use in bootstrapping 
  
  #subset by attribute
  attribute.list <- dlply(species.attributes, .(Attribute.Name)) #make list of data frames with one data frame per list
  a.scores <- lapply(attribute.list, attribute.score, bootstrap = bootstrap, samples = samples) #calculate scores for attributes via weighted mean
  out <- do.call(cbind, a.scores) #combine scores into vector
  
  #apply logic rule and add to vector 
  total <- logic.rule(out, bootstrap = bootstrap)
  
  #clean up vector
  sensitivity.attributes <- as.data.frame(cbind(out, total))
  colnames(sensitivity.attributes)[dim(sensitivity.attributes)[2]] <- 'Total Sensitivity'
  return(sensitivity.attributes)
}
