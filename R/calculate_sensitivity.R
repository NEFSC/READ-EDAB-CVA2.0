#' @title Calculate Sensitivity
#' @description Wrapper function that calculate species sensitivity using \code{attribute_score} and \code{apply_logic_rule} for individual attributes and apply logic rule to calculate total sensitivity. The option to perform bootstrapping is provided. Based on code from the Southeast Fisheries Science Center (SEFSC)
#'
#'
#' @param species_attributes a list of data frames, with each data frame corresponding to all of the expert scores for one species. See example for how to generate this list from the FCVA output.
#' @param bootstrap binary TRUE/FALSE to turn on/off bootstrapping of sensitivity scores. Will impact outputs.
#' @param samples number of samples to run for bootstrapping. Default is 10,000.

#' @return a data.frame where each row is a species, and columns represent attribute scores, as well as the Total Sensitivity score column returned by the logic rule

#' @examples
#' \dontrun{
#' data <-read.csv('expert_scores.csv') #sensitivity data from FSCVA Portal
#' species.data.list <- split(data, data$Stock.Name)
#' species.sensitivity <- lapply(species.data.list, calculate.sensitivity, bootstrap = F)
#' speciesDF <- do.call(rbind, species.sensitivity)
#' sensitivity.bootstrap <- lapply(species.data.list, calculate.sensitivity, bootstrap = T)
#' sensitivity.certainty <- mapply(bootstrap.certainty, sensitivity.bootstrap,
#'  species.sensitivity, SIMPLIFY = F)
#' #note the use of mapply rather than lapply here since both inputs are lists
#since both input lists are use species.data.list, they will be in the same order
#' sensitivityDF <- do.call(rbind, sensitivity.certainty)
#' #this combines the list of vectors into a data.frame
#' }

calculate_sensitivity <- function(species_attributes, bootstrap = TRUE, samples = 10000){
  #species_attributes is a list of dataframes, with one data frame per species
  #bootstrap is boolean TRUE/FALSE to trigger bootstrapping
  #samples is the number of samples to use in bootstrapping

  #subset by attribute
  attribute.list <- split(species_attributes, species_attributes$Attribute.Name) #make list of data frames with one data frame per list
  a.scores <- lapply(attribute.list, attribute_score, bootstrap = bootstrap, samples = samples) #calculate scores for attributes via weighted mean
  out <- do.call(cbind, a.scores) #combine scores into vector

  #apply logic rule and add to vector
  total <- apply_logic_rule(out, bootstrap = bootstrap)

  #clean up vector
  sensitivity.attributes <- as.data.frame(cbind(out, total))
  colnames(sensitivity.attributes)[dim(sensitivity.attributes)[2]] <- 'Total Sensitivity'
  return(sensitivity.attributes)
}
