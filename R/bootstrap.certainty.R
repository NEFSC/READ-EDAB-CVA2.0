#' @title Calculate Certainty on Bootstrapped Sensitivity Scores
#' @description Calculates certainty - the percentage of bootstrapped values that match expert scores - on sensitivity 

#' @param species.bootstrap output from \code{calculate.sensitivity(bootstrap = TRUE)}
#' @param match binary TRUE/FALSE to match bootstrap certainty (percentage of bootstrapped scores that match weighted averaged score) to output from \code{calculate.sensitivity(bootstrap = FALSE)}
#' @param sensitivity.dataframe The output from \code{calculate.sensitivity(bootstrap = FALSE)} to append bootstrap certainty to. Must be provided if \code{match = TRUE}. 

#' @return returns a data.frame or a vector. If \code{match = TRUE}, it returns the \code{sensitivity.dataframe} with a new column called 'Certainty', which includes the percentage of bootstrapped sensitivities that matched the weighted average final sensitivity. If \code{match = FALSE}, it returns a vector containing what percent of bootstrapped values were equal to each possible sensitivity score (1-4).


bootstrap.certainty <- function(species.bootstrap, match = T, sensitivity.dataframe){
  #species.bootstrap is the bootstrap matrix that results from calculate.sensitivity(bootstrap=T)
  
  percentages <- vector(length = 4)
  # Go through the data for each species and calculate the portion of scores in each category
  percentages[1] <- length(which(species.bootstrap[,ncol(species.bootstrap)]==1))/nrow(species.bootstrap)
  percentages[2] <- length(which(species.bootstrap[,ncol(species.bootstrap)]==2))/nrow(species.bootstrap)
  percentages[3] <- length(which(species.bootstrap[,ncol(species.bootstrap)]==3))/nrow(species.bootstrap)
  percentages[4] <- length(which(species.bootstrap[,ncol(species.bootstrap)]==4))/nrow(species.bootstrap)
  
  if(match){
    sensitivity.dataframe$Certainty <- 0
    for(x in 1:nrow(sensitivity.dataframe)){
      s <- sensitivity.dataframe$`Total Sensitivity`[x]
      sensitivity.dataframe$Certainty[x] <- certaintyDF[x,s]
      return(sensitivity.dataframe)
    }
  } else {
    return(percentages)
  }
}
