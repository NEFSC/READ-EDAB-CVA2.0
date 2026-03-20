#' @title Calculate Certainty on Bootstrapped Sensitivity Scores
#' @description Calculates certainty - the percentage of bootstrapped values that match expert scores - on sensitivity

#' @param bootstrap.scores output from \code{calculate.sensitivity(bootstrap = TRUE)}
#' @param expert.scores output from \code{calculate.sensitivity(bootstrap = FALSE)} to append bootstrap certainty to.

#' @return returns the \code{expert.scores} with a new column called 'Certainty', which includes the percentage of bootstrapped sensitivities that matched the weighted average final sensitivity.


bootstrap.certainty <- function(bootstrap.scores, expert.scores){

  s <- expert.scores$`Total Sensitivity` #pull the value of 1,2,3,4 from expert-derived scores
  expert.scores$Certainty <- length(which(bootstrap.scores[,ncol(bootstrap.scores)]==s))/nrow(bootstrap.scores) #calculate the percentage of scores that equal the expert-derived score

  return(expert.scores)
}
