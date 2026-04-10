#' @title Calculate Certainty on Bootstrapped Sensitivity Scores
#' @description Calculates certainty - the percentage of bootstrapped values that match expert scores - on sensitivity

#' @param bootstrap_scores output from \code{calculate_sensitivity(bootstrap = TRUE)}
#' @param expert_scores output from \code{calculate_sensitivity(bootstrap = FALSE)} to append bootstrap certainty to.

#' @return returns the \code{expert_scores} with a new column called 'Certainty', which includes the percentage of bootstrapped sensitivities that matched the weighted average final sensitivity.

calculate_sensitivity_certainty <- function(bootstrap_scores, expert_scores) {
  s <- expert_scores$`Total Sensitivity` #pull the value of 1,2,3,4 from expert-derived scores
  expert_scores$Certainty <- length(which(
    bootstrap_scores[, ncol(bootstrap_scores)] == s
  )) /
    nrow(bootstrap_scores) #calculate the percentage of scores that equal the expert-derived score

  return(expert_scores)
}
