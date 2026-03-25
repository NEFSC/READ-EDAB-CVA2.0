#' @title Calculate Certainty on Bootstrapped Directionality Scores
#' @description Calculates certainty - the percentage of bootstrapped values that match expert scores - on directionality

#' @param bootstrap_scores output from \code{calculate_directionality(bootstrap = TRUE)}
#' @param directionality_scores output from \code{calculate_directionality(bootstrap = FALSE)} to append bootstrap certainty to.

#' @return returns a vector with a length of two: 1) the output from \code{calculate_directionality(bootstrap = FALSE)} and 2) the certainty value, which is the percentage of bootstrapped sensitivities that matched the weighted average final directionality


directionality_certainty <- function(bootstrap_scores, directionality_scores){

  cert <- length(which(bootstrap_scores==directionality_scores))/length(bootstrap_scores) #calculate the percentage of scores that equal the expert-derived score
  directionality_scores <- c(directionality_scores, cert)
  names(directionality_scores) <- c("Directionality", 'Certainty')

  return(directionality_scores)
}
