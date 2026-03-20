#' @title Calculate Certainty on Bootstrapped Directionality Scores
#' @description Calculates certainty - the percentage of bootstrapped values that match expert scores - on directionality

#' @param bootstrap.scores output from \code{directionality(bootstrap = TRUE)}
#' @param directionality.scores output from \code{directionality(bootstrap = FALSE)} to append bootstrap certainty to.

#' @return returns a vector with a length of two: 1) the output from \code{directionality(bootstrap = FALSE)} and 2) the certainty value, which is the percentage of bootstrapped sensitivities that matched the weighted average final directionality


directionality.certainty <- function(bootstrap.scores, directionality.scores){

  cert <- length(which(bootstrap.scores==directionality.scores))/length(bootstrap.scores) #calculate the percentage of scores that equal the expert-derived score
  directionality.scores <- c(directionality.scores, cert)
  names(directionality.scores) <- c("Directionality", 'Certainty')

  return(directionality.scores)
}
