#' @title calculate Pseudo-R2 for BRT
#' @description
#' Used as part of \code{eval_brt}. From Camrin Brawn (WHOI): https://zenodo.org/records/7971532. Included with CVA2.0 package to ensure functionality.
#'
#' @param x is fitted model of class xxx
#' @return value representing Pseudo-R2
#'
#' @source https://github.com/elhazen/EcoCast-SciAdv

pseudo_r2_brt <- function(x) {
  if ("null.deviance" %in% names(x$self.statistics)){
    d2 <- 1 - (x$self.statistics$resid.deviance / x$self.statistics$null.deviance)
  } else{
    d2 <- 1 - (x$self.statistics$mean.resid / x$self.statistics$mean.null)
  }
  return(d2)
}
