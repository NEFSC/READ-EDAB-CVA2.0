#' calculate Pseudo-R2 for BRT
#' @param x is fitted model of class xxx
#' @return 
#' @source https://github.com/elhazen/EcoCast-SciAdv

pseudoR2.brt <- function(x) {
  if ("null.deviance" %in% names(x$self.statistics)){
    d2 <- 1 - (x$self.statistics$resid.deviance / x$self.statistics$null.deviance)
  } else{
    d2 <- 1 - (x$self.statistics$mean.resid / x$self.statistics$mean.null)
  }
  return(d2)
}
