#' @title Make Timeseries of Total Species-Specific Exposure
#'
#' @description Uses the logic rule from CVA1.0 to combine variable-specific exposures to create a total exposure timeseries
#'
#' @param matExp matrix output from \code{makeExposureTimeseries}
#' @param countAll TRUE/FALSE to use \code{weights} and \code{wThreshold} to subset variables to only important variables
#' @param weights output from \code{combineWeights} - a vector of variable weights in ensemble SDM
#' @param wThreshold numeric value used to subset weights, variables with weights less to or equal to this value will be excluded from total exposure calculation
#'
#' @return A single vector representing total exposure across time

combineTimeseries <- function(matExp, countAll, weights, wThreshold){
  if(countAll){ #if counting all included factors and not taking weight into account
    matSub <- matExp
  } else {
    wi <- which(weights > wThreshold)
    matSub <- matExp[wi,]
  }

  #count each rank in each column
  hr <- hh <- md <- vector(length = ncol(matSub))
  for(x in 1:ncol(matSub)){
    hr[x] <- length(which(matSub[,x] >= 3.5))
    hh[x] <- length(which(matSub[,x] >= 3))
    md[x] <- length(which(matSub[,x] >= 2.5))
  }

  #apply logic rule
  expV <- rep(1, times = ncol(matSub))
  expV <- replace(expV,  md >= 2, 2)
  expV <- replace(expV,  hh >= 2, 3)
  expV <- replace(expV,  hr >= 3, 4)

  return(expV)
}
