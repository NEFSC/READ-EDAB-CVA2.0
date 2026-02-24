#' @title Calculate Raw Exposure from Environmental Data
#' @description
#' Calculates the raw variables exposure as the ratio of the difference between future and present conditions, divided by present standard deviations (a z-score).
#'
#'
#' @param presentRasts,futureRasts lists of rasterStacks, representing the present and future timeseries to calculate exposure with, with each item in the list corresponding to an environmental variable. Both lists need to have the same length (aka number of variables) but the lengths of the timeseries (equal to the number of layers in rasterStacks) can differ. Function assumes that the rasters are monthly timeseries and will average monthly.
#' @return a list of rasterStacks with each layer representing monthly raw exposure. The length of the list is equal to the length of the lists supplied as \code{presentRasts} and \code{futureRasts}

calcExposure <- function(presentRasts, futureRasts){
  EXP <- vector(mode = 'list', length = length(futureRasts))
  for(v in 1:length(futureRasts)){
    ex <- vector(mode = 'list', length = 12)
    for(x in 1:12){
      #take mean of 'present' and 'future'
      mPres <- raster::calc(raster::subset(presentRasts[[v]][[1]], seq(x, raster::nlayers(presentRasts[[v]][[1]]), by = 12)), mean)
      mFut <- raster::calc(raster::subset(futureRasts[[v]][[1]], seq(x, raster::nlayers(futureRasts[[v]][[1]]), by = 12)), mean)

      #difference in means
      mDiff <- mFut - mPres

      ##calculate SD
      sdPres <- raster::calc(raster::subset(presentRasts[[v]][[1]], seq(x, raster::nlayers(presentRasts[[v]][[1]]), by = 12)), sd)
      sdFut <- raster::calc(raster::subset(futureRasts[[v]][[1]], seq(x, raster::nlayers(futureRasts[[v]][[1]]), by = 12)), sd)

      #calculate exposure
      ex[[x]] <- raster::stack(mDiff) / raster::stack(sdPres)
    }

    #stack
    EXP[[v]] <- raster::stack(ex)
    print(v)
  }
  names(EXP) <- names(futureRasts)
  return(EXP)
}
