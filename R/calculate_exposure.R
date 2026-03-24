#' @title Calculate Raw Exposure from Environmental Data
#' @description
#' Calculates the raw variables exposure as the ratio of the difference between future and present conditions, divided by present standard deviations (a z-score).
#'
#' @param present_rasters,future_rasters lists of rasterStacks, representing the present and future timeseries to calculate exposure with, with each item in the list corresponding to an environmental variable. Both lists need to have the same length (aka number of variables) but the lengths of the timeseries (equal to the number of layers in rasterStacks) can differ. Function assumes that the rasters are monthly timeseries and will average monthly.
#' @return a list of rasterStacks with each layer representing monthly raw exposure. The length of the list is equal to the length of the lists supplied as \code{present_rasters} and \code{future_rasters}

calculate_exposure <- function(present_rasters, future_rasters){
  EXP <- vector(mode = 'list', length = length(future_rasters))
  for(v in 1:length(future_rasters)){
    ex <- vector(mode = 'list', length = 12)
    for(x in 1:12){
      #take mean of 'present' and 'future'
      mPres <- raster::calc(raster::subset(present_rasters[[v]][[1]], seq(x, raster::nlayers(present_rasters[[v]][[1]]), by = 12)), mean)
      mFut <- raster::calc(raster::subset(future_rasters[[v]][[1]], seq(x, raster::nlayers(future_rasters[[v]][[1]]), by = 12)), mean)

      #difference in means
      mDiff <- mFut - mPres

      ##calculate SD
      sdPres <- raster::calc(raster::subset(present_rasters[[v]][[1]], seq(x, raster::nlayers(present_rasters[[v]][[1]]), by = 12)), sd)
      sdFut <- raster::calc(raster::subset(future_rasters[[v]][[1]], seq(x, raster::nlayers(future_rasters[[v]][[1]]), by = 12)), sd)

      #calculate exposure
      ex[[x]] <- raster::stack(mDiff) / raster::stack(sdPres)
    }

    #stack
    EXP[[v]] <- raster::stack(ex)
    print(v)
  }
  names(EXP) <- names(future_rasters)
  return(EXP)
}
