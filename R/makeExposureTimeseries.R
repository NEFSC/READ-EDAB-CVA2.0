#' @title Make Timeseries of Species-Specific Exposure for each variable
#' @description Combines raw variable exposures and SDM results, using the SDM results as weights to average the exposure across space to produce timeseries of species-specific exposure for each variable
#' @param rankExp output from \code{rankExposure}
#' @param sdmRast rasterStack with monthly averages of SDM results to use as weights for weighted average
#' @return A matrix with the number of rows equal to the number of variables supplied, and 12 columns (1 for each month) containing the weighted average exposure, weighted by SDM results, for each variable.

makeExposureTimeseries <- function(rankExp, sdmRast){
  matExp <- matrix(nrow = length(rankExp), ncol = 12)
  for(x in 1:length(rankExp)){
    s <- rankExp[[x]]

    meanExp <- vector(length = 12) #vector
    for(m in 1:12){
      r <- raster::subset(s, m)
      h <- raster::subset(sdmRast, m)
      r[is.na(h)] <- NA

      rDF <- raster::rasterToPoints(r)
      hDF <- raster::rasterToPoints(h)

      rh <- merge(hDF, rDF, by = c('x', 'y'), all.x = T)
      meanExp[m] <- weighted.mean(rh[,4], w = rh[,3], na.rm = T)
    }
    matExp[x,] <- meanExp
    #print(x)
  }
  rownames(matExp) <- names(rankExp)
  colnames(matExp) <- month.abb

  return(matExp)
}

