#' @title Make Maps or Timeseries of Variable-Specific Exposure
#'
#' @description Combines raw variable exposures and SDM results, using the SDM results as weights to average the exposure across time to produce maps of species-specific exposure for each variable or across space to produce timeseries of species-specific exposure for each variable
#'
#' @param type designates desired output, must equal 'map' or 'timeseries'
#' @param ranked_exposure output from \code{rank_exposure}
#' @param sdm_raster rasterStack with monthly averages of SDM results to use as weights for weighted average
#'
#' @return If \code{type == 'map'}, a rasterStack with the number of layers equal to the number of variables supplied, containing the weighted average exposure, weighted by SDM results, for each variable. If \code{type == 'timeseries'}, a matrix with the number of rows equal to the number of variables supplied, and 12 columns (1 for each month) containing the weighted average exposure, weighted by SDM results, for each variable.


make_variable_exposure <- function(type, ranked_exposure, sdm_raster){

  if(type == 'map'){
    mapExp <- vector(mode = 'list', length = length(ranked_exposure))
    for(x in 1:length(ranked_exposure)){
      s <- ranked_exposure[[x]]

      mapExp[[x]] <- weighted.mean(raster::stack(s), w = sdm_raster, na.rm = T)
      # print(x)
    }
    names(mapExp) <- names(ranked_exposure)
    return(raster::stack(mapExp))
  }

  if(type == 'timeseries'){
    matExp <- matrix(nrow = length(ranked_exposure), ncol = 12)
    for(x in 1:length(ranked_exposure)){
      s <- ranked_exposure[[x]]

      meanExp <- vector(length = 12) #vector
      for(m in 1:12){
        r <- raster::subset(s, m)
        h <- raster::subset(sdm_raster, m)
        r[is.na(h)] <- NA

        rDF <- raster::rasterToPoints(r)
        hDF <- raster::rasterToPoints(h)

        rh <- merge(hDF, rDF, by = c('x', 'y'), all.x = T)
        meanExp[m] <- weighted.mean(rh[,4], w = rh[,3], na.rm = T)
      }
      matExp[x,] <- meanExp
      #print(x)
    }
    rownames(matExp) <- names(ranked_exposure)
    colnames(matExp) <- month.abb

    return(matExp)
  }
}
