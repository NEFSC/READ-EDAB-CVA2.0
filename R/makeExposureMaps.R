#' @title Make Maps of Species-Specific Exposure for each variable
#' @description Combines raw variable exposures and SDM results, using the SDM results as weights to average the exposure across time to produce maps of species-specific exposure for each variable
#' @param rankExp output from \code{rankExposure}
#' @param sdmRast rasterStack with monthly averages of SDM results to use as weights for weighted average
#' @return a rasterStack with the number of layers equal to the number of variables supplied, containing the weighted average exposure, weighted by SDM results, for each variable.


makeExposureMaps <- function(rankExp, sdmRast){
  mapExp <- vector(mode = 'list', length = length(rankExp))
  for(x in 1:length(rankExp)){
    s <- rankExp[[x]]

    mapExp[[x]] <- weighted.mean(raster::stack(s), w = sdmRast, na.rm = T)
    # print(x)
  }
  names(mapExp) <- names(rankExp)
  return(raster::stack(mapExp))
}
