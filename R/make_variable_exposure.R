#' @title Make Maps or Timeseries of Variable-Specific Exposure
#'
#' @description Combines raw variable exposures and SDM results, using the SDM results as weights to average the exposure across time to produce maps of species-specific exposure for each variable or across space to produce timeseries of species-specific exposure for each variable
#'
#' @param type designates desired output, must equal 'map' or 'timeseries'
#' @param ranked_exposure a list of spatRasters of ranked exposure data. Output from \code{rank_exposure}
#' @param sdm_raster spatRaster with monthly averages of SDM results to use as weights for weighted average
#'
#' @return If \code{type == 'map'}, a spatRaster with the number of layers equal to the number of variables supplied, containing the weighted average exposure, weighted by SDM results, for each variable. If \code{type == 'timeseries'}, a matrix with the number of rows equal to the number of variables supplied, and 12 columns (1 for each month) containing the weighted average exposure, weighted by SDM results, for each variable.
#'
#'@export

make_variable_exposure <- function(type, ranked_exposure, sdm_raster, stock_polys = NULL) {
  if (type == 'map') {
    mapExp <- vector(mode = 'list', length = length(ranked_exposure))
    for (x in 1:length(ranked_exposure)) {
      #shouldn't need to but might need to project to exact same extent if they are off slightly
      mapExp[[x]] <- terra::weighted.mean(
        ranked_exposure[[x]],
        w = sdm_raster,
        na.rm = T
      )
      # print(x)
    }
    names(mapExp) <- names(ranked_exposure)
    return(terra::rast(mapExp))
  } #end if map 

  if (type == 'timeseries') {
    #global
    matExp <- matrix(nrow = length(ranked_exposure), ncol = 12)
    for (x in 1:length(ranked_exposure)) {
      s <- ranked_exposure[[x]]

      meanExp <- vector(length = 12) #vector
      for (m in 1:12) {
        r <- s[[m]]
        h <- sdm_raster[[m]]

        meanExp[m] <- terra::global(r*h, 'sum', na.rm = T) / terra::global(h, 'sum', na.rm = T) #calculate global weighted average by hand 
      }
      matExp[x, ] <- unlist(meanExp)
      #print(x)
    }
    rownames(matExp) <- names(ranked_exposure)
    colnames(matExp) <- month.abb
    
    if(!is.null(stock_polys)){ #if stock polygons are provided, 
      stock_mats <- vector(mode = 'list', length = length(stock_polys)) #confirm length is the correct call here
      for(p in 1:length(stock_polys)){
        pExp <- matrix(nrow = length(ranked_exposure), ncol = 12)
        for (x in 1:length(ranked_exposure)) {
          s <- ranked_exposure[[x]]
          
          meanExp <- vector(length = 12) #vector
          for (m in 1:12) {
            r <- s[[m]]
            h <- sdm_raster[[m]]
            
            r_stack <- c(r, h)
            names(r_stack) <- c('exp', 'sdm')
            
            vals <- terra::extract(r_stack, stock_polys[p]) #again check indexing here 
            meanExp[m] <- sum(vals$exp * vals$sdm, na.rm = T) / sum(vals$sdm, na.rm = T) #calculate weighted average 
          }
          pExp[x, ] <- meanExp
          #print(x)
        } #end x 
        rownames(pExp) <- names(ranked_exposure)
        colnames(pExp) <- month.abb
        stock_mats[[p]] <- pExp #add matrix to list 
      }
      matExp <- c(list(matExp), stock_mats) #make list of all matrices
      names(matExp) <- c('global', stock_polys$stock_area) #name each matrix 
    }

    return(matExp)
  } #end if timeseries
}
