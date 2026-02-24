#' @title Make Map of Total Species-Specific Exposure
#'
#' @description Uses the logic rule from CVA1.0 to combine variable-specific exposures to create a total exposure map
#'
#'
#' @param mapExp rasterStack output from \code{makeExposureMap}
#' @param countAll TRUE/FALSE to use \code{weights} and \code{wThreshold} to subset variables to only important variables
#' @param weights output from \code{combineWeights} - a vector of variable weights in ensemble SDM
#' @param wThreshold numeric value used to subset weights, variables with weights less to or equal to this value will be excluded from total exposure calculation
#'
#' @return A single raster representing total exposure across space


combineMaps <- function(mapExp, countAll, weights, wThreshold){
  if(countAll){
    mapSub <- mapExp
  } else {
    wi <- which(weights >= wThreshold)
    mapSub <- raster::subset(mapExp,wi)
  }

  #multiply mean exposure maps by variable weights to scale, and count how many layers have each rank within each cell
  hr <- raster::calc(mapSub, function(x){length(x[x >= 3.5])})
  hh <- raster::calc(mapSub, function(x){length(x[x >= 3])})
  md <- raster::calc(mapSub, function(x){length(x[x >= 2.5])})

  hr[is.na(raster::subset(mapSub,1))] <- NA
  hh[is.na(raster::subset(mapSub,1))] <- NA
  md[is.na(raster::subset(mapSub,1))] <- NA

  #apply logic rule from Hare et al 2015
  expL <- raster::subset(mapSub,1)
  expL[] <- 1
  expL <- replace(expL,  md >= 2, 2)
  expL <- replace(expL,  hh >= 2, 3)
  expL <- replace(expL,  hr >= 3, 4)
  expL[is.na(raster::subset(mapSub,1))] <- NA

  return(expL)
}
