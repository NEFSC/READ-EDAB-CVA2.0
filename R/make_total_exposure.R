#' @title Make Map or Timeseries of Total Species-Specific Exposure
#'
#' @description Uses the logic rule from CVA1.0 to combine variable-specific exposures to create a total exposure map or timeseries
#'
#' @param type designates desired output, must equal 'map' or 'timeseries'
#' @param variable_exposure If \code{type == 'map'}, a rasterStack output from \code{make_variable_exposure(type == 'map')}. If \code{type == 'timeseries'}, the matrix output from \code{make_variable_exposure(type == 'timeseries')}.
#' @param count_all TRUE/FALSE to use \code{weights} and \code{wThreshold} to subset variables to only important variables
#' @param weights output from \code{combine_weights} - a vector of variable weights in ensemble SDM
#' @param weight_threshold numeric value used to subset weights, variables with weights less to or equal to this value will be excluded from total exposure calculation
#'
#' @return If \code{type == 'map'}, the output is a raster representing total exposure across space. If \code{type == 'timeseries'}, the output is a vector representing total exposure across time.


make_total_exposure <- function(type, variable_exposure, count_all, weights, weight_threshold){

  if(type == 'map'){
    if(count_all){
      mapSub <- variable_exposure
    } else {
      wi <- which(weights >= weight_threshold)
      mapSub <- raster::subset(variable_exposure,wi)
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

  if(type == 'timeseries'){
    if(count_all){ #if counting all included factors and not taking weight into account
      matSub <- variable_exposure
    } else {
      wi <- which(weights > weight_threshold)
      matSub <- variable_exposure[wi,]
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
}
