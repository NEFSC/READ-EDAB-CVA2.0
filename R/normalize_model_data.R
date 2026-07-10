#' @title Normalize Model Data
#' @description
#' Normalize model data across the entire provided timeseries using a z-score
#'
#' @param raw Output spatRaster from \code{pull_mom6_hindcast} or \code{pull_mom6_forecast}. Alternatively, a spatRaster of data to be normalized
#' @param avg,sd Output spatRaster from \code{avg_model_data} and \code{sd_model_data}, respectively. Alternatively, averages and standard deviations of desired gridded datasets. 
#' @param spatial.temporal TRUE/FALSE to determine how normalization should occur. If TRUE, calculations should be spatially and temporally explicit. \code{avg} and \code{sd} should be spatRaster objects with averages and standard deviations of data in space and time. If FALSE, the calculation occurs across space and time; \code{avg} and \code{sd} should be single values representing the overall average and standard deviation for the given variable.
#'
#' @return a spatRaster of normalized environmental data


normalize_model_data <- function(raw, avg, sd, spatial_temporal){
  
  if(spatial_temporal){
    ##normalize data to monthly spatially and temporally explicit averages and sds 
    mth <- rep(1:12, times = raster::nlayers(v)/12) #creates repeating list of 1:12 for each year
    #normalize data
    norm <- NULL
    for(m in 1:raster::nlayers(raw)){
      #subset both rasterbricks
      subX <- raster::subset(raw,m)
      subA <- raster::subset(avg,mth[m])
      subS <- raster::subset(sd,mth[m])
      
      nm <- (subX - subA) / subS
      norm <- abind::abind(as.array(nm), norm, along = 3)
    } #end m
    norm <- terra::rast(norm)
    terra::ext(norm) <- terra::ext(raw) #make the extent the same as v
    terra::crs(norm) <- terra::crs(raw)
  } else {
    #use averages/standard deviations calculated across space and time (single values)
    norm <- terra::app(raw, fun = function(x){(x - avg) / sd})
  }
  
  return(norm)
  
} #end function
