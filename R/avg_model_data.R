#' @title Average Model Data
#' @description
#' Calculate monthly or global averages from raw model data across the entire provided timeseries. This is built specifically for MOM6 output, but would work on any spatRaster of gridded data.
#'
#' @param raw spatRaster from \code{pull_mom6_hindcast} or \code{pull_mom6_forecast}
#' @param spatial_temporal TRUE/FALSE to determine averaging method. If TRUE, averages are spatially and temporally explicit. If FALSE, averaging occurs across space and time
#'
#' @return If \code{spatial.temporal} is TRUE, a spatRaster is returned. If FALSE, a single value is returned. 
#'
#'@export

avg_model_data <- function(raw, spatial_temporal) {
  
  if(spatial_temporal){
    ## create monthly averages across space and time (months - should make this customizable at some point)
   avgs <- terra::tapp(raw, rep(1:12, times = terra::nlyr(raw)/12), fun = 'mean')
    names(avgs) <- month.abb
  } else {
    #calculate global average across layers and space
    avgs <- mean(terra::global(raw, 'mean', na.rm = T)$mean, na.rm = T)
  } #end if spatial.temporal
  
  return(avgs)
  
} #end function
