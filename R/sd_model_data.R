#' @title Calculate Standard Deviation on MOM6 Data
#' @description
#'
#' @param raw spatRaster from \code{pull_mom6_hindcast} or \code{pull_mom6_forecast}
#' @param spatial.temporal TRUE/FALSE to determine how standard deviations are calculated If TRUE, standard deviation calculations are spatially and temporally explicit. If FALSE, the calculation occurs across space and time
#' 
#' @return If \code{spatial.temporal} is TRUE, a spatRaster is returned. If FALSE, a single value is returned. 

sd_model_data <- function(raw, spatial.temporal){
  
  if(spatial.temporal){
    ## create monthly averages across space and time (months - should make this customizable at some point)
    sds <- NULL
    for (m in 1:12) {
      mn <- seq(m, terra::nlyr(v), by = 12) #grab all month xs from timeseries by creating a sequence
      MNS <- v[[mn]] #subset raster brick
      mm <- terra::stdev(MNS, na.rm = T) #average
      sds <- abind::abind(as.array(mm), sds, along = 3) #make array and bind together
      #print(m)
    } #end m
    #convert to rasterBrick
    
    sds <- terra::rast(sds)
    terra::ext(sds) <- terra::ext(v) #make the extent the same as v
    terra::crs(sds) <- terra::crs(v)
    names(sds) <- 1:12
  } else {
    #calculate global average across layers and space
    sds <- terra::global(v, 'sd', na.rm = T)
  } #end if spatial.temporal
  
  return(sds)
} #end function
