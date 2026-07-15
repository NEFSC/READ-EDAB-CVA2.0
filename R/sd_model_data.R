#' @title Calculate Standard Deviation on MOM6 Data
#' @description Calculate monthly or global standard deviations from raw model data across the entire provided timeseries. This is built specifically for MOM6 output, but would work on any spatRaster of gridded data.
#'
#' @param raw spatRaster from \code{pull_mom6_hindcast} or \code{pull_mom6_forecast}
#' @param spatial_temporal TRUE/FALSE to determine how standard deviations are calculated If TRUE, standard deviation calculations are spatially and temporally explicit. If FALSE, the calculation occurs across space and time
#' 
#' @return If \code{spatial.temporal} is TRUE, a spatRaster is returned. If FALSE, a single value is returned. 

sd_model_data <- function(raw, spatial_temporal){
  
  if(spatial_temporal){
    ## create monthly averages across space and time (months - should make this customizable at some point)
    sds <- NULL
    for (m in 1:12) {
      mn <- seq(m, terra::nlyr(raw), by = 12) #grab all month xs from timeseries by creating a sequence
      MNS <- raw[[mn]] #subset raster brick
      mm <- terra::stdev(MNS, na.rm = T) #average
      sds <- abind::abind(as.array(mm), sds, along = 3) #make array and bind together
      #print(m)
    } #end m
    #convert to rasterBrick
    
    sds <- terra::rast(sds)
    terra::ext(sds) <- terra::ext(raw) #make the extent the same as v
    terra::crs(sds) <- terra::crs(raw)
    names(sds) <- 1:12
  } else {
    #calculate global average across layers and space
    gMean <- mean(terra::global(raw, 'mean', na.rm = T)$mean) #first calculate global mean
    
    squared_deviations <- (raw - gMean)^2 #then calculate squared differences 
    
    #Get the sum of all these squared deviations
    sum_of_squares_df <- terra::global(squared_deviations, "sum", na.rm = TRUE)
    total_sum_of_squares <- sum(sum_of_squares_df$sum)
    
    #get n-1 
    N_minus_1 <- sum(terra::global(raw, "notNA")$notNA) - 1
    
    #Take the square root of the variance to get your single final SD value
    sds <- sqrt(total_sum_of_squares / N_minus_1)
  } #end if spatial.temporal
  
  return(sds)
} #end function
