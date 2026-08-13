#' @title Calculate Raw Exposure from Environmental Data
#' @description
#' Calculates the raw variables exposure as the ratio of the difference between future and present conditions, divided by present standard deviations (a z-score).
#'
#' @param present,future spatRasters, representing the present and future timeseries to calculate exposure with. Function assumes that the rasters are monthly timeseries and will average monthly.
#' @param spatial_temporal TRUE/FALSE to determine method for averaging present data ONLY. If TRUE, normalization is spatially and temporally explicit. If FALSE, normalization occurs using averages/standard deviations calculated across space and time. Future datasets are always averaged across time and space to ensure that exposure value is spatially explicit. 
#' @param mask_bathy TRUE/FALSE indicating whether or not to use bathymetry data as a mask for raw data
#' @param bathy a spatRaster of bathymetry data. Must be same extent as results from pull_mom6_hindcast/forecast
#' @param bathy_range a vector containing the minimum and maximum desired bathymetry values 
#' @return a spatRaster of monthly raw exposure 
#'
#'@export

calculate_raw_exposure <- function(present, future, spatial_temporal, mask_bathy, bathy, bathy_range) {
 
   if (mask_bathy) {
    #Unwrap the bathymetry raster inside the worker
    if (inherits(bathy, "PackedSpatRaster")) {
      bathy <- terra::unwrap(bathy)
    }
    #reproject because there are slight differences in resolution/extent for some reason, especially with the forecasts
    bathy_aligned <- terra::resample(bathy, present, method = "bilinear")
    
    present <- terra::ifel(bathy_aligned <= bathy_range[1] | bathy_aligned > bathy_range[2], NA, present)
    future <- terra::ifel(bathy_aligned <= bathy_range[1] | bathy_aligned > bathy_range[2], NA, future)
  }
  
  #calculate averages
  pAvg <- avg_model_data(present, spatial_temporal = spatial_temporal) #can change whether z-score is calculated relative to spatially/temporally explicit present averages or global averages 
  
  fAvg <- avg_model_data(future, spatial_temporal = TRUE) #always want future averages done across time and space (create monthly, spatially-explicit averages) so that exposure calculation is spatially/temporally explicit
  
  #calculate present SD 
  pSD <- sd_model_data(present, spatial_temporal)

  #calculate exposure
  EXP <- (fAvg - pAvg) / pSD
  
  return(EXP)
}
