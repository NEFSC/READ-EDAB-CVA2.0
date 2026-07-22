#' @title Build spatRaster for single timestep
#' @description Helper function used in \code{make_sdm_predictions} to combine all environmental rasters for a single timestep and necessary static variables to predict model
#'
#' @param x a numeric corresponding to the layer to pull. Generated within \code{make_sdm_predictions}.
#' @param rasts list of environmental spatRasters to use in predictions. Each spatRaster should have the name number of layers, corresponding to different timestamps
#' @param template_r a single-layer spatRaster to use as a template for other raster. Used to make sure \code{static_variables} have the correct resolution and projection, as well as make month/year rasters. Generated within \code{make_sdm_predictions}.
#' @param rlon,rlat rasters of longitude and latitude respectively. Generated within \code{make_sdm_predictions}.
#' @param static_variables spatRaster containing the static variables used in model. 
#' @param xy_col a vector with a length of 2 indicating the longitude and latitude column names
#' @param month_col,year_col column names for month and year columns respectively. 
#'
#' @return a spatRaster with a number of layers equal to the number of layers in \code{rasts} + number of layesr in \code{static_variables} + 4 (lon, lat, month, year)
#'
#'@export

prep_time_step_stack <- function(x, rasts, template_r, rlon, rlat, static_variables, month_col, year_col, xy_col) {
  # 1. Get the month/year name of layer x
  lyr_name <- names(rasts[[1]])[x]  # Cleaner: no nested [[1]] needed
  mm_yr <- strsplit(lyr_name, split = "[.]")[[1]]
  mm <- as.numeric(gsub("X", "", mm_yr[1]))
  yr <- as.numeric(mm_yr[2])
  
  # 2. Create constant rasters natively in terra
  rMonth <- terra::rast(template_r, vals = mm)
  rYear  <- terra::rast(template_r, vals = yr)
  
  # 3. Pull layer 'x' from every SpatRaster in your 'rasts' list
  # This dynamically grabs the same month across all your predictor variables
  dynamic_lyrs <- lapply(rasts, function(r) r[[x]]) 
  
  # 4. Force the static variables to align perfectly with your dynamic data grid
  # This matches extents, projections, and resolutions exactly
  static_aligned <- terra::resample(static_variables, template_r, method = "bilinear")
  
  # 5. Combine everything using terra::c()
  sr_terra <- c(rlon, rlat, rMonth, rYear, terra::rast(dynamic_lyrs), static_aligned)
  names(sr_terra)[1:4] <- c(xy_col, month_col, year_col)
  
  return(sr_terra)
}