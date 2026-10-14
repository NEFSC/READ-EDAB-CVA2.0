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
#' @param threshold upper and lower bound limit to define "bad" data. Defaults to c(-9999, 9999).
#'
#' @return a spatRaster with a number of layers equal to the number of layers in \code{rasts} + number of layesr in \code{static_variables} + 4 (lon, lat, month, year)
#'
#'@export

prep_time_step_stack <- function(
  x,
  rasts,
  template_r,
  rlon,
  rlat,
  static_variables,
  month_col,
  year_col,
  xy_col,
  threshold = c(-9999, 9999)
) {
  # 1. Get the month/year name of layer x
  lyr_name <- names(rasts[[1]])[x] # Cleaner: no nested [[1]] needed
  mm_yr <- strsplit(lyr_name, split = "[.]")[[1]]
  mm <- as.numeric(gsub("X", "", mm_yr[1]))
  yr <- as.numeric(mm_yr[2])

  # 2. Create constant rasters natively in terra
  rMonth <- terra::rast(template_r, vals = mm)
  rYear <- terra::rast(template_r, vals = yr)

  # 3. Pull layer 'x' from every SpatRaster in your 'rasts' list
  dynamic_lyrs <- lapply(rasts, function(r) r[[x]])
  dynamic_stack <- terra::rast(dynamic_lyrs)
  
  # 4. Layer-specific Data Quality Check
  layer_ranges <- terra::minmax(dynamic_stack)
  
  # Evaluate conditions per layer
  is_inf <- is.infinite(layer_ranges["min", ]) | is.infinite(layer_ranges["max", ])
  is_extreme <- layer_ranges["max", ] > threshold[2] | layer_ranges["min", ] < threshold[1]
  
  # Identify indices of layers that fail the check (using which() safely handles NAs)
  bad_idx <- which(is_inf | is_extreme)
  
  if (length(bad_idx) > 0) {
    bad_names <- names(dynamic_stack)[bad_idx]
    message(paste("Corrupted data detected in variable(s):", paste(bad_names, collapse = ", "), 
                  "at timestep", lyr_name, "- converting only these layers to NA."))
    
    # Replace ONLY the corrupted layers with NA
    for (i in bad_idx) {
      dynamic_stack[[i]] <- terra::init(dynamic_stack[[i]], NA)
    }
  }
  
  # 4. Force the static variables to align perfectly with your dynamic data grid
  # This matches extents, projections, and resolutions exactly
  static_aligned <- terra::resample(
    static_variables,
    template_r,
    method = "bilinear"
  )

  # 5. Combine everything using terra::c()
  sr_terra <- c(
    rlon,
    rlat,
    rMonth,
    rYear,
    dynamic_stack,
    static_aligned
  )
  names(sr_terra)[1:4] <- c(xy_col, month_col, year_col)

  return(sr_terra)
}
