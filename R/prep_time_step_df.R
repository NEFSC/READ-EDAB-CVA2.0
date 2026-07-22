#' @title Build data.frame for single timestep
#' @description Helper function used in \code{make_sdm_predictions} to combine all environmental rasters for a single timestep and necessary static variables to predict model into a single data.frame
#'
#' @param x a numeric corresponding to the layer to pull. Generated within \code{make_sdm_predictions}.
#' @param rasts list of environmental spatRasters to use in predictions. Each spatRaster should have the name number of layers, corresponding to different timestamps
#' @param static_variables spatRaster containing the static variables used in model. 

#' @return a data.frame containing all of the data in \code{rasts} and \code{static_variables} for a single timestep. 
#'
#'@export

prep_time_step_df <- function(x, rasts, static_variables){
  
  # 0. Get the month/year name of layer x
  lyr_name <- names(rasts[[1]])[x]  # Cleaner: no nested [[1]] needed
  mm_yr <- strsplit(lyr_name, split = "[.]")[[1]]
  mm <- as.numeric(gsub("X", "", mm_yr[1]))
  yr <- as.numeric(mm_yr[2])
  
  # 1. Isolate just the dynamic layers for this time step
  dynamic_lyrs <- lapply(rasts, function(r) r[[x]])
  dynamic_stack <- terra::rast(dynamic_lyrs)
  
  # 2. Extract the data frame from the environmental variables ONLY 
  # This forces terra to stream the data from the NetCDF without coordinate friction
  env_df <- terra::as.data.frame(dynamic_stack, xy = TRUE, na.rm = FALSE)
  
  # 3. Extract the static variables as a data frame
  static_df <- terra::as.data.frame(static_variables, xy = FALSE, na.rm = FALSE)
  
  # 4. Bind them together side-by-side (Since grids match, rows align perfectly)
  sr_df <- cbind(
    x = env_df$x,
    y = env_df$y,
    month = mm,  # Use the scalar month parsed earlier
    year = yr,   # Use the scalar year parsed earlier
    env_df[, -c(1,2), drop = FALSE], # Drop the x/y columns from env_df to avoid duplication
    static_df
  )
  
  # 5. Clean up any rows where actual environmental data is missing
  sr_df <- sr_df[stats::complete.cases(sr_df), ]
  
  return(sr_df)
}