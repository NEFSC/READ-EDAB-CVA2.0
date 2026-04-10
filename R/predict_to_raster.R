#' @title Convert Predicted Values in Data.Frame to Rasters
#' @description Make rasters from predicted values in a data.frame. Helps speed up sdmTMB predictions.
#'
#' @param df the output from \code{raster_to_df}, with a column \code{est} containing the predicted values. See example.
#' @param static_variables list of rasters containing the static variables used in model. Should be the same object used in \code{merge_spp_env}. Used to set resolution of output rasters
#'
#' @return a rasterStack with the same number of layers as in \code{rasts} with probability of occurance predicted for each day with available data. Values will range from 0 to 1.
#'
#' @examples
#' \dontrun{
#' #create data frame of environmental data
#' allDF <- raster_to_df(rasts = rasts, static_variables = staticVars,
#' bathy_raster = bathyR, bathy_max = 1000, mask = T)
#'
#' #predict sdmTMB model for all data; not for each timestep as in \code{make_sdm_predictions}
#' #this generates the \code{est} column
#' pred <- predict(mod, newdata = allDF, type = 'response')
#'
#' #add appropriate month.year (my) column to create rasters
#' pred$my <- paste(pred$month, pred$year, sep = '.')
#' abund <- predict_to_raster(df = pred, staticData = staticVars) #make into rasters
#' }

predict_to_raster <- function(df, static_variables) {
  hsm <- vector(mode = 'list', length = length(unique(df$my)))
  for (x in 1:length(unique(df$my))) {
    sub <- df[df$my == unique(df$my)[x], ]
    sp::coordinates(sub) <- ~ x + y
    sp::proj4string(sub) <- sp::CRS("+proj=longlat +datum=WGS84 +no_defs ")
    hsm[[x]] <- raster::rasterize(
      x = sub,
      y = static_variables[[1]],
      field = sub$est
    )
    print(x)
  }

  return(hsm)
}
