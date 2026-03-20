#' @title Convert Predicted Values in Data.Frame to Rasters
#' @description Make rasters from predicted values in a data.frame. Helps speed up sdmTMB predictions.
#'
#' @param df the output from \code{makePredDF}
#' @param staticVars list of rasters containing the static variables used in model. Should be the same object used in \code{merge_spp_env}. Used in this instance to set resolution of rasters
#'
#' @return a rasterStack with the same number of layers as in \code{rasts} with probability of occurance predicted for each day with available data. Values will range from 0 to 1.
#'

makePredRasters <- function(df, staticVars){

  hsm <- vector(mode = 'list', length = length(unique(df$my)))
  for(x in 1:length(unique(df$my))){
    sub <- df[df$my == unique(df$my)[x],]
    sp::coordinates(sub) <- ~x + y
    sp::proj4string(sub) <- sp::CRS("+proj=longlat +datum=WGS84 +no_defs ")
    hsm[[x]] <- raster::rasterize(x = sub, y = staticVars[[1]], field = sub$est)
    print(x)
  }

  return(hsm)
}
