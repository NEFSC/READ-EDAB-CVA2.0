#' @title Convert Environmental Rasters to a Data.Frame
#' @description Make a data frame out of environmental data to be used for predictions. Helps speed up sdmTMB predictions.
#'
#' @param rasts list of rasterStacks corresponding to the environmental covariates used to build the models. The number of layers in each rasterStack should be the same and correspond to the length of the timeseries for the models to be predicted on
#' @param static_variables list of rasters containing the static variables used in model. Should be the same object used in \code{merge_spp_env}.
#' @param mask TRUE/FALSE indicating whether to mask off certain depths (i.e. waters deeper than 1000 m)
#' @param bathy_raster,bathy_max a raster of bathymetry with the same extent and resolution as \code{rasts} and the maximum depth you want included. For example, if you want to mask off waters deeper than 1000 m, \code{bathy_max} would be set to 1000. The value should be positive regardless of the sign of your bathymetry data.
#'
#' @return a data.frame containing all of the environmental data in \code{rasts} for all timesteps. It will be big depending on the length of your time series. This is meant to be used to help predict sdmTMB models.
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

raster_to_df <- function(
  rasts,
  static_variables,
  bathy_raster,
  bathy_max,
  mask
) {
  #make static vars (month/year) into rasters
  r <- raster::subset(rasts[[1]][[1]], 1)
  rlon <- rlat <- r #copy r to rlon and rlat rasters [1]][1]which will contain the longitude and latitude
  xy <- raster::xyFromCell(r, 1:length(r)) #matrix of longitudes (x) and latitudes(y)
  rlon[] <- xy[, 1] #raster of longitudes
  rlat[] <- xy[, 2] #raster of latitides
  rMonth <- rYear <- r

  allDF <- NULL
  for (x in 1:raster::nlayers(rasts[[1]][[1]])) {
    #all the rasters in rasts have the same number of layers so it doesn't matter which one we call

    #may need to change if names aren't always going to be month.year
    mm.year <- strsplit(names(rasts[[1]][[1]])[x], split = '[.]') #all the rasterbricks in normVars also have the same names so again, doesn't matter which one we call

    mm <- as.numeric(gsub('X', '', mm.year[[1]][1]))
    rMonth[] <- mm
    yr <- as.numeric(mm.year[[1]][2])
    rYear[] <- yr

    nStack <- vector(mode = 'list', length = length(rasts))
    for (n in 1:length(rasts)) {
      nStack[[n]] <- raster::subset(rasts[[n]][[1]], x)
    }
    nStack <- c(nStack, static_variables)
    nStack <- raster::stack(nStack)
    names(nStack) <- c(names(rasts), names(static_variables))
    raster::crs(nStack) <- raster::crs(rasts[[1]][[1]])
    raster::extent(nStack) <- raster::extent(rasts[[1]][[1]])

    sr <- raster::stack(rlon, rlat, rMonth, rYear, nStack)
    names(sr)[1:4] <- c("x", "y", "month", "year")

    if (mask) {
      #if mask == T
      #mask off waters deeper than 1000 m
      #i <- which(names(sr) == bathy_nm)
      sr <- replace(sr, abs(bathy_raster) > bathy_max, NA) #replace values with an absolute value greater than bathy_max with NA
    }

    ##convert rasterStack to dataframe to play well with model
    srDF <- as.data.frame(raster::rasterToPoints(sr)[, -c(1:2)])
    srDF <- srDF[complete.cases(srDF), ]
    allDF <- rbind(allDF, srDF)
  }

  return(allDF)
} #end function
