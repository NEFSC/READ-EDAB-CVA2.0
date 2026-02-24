#' @title Calculate Change in Distribution Metrics
#' @description Calculate two metrics to help describe changes in species distributions - Center of Gravity and Area of Preferred Habitat
#'
#'
#' @param abund model predictions, output from \code{make_predictions} as a rasterStack. Layers need names that correspond to the timestamps represented.
#' @param area.threshold a value between 0 and 1, defaults to 0.75. The total area of grid cells with the probability of occurance of a species equal to or above this threshold will be calculated.
#' @param cell.area a number representing the area of a single model grid cell on which the model is predicted.
#'
#' @return a data.frame with three columns: 1) timestamp, equal to the names of the layers in \code{abund}; 2) COG; 3) Area of probabilities greater than the \code{area.threshold}

change_in_distribution <- function(abund, area.threshold = 0.75, cell.area = 8*8){
  #calculates both center of gravity and area for each time step in abund
  #returns a data.frame

  #center of gravity - modeled after DiSMAP
  mDF <- as.data.frame(raster::rasterToPoints(abund))
  cog <- vector(length = ncol(mDF))
  for(x in 3:ncol(mDF)){
    wt <- mDF$y * mDF[,x]
    cog[x] <- sum(wt, na.rm = T) / sum(mDF[,x], na.rm = T)
  }
  cog <- cog[-c(1:2)]

  #area at each timestamp
  val <- raster::values(abund) #extract all values
  area <- apply(val, 2, FUN = function(x){length(which(x >= area.threshold))}) * cell.area

  #make data.frame
  df <- as.data.frame(cbind(names(abund), as.numeric(cog), as.numeric(area)))
  names(df) <- c('timestamp', 'cog', 'area')
  df$cog <- as.numeric(df$cog)
  df$area <- as.numeric(df$area)

  return(df)
}






