#' @title Average MOM6 Data
#' @description
#' Calculate monthly average from raw MOM6 data across the entire provided timeseries
#'
#'
#' @param rawList List of output rasters from \code{pull_hind} or \code{pull_forecast}
#'
#' @return a list whose length is equal to the number of variables supplied, where each item in the list is a rasterStack of data associated with that variable


avg_env <- function(rawList){

  avgList <- vector(mode = 'list', length = length(rawList))

  for(x in 1:length(rawList)){
    v <- rawList[[x]]

    ## create monthly averages
    avgs <- NULL
    for(m in 1:12){
      mn <- seq(m, raster::nlayers(v), by = 12) #grab all month xs from timeseries by creating a sequence
      MNS <- raster::subset(v, mn) #subset raster brick
      mm <- raster::mean(MNS) #average
      avgs <- abind::abind(as.array(mm), avgs, along = 3) #make array and bind together
      #print(m)
    } #end m
    #convert to rasterBrick

    avgs <- raster::brick(avgs)
    raster::extent(avgs) <- raster::extent(v) #make the extent the same as v
    raster::crs(avgs) <- raster::crs(v)
    names(avgs) <- 1:12
    avgList[[x]] <- avgs

  } #end x

  names(avgList) <- names(rawList)
  return(avgList)
} #end function
