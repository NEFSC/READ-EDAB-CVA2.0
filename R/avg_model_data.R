#' @title Average Model Data
#' @description
#' Calculate monthly average from raw model data across the entire provided timeseries. This is built specificially for MOM6 output, but would work on any rasterStack of gridded data.
#'
#' @param raw_list List of output rasters from \code{pull_hind} or \code{pull_forecast}
#'
#' @return a list whose length is equal to the number of variables supplied, where each item in the list is a rasterStack of data associated with that variable


avg_model_data <- function(raw_list){

  avgList <- vector(mode = 'list', length = length(raw_list))

  for(x in 1:length(raw_list)){
    v <- raw_list[[x]]

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

  names(avgList) <- names(raw_list)
  return(avgList)
} #end function
