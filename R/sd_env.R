#' @title Calculate Standard Deviation on MOM6 Data
#' @description
#' Calculate monthly standard deviation from raw MOM6 data across the entire provided timeseries
#'
#' @param rawList List of output rasters from \code{pull_hind} or \code{pull_forecast}
#'
#' @return a list whose length is equal to the number of variables supplied, where each item in the list is a rasterStack of data associated with that variable

sd_env <- function(rawList){

  sdList <- vector(mode = 'list', length = length(rawList))

  for(x in 1:length(rawList)){
    v <- rawList[[x]]

    ## create monthly averages
    sds <- NULL
    for(m in 1:12){
      mn <- seq(m, raster::nlayers(v), by = 12) #grab all month xs from timeseries by creating a sequence
      MNS <- raster::subset(v, mn) #subset raster brick
      mm <- raster::calc(MNS, fun = sd, na.rm = T) #average
      sds <- abind::abind(as.array(mm), sds, along = 3) #make array and bind together
      #print(m)
    } #end m
    #convert to rasterBrick

    sds <- raster::brick(sds)
    raster::extent(sds) <- raster::extent(v) #make the extent the same as v
    raster::crs(sds) <- raster::crs(v)
    names(sds) <- 1:12
    sdList[[x]] <- sds

  } #end x

  #names(sdList) <- names(rastList)
  return(sdList)
} #end function
