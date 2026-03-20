#' @title Calculate Standard Deviation on MOM6 Data
#' @description
#' Calculate monthly standard deviation from raw MOM6 data across the entire provided timeseries
#'
#' @param rawList List of output rasters from \code{pull_hind} or \code{pull_forecast}
#' @param avgList,sdList List of output rasters from \code{avg_env} and \code{sd_env}, respectively
#' @param shortNames vector of simplified variable names to help name resulting raster files.
#'
#' @return a list whose length is equal to the number of variables supplied, where each item in the list is a rasterStack of data associated with that variable


norm_env <- function(rawList, avgList, sdList, shortNames){
  normList <- vector(mode = 'list', length = length(rawList))

  for(x in 1:length(rawList)){
    v <- rawList[[x]]
    avgs <- avgList[[x]]
    sds <- sdList[[x]]

    ##normalize data to monthly averages
    mth <- rep(1:12, times = raster::nlayers(v)/12) #creates repeating list of 1:12 for each year
    #normalize data
    norm <- NULL
    for(m in 1:raster::nlayers(v)){
      #subset both rasterbricks
      subX <- raster::subset(v,m)
      subA <- raster::subset(avgs,mth[m])
      subS <- raster::subset(sds,mth[m])

      nm <- (subX - subA) / subS
      norm <- abind::abind(as.array(nm), norm, along = 3)
    } #end m
    norm <- raster::brick(norm)
    raster::extent(norm) <- raster::extent(v) #make the extent the same as v
    raster::crs(norm) <- raster::crs(v)
    names(norm) <- names(v) #set names
    normList[[x]] <- norm #save
  } #end x

  names(normList) <- shortNames #set names of lists

  return(normList)

} #end function
