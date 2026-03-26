#' @title Normalize Model Data
#' @description
#' Normalize model data across the entire provided timeseries using a z-score
#'
#' @param raw_list List of output rasters from \code{pull_hind} or \code{pull_forecast}. Alternatively, a list of rasters containing gridded data to be normalized, with each entry in the list containing a rasterStack of a timeseries of environmental data to be normalized.
#' @param avg_list,sd_list List of output rasters from \code{avg_model_data} and \code{sd_env}, respectively. Alternatively, averages and standard deviations of desired gridded datasets. Should be a list of rasterStacks, with each entry in the list containing the rasterStack for a different variable.
#' @param short_names vector of simplified variable names to help name resulting raster files.
#'
#' @return a list whose length is equal to the length of \code{raw_list}, where each item in the list is a rasterStack of normalized data associated with that variable


normalize_model_data <- function(raw_list, avg_list, sd_list, short_names){
  normList <- vector(mode = 'list', length = length(raw_list))

  for(x in 1:length(raw_list)){
    v <- raw_list[[x]]
    avgs <- avg_list[[x]]
    sds <- sd_list[[x]]

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

  names(normList) <- short_names #set names of lists

  return(normList)

} #end function
