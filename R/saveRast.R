#' @title saveRast
#' @description A wrapper function for the \code{create_rast} function. This function adds log and skip functionality to help with batch runs and running \code{create_rast} in parallel.

#' @param csvName character string indicating which csv to use to create raster without extension (i.e. 'example' NOT 'example.csv')
#' @param isObs TRUE/FALSE indicating whether or not the data is observer or similar fisheries-dependent data. Will force check to determine if species is observed at least 30 times throughout timeseries before creating raster
#' @param grid static link to a ncdcf object with the variables lon, lat, time - can be link to remote data - must be able to be read with nc_open
#' @param spp Species name to add to log files and save data to correct directory (see vignette for recommended directory set up)
#' @param sppNames a vector containing all possible names for the target species. Must have a length >= 1
#' @param skip TRUE/FALSE indicating whether to skip creating the raster file if file already exists

#' @return \code{saveRast} returns the range of the rasterBrick returned by \code{create_rast}. This should be equal to 0 2 if species is present in dataset, and 0 1 if species is not caught in dataset. This function will also save the resulting rasterBrick as a netcdf file in the species' input_rasters folder

saveRast <- function(csvName, isObs, grid, spp, sppNames, skip){

  data <- read.csv(paste0('./Data/csvs/standardized/', csvName, '.csv'))

  #open log file
  sink(file = file.path(getwd(), 'logs', paste0(csvName, '.log')), append = T)
  #sink(file = file.path(getwd(), 'logs', paste0(csvName, '.log')), append = T, type = 'message')

  # Ensure the sinks are closed when the function exits, regardless of how it exits.
  on.exit({
    #sink(type = "message")
    sink()
  })

  print(paste(csvName , spp, sep = '-'))
  print(Sys.time())

  if(skip){
    if(file.exists(paste0(file.path(getwd(), spp, 'input_rasters'), '/', csvName, '.nc'))){
      print('file exists and skip == T, so skipping this file!')
      return(NA)
    } else {

      nms <- strsplit(sppNames, split = ',')[[1]]
      rast <- create_rast(data = data, isObs = isObs, grid = grid, targetVec = nms)

      if(is.null(rast)){
        print('rast is NULL - minimum conditions not met')
      } else {
        print(range(rast[]))
       raster::writeRaster(rast, filename = paste0(file.path(getwd(),spp, 'input_rasters'), '/', csvName, '.nc'), bylayer = F,overwrite = T)
      }

      return(range(rast[]))

    }#end skip && if file is present
  } else { #if skip = F, just run it without checking
    nms <- strsplit(sppNames, split = ',')[[1]]
    rast <- create_rast(data = data, isObs = isObs, grid = grid, targetVec = nms)

    if(is.null(rast)){
      print('rast is NULL - minimum conditions not met')
    } else {
      print(range(rast[]))
      raster::writeRaster(rast, filename = paste0(file.path(getwd(),spp, 'input_rasters'), '/', csvName, '.nc'), bylayer = F,overwrite = T)
    }

    return(range(rast[]))
  }

  sink() #close log file
} #end function
