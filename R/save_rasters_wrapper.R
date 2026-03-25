#' @title Save and Build Fisheries Rasters
#' @description A wrapper function for the \code{build_fisheries_raster} function. This function adds log and skip functionality to help with batch runs and running \code{build_fisheries_raster} in parallel.

#' @param csv_name character string indicating which csv to use to create raster without extension (i.e. 'example' NOT 'example.csv')
#' @param is_obs TRUE/FALSE indicating whether or not the data is observer or similar fisheries-dependent data. Will force check to determine if species is observed at least 30 times throughout timeseries before creating raster
#' @param grid static link to a ncdcf object with the variables lon, lat, time - can be link to remote data - must be able to be read with nc_open
#' @param spp Species name to add to log files and save data to correct directory (see vignette for recommended directory set up)
#' @param spp_names a vector containing all possible names for the target species. Must have a length >= 1
#' @param skip TRUE/FALSE indicating whether to skip creating the raster file if file already exists

#' @return \code{saveRast} returns the range of the rasterBrick returned by \code{build_fisheries_raster}. This should be equal to 0 2 if species is present in dataset, and 0 1 if species is not caught in dataset. This function will also save the resulting rasterBrick as a netcdf file in the species' input_rasters folder

save_rasters_wrapper <- function(csv_name, is_obs, grid, spp, spp_names, skip){

  data <- read.csv(paste0('./Data/csvs/standardized/', csv_name, '.csv'))

  #open log file
  sink(file = file.path(getwd(), 'logs', paste0(csv_name, '.log')), append = T)
  #sink(file = file.path(getwd(), 'logs', paste0(csv_name, '.log')), append = T, type = 'message')

  # Ensure the sinks are closed when the function exits, regardless of how it exits.
  on.exit({
    #sink(type = "message")
    sink()
  })

  print(paste(csv_name , spp, sep = '-'))
  print(Sys.time())

  if(skip){
    if(file.exists(paste0(file.path(getwd(), spp, 'input_rasters'), '/', csv_name, '.nc'))){
      print('file exists and skip == T, so skipping this file!')
      return(NA)
    } else {

      nms <- strsplit(spp_names, split = ',')[[1]]
      rast <- build_fisheries_raster(data = data, is_obs = is_obs, grid = grid, all_names = nms)

      if(is.null(rast)){
        print('rast is NULL - minimum conditions not met')
      } else {
        print(range(rast[]))
       raster::writeRaster(rast, filename = paste0(file.path(getwd(),spp, 'input_rasters'), '/', csv_name, '.nc'), bylayer = F,overwrite = T)
      }

      return(range(rast[]))

    }#end skip && if file is present
  } else { #if skip = F, just run it without checking
    nms <- strsplit(spp_names, split = ',')[[1]]
    rast <- build_fisheries_raster(data = data, is_obs = is_obs, grid = grid, all_names = nms)

    if(is.null(rast)){
      print('rast is NULL - minimum conditions not met')
    } else {
      print(range(rast[]))
      raster::writeRaster(rast, filename = paste0(file.path(getwd(),spp, 'input_rasters'), '/', csv_name, '.nc'), bylayer = F,overwrite = T)
    }

    return(range(rast[]))
  }

  sink() #close log file
} #end function
