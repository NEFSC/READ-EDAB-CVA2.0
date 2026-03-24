#' @title Combine Multiple Rasters of Fisheries Data
#' @description A wrapper function for the \code{merge_rasts} function. This function adds log and skip functionality to help with batch runs and running \code{merge_rasts} in parallel.
#'
#' @param name Species name to add to log files and save data to correct directory (see vignette for recommended directory set up)
#' @param skip TRUE/FALSE indicating whether to skip creating the raster file if file already exists
#'
#' @return \code{combine_rasters_wrapper} returns the range of the rasterBrick returned by \code{merge_rasts}. This should be equal to 0 2, or else there are no presences in the dataset and the models will fail. This function will also save the resulting rasterBrick as a netcdf file in the species' input_rasters folder

combine_rasters_wrapper <- function(name, skip){
  sink(file.path(getwd(), 'logs', 'combineRasters.log'), append = T)
  # Ensure the sinks are closed when the function exits, regardless of how it exits.
  on.exit({
    #sink(type = "message")
    sink()
  })
  print(Sys.time())
  print(name)

  if(skip){
    if(file.exists(paste(file.path(getwd(),name, 'input_rasters'), 'combined_rasters.nc', sep = '/'))){
      print('file exists and skip == T, so skipping this file!')
      return(NA)
    } else {
      flist <-  dir(file.path(getwd(),name, 'input_rasters'), full.names = T)
      rasts <- vector(mode = 'list', length = length(flist))
      for(r in 1:length(flist)){
        rast <- raster::brick(flist[r]) #rast
        rasts[[r]] <- rast
        #print(r)
      }

      combinedRasts <- merge_rasts(rasts)
      print(range(combinedRasts[]))
      raster::writeRaster(combinedRasts, filename = paste(file.path(getwd(),name, 'input_rasters'), 'combined_rasters.nc', sep = '/'), bylayer = F, overwrite = T)
    } #end else
  } else { #if skip = F, do it anyway

    flist <-  dir(file.path(getwd(),name, 'input_rasters'), full.names = T)
    if(file.exists(paste(file.path(getwd(),name, 'input_rasters'), 'combined_rasters.nc', sep = '/'))){
      i <- which(flist == paste(file.path(getwd(),name, 'input_rasters'), 'combined_rasters.nc', sep = '/'))
      flist <- flist[-i]
    }
    rasts <- vector(mode = 'list', length = length(flist))
    for(r in 1:length(flist)){
      rast <- raster::brick(flist[r]) #rast
      rasts[[r]] <- rast
      #print(r)
    }

    combinedRasts <- merge_rasts(rasts)
    print(range(combinedRasts[]))
    raster::writeRaster(combinedRasts, filename = paste(file.path(getwd(),name, 'input_rasters'), 'combined_rasters.nc', sep = '/'), bylayer = F, overwrite = T)
    #sink()
  } #end else
}

