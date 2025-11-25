#' @title Wrapper Functions
#' @description Functions used to log and save outputs from processing steps. These functions rely on directory set up following the guidelines in the CVA2.0 Github Manual.

#' \itemize {
#' \item \code{combineRasters} is a wrapper function for \code{merge_rasts} that creates log files to track progress and with a skip functionality 
 #' \item \code{saveRast} is a wrapper function for \code{create_rast} that creates log files to track progress and with a skip functionality 
#' \item \code{getMOM6Hindcast} and \code{getMOM6Forecast} are wrapper functions that pull, average, and normalize the data using the MOM6 Data functions and save the output for each step.
}

#' @param csvName character string indicating which csv to use to create raster without extension (i.e. 'example' NOT 'example.csv')
#' @param isObs TRUE/FALSE indicating whether or not the data is observer or similar fisheries-dependent data. Will force check to determine if species is observed at least 30 times throughout timeseries before creating raster
#' @param grid static link to a ncdcf object with the variables lon, lat, time - can be link to remote data - must be able to be read with nc_open
#' @param spp,name Species name to add to log files and save data to correct directory (see vignette for recommended directory set up)
#' @param sppNames a vector containing all possible names for the target species. Must have a length >= 1
#' @param skip TRUE/FALSE indicating whether to skip creating the raster file if file already exists

#' @param varDF a data.frame with the columns Long.Name and Short.Name to be used as reqVars and shortNames, respectively  
#' @param inPar TRUE/FALSE to determine if data pulls and averaging such be conducted in parallel with dopar 
#' @param jsonURL URL pointing to JSON table variable lists for desired MOM6 run type (hindcast, decadal forecast) and domain 
#' @param release release code. Must match one of the options in the 'cefi_release' column in provided JSON table 
#' @param init initialization code. Must match one of the options in the 'cefi_init_date' column in provided JSON table. For forecast only
#' @param ens ensemble member. Must be equal to 1-10. For decadal forecasts, different ensemble members represent slightly different forcing scenarios. For forecast only.

#' @return \code{saveRast} and \code{combineRasters} returns the range of the rasterBrick returned by \code{create_rast}. For \code{saveRast}, this should be equal to 0 2 if species is present in dataset, and 0 1 if species is not caught in dataset. For \code{combineRasters}, this should be equal to 0 2 or else there are no presences in the dataset and the models will fail. Both functions will also save the resulting rasterBrick as a netcdf file in the species' input_rasters folder
#' @return Both \code{getMOM6Hindcast} and \code{getMOM6Forecast} return a list whose length is equal to the number of variables supplied, where each item in the list is a rasterStack of data associated with that variable, with the number of layers equal to the number of time steps available 

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

combineRasters <- function(name, skip){
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
        rast <- brick(flist[r]) #rast
        rasts[[r]] <- rast
        #print(r)
      }
      
      combinedRasts <- merge_rasts(rasts)
      print(range(combinedRasts[]))
      writeRaster(combinedRasts, filename = paste(file.path(getwd(),name, 'input_rasters'), 'combined_rasters.nc', sep = '/'), bylayer = F, overwrite = T)
    } #end else
  } else { #if skip = F, do it anyway
    
    flist <-  dir(file.path(getwd(),name, 'input_rasters'), full.names = T)
    if(file.exists(paste(file.path(getwd(),name, 'input_rasters'), 'combined_rasters.nc', sep = '/'))){
      i <- which(flist == paste(file.path(getwd(),name, 'input_rasters'), 'combined_rasters.nc', sep = '/'))
      flist <- flist[-i]
    }
    rasts <- vector(mode = 'list', length = length(flist))
    for(r in 1:length(flist)){
      rast <- brick(flist[r]) #rast
      rasts[[r]] <- rast
      #print(r)
    }
    
    combinedRasts <- merge_rasts(rasts)
    print(range(combinedRasts[]))
    writeRaster(combinedRasts, filename = paste(file.path(getwd(),name, 'input_rasters'), 'combined_rasters.nc', sep = '/'), bylayer = F, overwrite = T)
    #sink()
  } #end else 
}


getMOM6Hindcast <- function(varDF, inPar = TRUE, jsonURL, release){
  raw <- avg <- sds <- norm <- vector(mode = 'list', length = nrow(varDF))
  
  if(inPar == TRUE){
    cluster <- makeCluster(10, type='PSOCK')
    registerDoParallel(cluster)
    raw <- foreach(x = 1:nrow(var.list), .packages = c("ncdf4", 'raster', 'jsonlite')) %dopar% {
      #for(x in 1:nrow(var.list)){
      r <- pull_hind(jsonURL = jsonURL, reqVars = varDF$Long.Name[x], shortNames = varDF$Short.Name[x], release = release)
      raw[[x]] <- r
    }
    stopCluster(cluster)
  } else {
    for(x in 1:nrow(var.list)){
      r <- pull_hind(varURL = jsonURL, reqVars = varDF$Long.Name[x], shortNames = varDF$Short.Name[x], release = release)
      raw[[x]] <- r
    }
  }
  names(raw) <- varDF$Short.Name
  save(raw, file = './Data/MOM6/raw_MOM6_hindcast.RData')
  
  if(inPar == TRUE){
    cluster <- makeCluster(10, type='PSOCK')
    registerDoParallel(cluster)
    avg <- foreach(x = 1:nrow(var.list), .packages = c("ncdf4", 'raster', 'jsonlite')) %do% {
      #for(x in 1:nrow(var.list)){
      a <- avg_env(raw[[x]])
      avg[[x]] <- a
    }
    stopCluster(cluster)
  } else {
    for(x in 1:nrow(var.list)){
      a <- avg_env(raw[[x]])
      avg[[x]] <- a
    }
  }
  names(avg) <- varDF$Short.Name
  save(avg, file = './Data/MOM6/avg_MOM6_hindcast.RData')
  

  if(inPar == TRUE){
    cluster <- makeCluster(10, type='PSOCK')
    registerDoParallel(cluster)
    sds <- foreach(y = 1:nrow(var.list), .packages = c("ncdf4", 'raster', 'jsonlite')) %do% {
    #for(y in 1:nrow(var.list)){
      s <- sd_env(raw[[y]])
      sds[[y]] <- s
    }
    stopCluster(cluster)
  } else {
    for(y in 1:nrow(var.list)){
      s <- sd_env(raw[[y]])
      sds[[y]] <- s
    }
  }
  names(sds) <- varDF$Short.Name
  save(sds, file = './Data/MOM6/sd_MOM6_hindcast.RData')
  
  if(inPar == TRUE){
    cluster <- makeCluster(10, type='PSOCK')
    registerDoParallel(cluster)
    norm <- foreach(x = 1:nrow(var.list), .packages = c("ncdf4", 'raster', 'jsonlite', 'abind')) %dopar% {
    #for(y in 1:nrow(var.list)){
      n <- norm_env(rawList = raw[[x]], avgList = avg[[x]], sdList = sds[[x]], shortNames = varDF$Short.Name[x])
      norm[[x]] <- n
    }
    stopCluster(cluster)
  } else {
    for(y in 1:nrow(var.list)){
      n <- norm_env(rawList = raw[[y]], avgList = avg[[y]], sdList = sds[[y]], shortNames = varDF$Short.Name[y])
      norm[[y]] <- n
    }
  }
  names(norm) <- varDF$Short.Name
  save(norm, file = './Data/MOM6/norm_MOM6_hindcast.RData')
  
  return(norm)
} #end function

getMOM6Forecast <- function(varDF, inPar = TRUE, jsonURL, release, init, ens){
  raw <- avg <- sds <- norm <- vector(mode = 'list', length = nrow(varDF))
  
  if(inPar == TRUE){
    cluster <- makeCluster(10, type='PSOCK')
    registerDoParallel(cluster)
    raw <- foreach(x = 1:nrow(varDF), .packages = c("ncdf4", 'raster', 'jsonlite')) %dopar% {
      #for(x in 1:nrow(var.list)){
      r <- pull_forecast(jsonURL = jsonURL, reqVars = varDF$Long.Name[x], shortNames = varDF$Short.Name[x], release = release, init = init, ens = ens)
      raw[[x]] <- r
    }
    stopCluster(cluster)
  } else {
    for(x in 1:nrow(varDF)){
      r <- pull_forecast(jsonURL = jsonURL, reqVars = varDF$Long.Name[x], shortNames = varDF$Short.Name[x], release = release, init = init, ens = ens)
      raw[[x]] <- r
    }
  }
  names(raw) <- varDF$Short.Name
  save(raw, file = './Data/MOM6/raw_MOM6_forecast.RData')
  
  if(inPar == TRUE){
    cluster <- makeCluster(10, type='PSOCK')
    registerDoParallel(cluster)
    avg <- foreach(x = 1:nrow(varDF), .packages = c("ncdf4", 'raster', 'jsonlite')) %do% {
      #for(x in 1:nrow(varDF)){
      a <- avg_env(raw[[x]])
      avg[[x]] <- a
    }
    stopCluster(cluster)
  } else {
    for(x in 1:nrow(varDF)){
      a <- avg_env(raw[[x]])
      avg[[x]] <- a
    }
  }
  names(avg) <- varDF$Short.Name
  save(avg, file = './Data/MOM6/avg_MOM6_forecast.RData')
  
  
  if(inPar == TRUE){
    cluster <- makeCluster(10, type='PSOCK')
    registerDoParallel(cluster)
    sds <- foreach(y = 1:nrow(varDF), .packages = c("ncdf4", 'raster', 'jsonlite')) %do% {
      #for(y in 1:nrow(varDF)){
      s <- sd_env(raw[[y]])
      sds[[y]] <- s
    }
    stopCluster(cluster)
  } else {
    for(y in 1:nrow(varDF)){
      s <- sd_env(raw[[y]])
      sds[[y]] <- s
    }
  }
  names(sds) <- varDF$Short.Name
  save(sds, file = './Data/MOM6/sd_MOM6_forecast.RData')
  
  if(inPar == TRUE){
    cluster <- makeCluster(10, type='PSOCK')
    registerDoParallel(cluster)
    norm <- foreach(x = 1:nrow(varDF), .packages = c("ncdf4", 'raster', 'jsonlite', 'abind')) %dopar% {
      #for(y in 1:nrow(varDF)){
      n <- norm_env(rawList = raw[[x]], avgList = avg[[x]], sdList = sds[[x]], shortNames = varDF$Short.Name[x])
      norm[[x]] <- n
    }
    stopCluster(cluster)
  } else {
    for(y in 1:nrow(varDF)){
      n <- norm_env(rawList = raw[[y]], avgList = avg[[y]], sdList = sds[[y]], shortNames = varDF$Short.Name[y])
      norm[[y]] <- n
    }
  }
  names(norm) <- varDF$Short.Name
  save(norm, file = './Data/MOM6/norm_MOM6_forecast.RData')
  
  return(norm)
} #end function