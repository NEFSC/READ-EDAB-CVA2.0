#' @title Wrapper Functions
#' @description Functions used to log and save outputs from processing steps. These functions rely on directory set up following the guidelines in the CVA2.0 Github Manual.

#' \itemize {
#' \item \code{combineRasters} is a wrapper function for \code{merge_rasts} that creates log files to track progress and with a skip functionality 
 #' \item \code{saveRast} is a wrapper function for \code{create_rast} that creates log files to track progress and with a skip functionality 
#' \item \code{getMOM6Hindcast} and \code{getMOM6Forecast} are wrapper functions that pull, average, and normalize the data using the MOM6 Data functions and save the output for each step.
#' \item \code{makeDF} is the wrapper function for \code{merge_spp_env}, \code{match_guilds}, \code{remove_corr}, and \code{clean_data} that produces logs and has skip functionality
#' \item \code {makeMods} is the wrapper function that saves the outputs from each of the model building functions (\code{make_sdm}, \code{sdm_cv}, \code{sdm_preds}, \code{sdm_eval}, and \code{sdm_importance}), logs progress, and has skip functionality
#' \item \code{predictMods} is the wrapper function for \code{make_predictions} that logs progress and has skip functionality
#' \item \code{makeEns} is a wrapper function to build and predict the ensemble models using \code{make_sdm} and \code{make_predictions}. Progress is tracked via a log file. NOTE that this function does NOT include skip functionality. This is because it is quicker than the other predictions since it is performing a weighted average of the individual model predictions. 
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
#' @param mMin,mMax,yMin,yMax (m)onth and (y)ear ranges to help create names for combined rasterBrick since they aren't always retained and these are necessary for matching fisheries and environmental data
#' @param wMetric a character string containing one of two evaluation metrics that can be used to build the model weights: rmse or auc

#' @return \code{saveRast} and \code{combineRasters} returns the range of the rasterBrick returned by \code{create_rast}. For \code{saveRast}, this should be equal to 0 2 if species is present in dataset, and 0 1 if species is not caught in dataset. For \code{combineRasters}, this should be equal to 0 2 or else there are no presences in the dataset and the models will fail. Both functions will also save the resulting rasterBrick as a netcdf file in the species' input_rasters folder
#' @return Both \code{getMOM6Hindcast} and \code{getMOM6Forecast} return a list whose length is equal to the number of variables supplied, where each item in the list is a rasterStack of data associated with that variable, with the number of layers equal to the number of time steps available 
#' @return \code{makeDF} does not return anything to the environment; it only saves the final data set (the output of \code{clean_data}) and writes to the log file
#' @return \code{make_predictions} returns a rasterStack with the same number of layers as in \code{rasts} with probability of occurance predicted for each day with available data. Values will range from 0 to 1. 
#' @return \code{makeMods} returns a summary object of the model built. Outputs from the subsequent functions called within are saved within specific directories. See the vignette for recommended directory set up.
#' @return \code{predictMods},\code{makeEns} each return the file path for the predicted model and ensemble predictions. Outputs from the functions called within are saved in specific directories. See the vignette for recommended directory set up.

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

makeDF <- function(name, skip, mMin, mMax, yMin, yMax){
  sink(file.path(getwd(), 'logs', 'build_data_frames.log'), append = T)
  # Ensure the sinks are closed when the function exits, regardless of how it exits.
  on.exit({
    #sink(type = "message")
    sink()
  })
  
  print(Sys.time())
  print(name)
  
  if(skip){
    if(file.exists(paste(file.path(getwd(),name), 'pa_clean.RData', sep = '/'))){
      print('file exists and skip == T, so skipping this species!')
      return(NA)
    } else {
      combinedRasts <- brick(paste(file.path(getwd(),name, 'input_rasters'), 'combined_rasters.nc', sep = '/')) #combinedRasts
      
      my <- expand.grid(mMin:mMax, yMin:yMax)
      names(combinedRasts) <- paste(sprintf("%02d", my$Var1), my$Var2, sep = '.')
      
      print(paste0(name, '- merging spp and env data -', Sys.time()))
      df <- merge_spp_env(rastStack = combinedRasts, envData = norm, addStatic = TRUE, staticData = './Data/staticVariables_cropped_normZ.RData')
      save(df, file = paste(file.path(getwd(),name), 'presence_absence_environment.RData', sep = '/'))
      
      print(paste0(name, '- matching guild -', Sys.time()))
      
      dfG <- match_guilds(spp_env = df, spp = c('Summer Flounder', 'Fluke', "PARALICHTHYS DENTATUS"), spp_col = 'SCI_NAME', spp_guild = 'spp_list.csv', feeding_key = 'feeding_guilds.csv', feeding_col = 'Feeding.Guild', habitat_key = 'habitat_guilds.csv',  habitat_col = 'Habitat.Guild', static_vars = c('x', 'y', 'month', 'year', 'bathy', 'rugosity', 'dist2coast'), pa_col = 'value')
      save(dfG, file = paste(file.path(getwd(),name), 'pa_guild.RData', sep = '/'))
      
      print(paste0(name, '- converting to binary -', Sys.time()))
      dfB <- clean_data(dfG, pa_col = 'value')
      save(dfB, file = paste(file.path(getwd(),name), 'pa_binary.RData', sep = '/'))
      
      print(paste0(name, ' - removing correlated variables -', Sys.time()))
      dfC <- remove_corr(dfB, pa_col = 'value', xy_col = c('x', 'y'), month_col = 'month', year_col = 'year')
      save(dfC, file = paste(file.path(getwd(),name), 'pa_clean.RData', sep = '/'))
    }
  } else {
    combinedRasts <- brick(paste(file.path(getwd(),name, 'input_rasters'), 'combined_rasters.nc', sep = '/')) #combinedRasts
    
    my <- expand.grid(mMin:mMax, yMin:yMax)
    names(combinedRasts) <- paste(sprintf("%02d", my$Var1), my$Var2, sep = '.')
    
    print(paste0(name, '- merging spp and env data -', Sys.time()))
    df <- merge_spp_env(rastStack = combinedRasts, envData = norm, addStatic = TRUE, staticData = './Data/staticVariables_cropped_normZ.RData')
    save(df, file = paste(file.path(getwd(),name), 'presence_absence_environment.RData', sep = '/'))
    
    print(paste0(name, '- matching guild -', Sys.time()))
    
    dfG <- match_guilds(spp_env = df, spp = c('Summer Flounder', 'Fluke', "PARALICHTHYS DENTATUS"), spp_col = 'SCI_NAME', spp_guild = 'spp_list.csv', feeding_key = 'feeding_guilds.csv', feeding_col = 'Feeding.Guild', habitat_key = 'habitat_guilds.csv',  habitat_col = 'Habitat.Guild', static_vars = c('x', 'y', 'month', 'year', 'bathy', 'rugosity', 'dist2coast'), pa_col = 'value')
    save(dfG, file = paste(file.path(getwd(),name), 'pa_guild.RData', sep = '/'))
    
    print(paste0(name, '- converting to binary -', Sys.time()))
    dfB <- clean_data(dfG, pa_col = 'value')
    save(dfB, file = paste(file.path(getwd(),name), 'pa_binary.RData', sep = '/'))
    
    print(paste0(name, ' - removing correlated variables -', Sys.time()))
    dfC <- remove_corr(dfB, pa_col = 'value', xy_col = c('x', 'y'), month_col = 'month', year_col = 'year')
    save(dfC, file = paste(file.path(getwd(),name), 'pa_clean.RData', sep = '/'))
  }
}

makeMods <- function(spp, model, skip){
  #open log file
  sink(file = file.path(getwd(), 'logs', paste0(model, '.log')), append = T)
  #sink(file = file.path(getwd(), 'logs', paste0(csvName, '.log')), append = T, type = 'message')
  
  # Ensure the sinks are closed when the function exits, regardless of how it exits.
  on.exit({
    #sink(type = "message")
    sink()
  })
  
  print(Sys.time())
  print(spp)
  
  load(paste(file.path(getwd(),spp), 'pa_clean.RData', sep = '/')) #load data - dfC
  
  #for(y in models){ 
  #print(y)
  
  if(skip){
    if(file.exists(paste0(file.path(getwd(),spp, 'model_output', 'models'), '/', toupper(model), '.RData'))){
      print('model exists and skip == T, so skipping this file!')
      #return(NA)
    } else {
      print(paste0(spp, '- making model - ', Sys.time()))
      mod <- make_sdm(se = dfC, pa_col = 'value', xy_col = c('x', 'y'), month_col = 'month', year_col = 'year', model = model)
      save(mod, file = paste0(file.path(getwd(),spp, 'model_output', 'models'), '/', toupper(model), '.RData'))
    }
  } else {
    print(paste0(spp, '- making model - ', Sys.time()))
    mod <- make_sdm(se = dfC, pa_col = 'value', xy_col = c('x', 'y'), month_col = 'month', year_col = 'year', model = model)
    save(mod, file = paste0(file.path(getwd(),spp, 'model_output', 'models'), '/', toupper(model), '.RData'))
  }
  
  if(skip){
    if(file.exists(paste0(file.path(getwd(),spp, 'model_output', 'cvs'), '/', toupper(model), '.RData'))){
      print('cv exist and skip == T, so skipping this file!')
      #return(NA)
    } else {
      print(paste0(spp, '- performing CV - ', Sys.time()))
      load(paste0(file.path(getwd(),spp, 'model_output', 'models'), '/', toupper(model), '.RData'))
      cv <- sdm_cv(mod = mod, se = dfC, pa_col = 'value', xy_col = c('x', 'y'), month_col = 'month', year_col = 'year', model = model)
      save(cv, file = paste0(file.path(getwd(),spp, 'model_output', 'cvs'), '/', toupper(model), '.RData'))
    }
  } else {
    print(paste0(spp, '- performing CV - ', Sys.time()))
    load(paste0(file.path(getwd(),spp, 'model_output', 'models'), '/', toupper(model), '.RData'))
    cv <- sdm_cv(mod = mod, se = dfC, pa_col = 'value', xy_col = c('x', 'y'), month_col = 'month', year_col = 'year', model = model)
    save(cv, file = paste0(file.path(getwd(),spp, 'model_output', 'cvs'), '/', toupper(model), '.RData'))
  }
  
  if(skip){
    if(file.exists(paste0(file.path(getwd(),spp, 'model_output', 'preds'), '/',toupper(model), '.RData'))){
      print('preds exist and skip == T, so skipping this file!')
      #return(NA)
    } else {
      print(paste0(spp, '- Getting Preds - ', Sys.time()))
      load(paste0(file.path(getwd(),spp, 'model_output', 'cvs'), '/', toupper(model), '.RData'))
      preds <- sdm_preds(cv = cv, model = model)
      save(preds, file = paste0(file.path(getwd(),spp, 'model_output', 'preds'), '/',toupper(model), '.RData'))
    }
  } else {
    print(paste0(spp, '- Getting Preds - ', Sys.time()))
    load(paste0(file.path(getwd(),spp, 'model_output', 'cvs'), '/', toupper(model), '.RData'))
    preds <- sdm_preds(cv = cv, model = model)
    save(preds, file = paste0(file.path(getwd(),spp, 'model_output', 'preds'), '/',toupper(model), '.RData'))
  }
  
  if(skip){
    if(file.exists(paste0(file.path(getwd(),spp, 'model_output', 'eval_metrics'), '/',toupper(model), '.RData'))){
      print('eval metrics exist and skip == T, so skipping this file!')
      #return(NA)
    } else {
      print(paste0(spp, '- Evaluating Model - ', Sys.time()))
      load(paste0(file.path(getwd(),spp, 'model_output', 'preds'), '/',toupper(model), '.RData'))
      ev <- sdm_eval(preds = preds, metric = 'auc', model = model)
      save(ev, file = paste0(file.path(getwd(),spp, 'model_output', 'eval_metrics'), '/',toupper(model), '.RData'))
    }
  } else {
    print(paste0(spp, '- Evaluating Model - ', Sys.time()))
    load(paste0(file.path(getwd(),spp, 'model_output', 'preds'), '/',toupper(model), '.RData'))
    ev <- sdm_eval(preds = preds, metric = 'auc', model = model)
    save(ev, file = paste0(file.path(getwd(),spp, 'model_output', 'eval_metrics'), '/',toupper(model), '.RData'))
  }
  
  if(skip){
    if(file.exists(paste0(file.path(getwd(),spp, 'model_output', 'importance'), '/',toupper(model), '.RData'))){
      print('importance exists and skip == T, so skipping this file!')
      return(NA)
    } else {
      print(paste0(spp, '- Getting Variable Importance - ', Sys.time()))
      load(paste0(file.path(getwd(),spp, 'model_output', 'models'), '/', toupper(model), '.RData'))
      imp <- sdm_importance(mod = mod, se = dfC, pa_col = 'value', xy_col = c('x', 'y'), month_col = 'month', year_col = 'year', model = model)
      save(imp, file = paste0(file.path(getwd(),spp, 'model_output', 'importance'), '/',toupper(model), '.RData'))
    }
  } else {
    print(paste0(spp, '- Getting Variable Importance - ', Sys.time()))
    load(paste0(file.path(getwd(),spp, 'model_output', 'models'), '/', toupper(model), '.RData'))
    imp <- sdm_importance(mod = mod, se = dfC, pa_col = 'value', xy_col = c('x', 'y'), month_col = 'month', year_col = 'year', model = model)
    save(imp, file = paste0(file.path(getwd(),spp, 'model_output', 'importance'), '/',toupper(model), '.RData'))
  }
  
  
  #} #end y 
  return(summary(mod))
}

predictMods <- function(spp, model, skip){
  #open log file
  sink(file = file.path(getwd(), 'logs', paste0(model, '_prediction.log')), append = T)
  #sink(file = file.path(getwd(), 'logs', paste0(csvName, '.log')), append = T, type = 'message')
  
  # Ensure the sinks are closed when the function exits, regardless of how it exits.
  on.exit({
    #sink(type = "message")
    sink()
  })
  
  print(Sys.time())
  print(spp)
  
  load(paste0(file.path(getwd(),spp, 'model_output', 'models'), '/', toupper(model), '.RData')) #load model - mod
  load(paste(file.path(getwd(),spp), 'pa_clean.RData', sep = '/')) #load data - dfC
  
  if(skip){
    if(file.exists(paste(file.path(getwd(),spp, 'output_rasters'), '/', toupper(model), '.RData', sep = ''))){
      print('prediction already exists & skip == T, so skipping')
    } else {
      #predict model
      print(paste0(spp, '- predicting ', model, ' - ', Sys.time()))
      abund <- make_predictions(mod = mod, model = model, rasts = norm, mask = T, bathy_nm = 'bathy', bathy_max = 1000, se = dfC, staticData = './Data/staticVariables_cropped_normZ.RData', xy_col = c('x', 'y'), month_col = 'month', year_col = 'year')
      save(abund, file = paste(file.path(getwd(),spp, 'output_rasters'), '/', toupper(model), '.RData', sep = ''))
    } 
  } else {
    #predict model
    print(paste0(spp, '- predicting ', model, ' - ', Sys.time()))
    abund <- make_predictions(mod = mod, model = model, rasts = norm, mask = T, bathy_nm = 'bathy', bathy_max = 1000, se = dfC, staticData = './Data/staticVariables_cropped_normZ.RData', xy_col = c('x', 'y'), month_col = 'month', year_col = 'year')
    save(abund, file = paste(file.path(getwd(),spp, 'output_rasters'), '/', toupper(model), '.RData', sep = ''))
  }
  
  return(paste(file.path(getwd(),spp, 'output_rasters'), '/', toupper(model), '.RData', sep = ''))
}

makeEns <- function(spp, wMetric){
  #open log file
  sink(file = file.path(getwd(), 'logs', 'ensemble.log'), append = T)
  #sink(file = file.path(getwd(), 'logs', paste0(csvName, '.log')), append = T, type = 'message')
  
  # Ensure the sinks are closed when the function exits, regardless of how it exits.
  on.exit({
    #sink(type = "message")
    sink()
  })
  
  print(Sys.time())
  print(spp)
  
  print('generating weights...')
  #load in evaluation metrics
  evalFlist <- dir(file.path(getwd(),spp, 'model_output', 'eval_metrics'), pattern = '.RData', full.names = T)
  eval <- vector(length = length(evalFlist))
  for(y in 1:length(evalFlist)){
    load(evalFlist[y])
    eval[y] <- ev
  }
  
  #generate weights
  if(wMetric == 'rmse'){
    weights <- EFHSDM::MakeEnsemble(rmse = eval)
  }
  if(wMetric == 'auc'){
    weights <- eval/sum(eval) #we need to make weights like this since AUC bigger = better; whereas RMSE smaller = better
  }
  save(weights, file = paste(file.path(getwd(),spp, 'model_output'), 'ensemble_weights.RData', sep = '/'))
  
  print('making ensemble...')
  #pull preds 
  predFlist <- dir(file.path(getwd(),spp, 'model_output', 'preds'), pattern = '.RData', full.names = T)
  pds <- vector(mode = 'list', length = length(predFlist))
  for(y in 1:length(predFlist)){
    load(predFlist[y])
    pds[[y]] <- preds[complete.cases(preds),]
  }
  
  #make ensemble 
  ens <- make_sdm(model = 'ens', ensembleWeights = weights, ensemblePreds = pds)
  save(ens, file = paste(file.path(getwd(),spp, 'model_output', 'models'), 'ENSEMBLE.RData', sep = '/'))
  
  print('predicting ensemble...')
  
  abundFlist <- dir(file.path(getwd(),spp, 'output_rasters'), pattern = '.RData', full.names = T)
  #test if an ensemble object exists because we don't want to include that
  i <- grep('ENSEMBLE', abundFlist)
  if(length(i) != 0){
    abundFlist <- abundFlist[-i]
  }
  abds <- vector(mode = 'list', length = length(abundFlist))
  for(y in 1:length(abundFlist)){
    load(abundFlist[y])
    abds[[y]] <- abund
  }
  
  
  abund <- make_predictions(model = 'ens', rasts = abds, weights = weights, staticData = NULL, mask = F, bathy_nm = NULL, bathy_max = NULL, se = NULL, month_col = NULL, year_col = NULL, xy_col = NULL)
  nms <- expand.grid(month.abb, 1993:2019)
  names(abund) <- paste(nms$Var1, nms$Var2, sep = '.')
  save(abund, file = paste(file.path(getwd(),spp, 'output_rasters'), 'ENSEMBLE.RData', sep = '/'))
  print(Sys.time())
  
  return(paste(file.path(getwd(),spp, 'output_rasters'), 'ENSEMBLE.RData', sep = '/'))
}
