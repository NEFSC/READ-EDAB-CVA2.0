#' @title Pull MOM6 Forecast Data
#' @description \code{getMOM6Forecast} is a wrapper function that pulls, averages, and normalizes the MOM6 forecast data using the \code{pull_forecast}, \code{avg_env}, \code{sd_env}, and \code{norm_env} functions, and saves the output from each step. This function adds log functionality to help with batch runs and running the encompassing functions in parallel. Note that unlike similar functions, there is no skip functionality since this function is meant to only be run once.

#' @param varDF a data.frame with the columns Long.Name and Short.Name to be used as reqVars and shortNames, respectively
#' @param inPar TRUE/FALSE to determine if data pulls and averaging should be conducted in parallel with dopar
#' @param jsonURL URL pointing to JSON table variable lists for desired MOM6 run type and domain
#' @param release release code. Must match one of the options in the 'cefi_release' column in provided JSON table
#' @param init initialization code. Must match one of the options in the 'cefi_init_date' column in provided JSON table. For forecast only
#' @param ens ensemble member. Must be equal to 1-10. For decadal forecasts, different ensemble members represent slightly different forcing scenarios. For forecast only.

#' @return the output from \code{norm_env} - a list whose length is equal to the number of variables supplied, where each item in the list is a rasterStack of data associated with that variable, with the number of layers equal to the number of time steps available. The function also saves the output from each step - see the package website for necessary directory set up.

getMOM6Forecast <- function(varDF, inPar = TRUE, jsonURL, release, init, ens){
  raw <- avg <- sds <- norm <- vector(mode = 'list', length = nrow(varDF))

  if(inPar == TRUE){
    cluster <- parallel::makeCluster(10, type='PSOCK')
    doParallel::registerDoParallel(cluster)
    raw <- foreach::foreach(x = 1:nrow(varDF), .packages = c("ncdf4", 'raster', 'jsonlite')) %dopar% {
      #for(x in 1:nrow(var.list)){
      r <- pull_forecast(varURL = jsonURL, reqVars = varDF$Long.Name[x], shortNames = varDF$Short.Name[x], release = release, init = init, ens = ens)
      raw[[x]] <- r
    }
    parallel::stopCluster(cluster)
  } else {
    for(x in 1:nrow(varDF)){
      r <- pull_forecast(varURL = jsonURL, reqVars = varDF$Long.Name[x], shortNames = varDF$Short.Name[x], release = release, init = init, ens = ens)
      raw[[x]] <- r
    }
  }
  names(raw) <- varDF$Short.Name
  save(raw, file = './Data/MOM6/raw_MOM6_forecast.RData')

  if(inPar == TRUE){
    cluster <- parallel::makeCluster(10, type='PSOCK')
    doParallel::registerDoParallel(cluster)
    avg <- foreach::foreach(x = 1:nrow(varDF), .packages = c("ncdf4", 'raster', 'jsonlite')) %do% {
      #for(x in 1:nrow(varDF)){
      a <- avg_env(raw[[x]])
      avg[[x]] <- a
    }
    parallel::stopCluster(cluster)
  } else {
    for(x in 1:nrow(varDF)){
      a <- avg_env(raw[[x]])
      avg[[x]] <- a
    }
  }
  names(avg) <- varDF$Short.Name
  save(avg, file = './Data/MOM6/avg_MOM6_forecast.RData')


  if(inPar == TRUE){
    cluster <- parallel::makeCluster(10, type='PSOCK')
    doParallel::registerDoParallel(cluster)
    sds <- foreach::foreach(y = 1:nrow(varDF), .packages = c("ncdf4", 'raster', 'jsonlite')) %do% {
      #for(y in 1:nrow(varDF)){
      s <- sd_env(raw[[y]])
      sds[[y]] <- s
    }
    parallel::stopCluster(cluster)
  } else {
    for(y in 1:nrow(varDF)){
      s <- sd_env(raw[[y]])
      sds[[y]] <- s
    }
  }
  names(sds) <- varDF$Short.Name
  save(sds, file = './Data/MOM6/sd_MOM6_forecast.RData')

  if(inPar == TRUE){
    cluster <- parallel::makeCluster(10, type='PSOCK')
    doParallel::registerDoParallel(cluster)
    norm <- foreach::foreach(x = 1:nrow(varDF), .packages = c("ncdf4", 'raster', 'jsonlite', 'abind')) %dopar% {
      #for(y in 1:nrow(varDF)){
      n <- norm_env(rawList = raw[[x]], avgList = avg[[x]], sdList = sds[[x]], shortNames = varDF$Short.Name[x])
      norm[[x]] <- n
    }
    parallel::stopCluster(cluster)
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
