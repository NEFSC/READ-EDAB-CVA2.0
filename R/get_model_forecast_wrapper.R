#' @title Pull, Average, and Normalize MOM6 Forecast Data
#' @description A wrapper function that pulls, averages, and normalizes the MOM6 forecast data using the \code{pull_forecast}, \code{avg_env}, \code{sd_env}, and \code{norm_env} functions, and saves the output from each step. This function adds log functionality to help with batch runs and running the encompassing functions in parallel. Note that unlike similar functions, there is no skip functionality since this function is meant to only be run once.

#' @param var_df a data.frame with the columns Long.Name and Short.Name to be used as reqVars and shortNames, respectively
#' @param in_par TRUE/FALSE to determine if data pulls and averaging should be conducted in parallel with dopar
#' @param json_url URL pointing to JSON table variable lists for desired MOM6 run type and domain
#' @param release release code. Must match one of the options in the 'cefi_release' column in provided JSON table
#' @param init initialization code. Must match one of the options in the 'cefi_init_date' column in provided JSON table. For forecast only
#' @param ens ensemble member. Must be equal to 1-10. For decadal forecasts, different ensemble members represent slightly different forcing scenarios. For forecast only.

#' @return the output from \code{norm_env} - a list whose length is equal to the number of variables supplied, where each item in the list is a rasterStack of data associated with that variable, with the number of layers equal to the number of time steps available. The function also saves the output from each step - see the package website for necessary directory set up.

get_model_forecast_wrapper <- function(var_df, in_par = TRUE, json_url, release, init, ens){
  raw <- avg <- sds <- norm <- vector(mode = 'list', length = nrow(var_df))

  if(in_par == TRUE){
    cluster <- parallel::makeCluster(10, type='PSOCK')
    doParallel::registerDoParallel(cluster)
    raw <- foreach::foreach(x = 1:nrow(var_df), .packages = c("ncdf4", 'raster', 'jsonlite')) %dopar% {
      #for(x in 1:nrow(var.list)){
      r <- pull_forecast(varURL = json_url, reqVars = var_df$Long.Name[x], shortNames = var_df$Short.Name[x], release = release, init = init, ens = ens)
      raw[[x]] <- r
    }
    parallel::stopCluster(cluster)
  } else {
    for(x in 1:nrow(var_df)){
      r <- pull_forecast(varURL = json_url, reqVars = var_df$Long.Name[x], shortNames = var_df$Short.Name[x], release = release, init = init, ens = ens)
      raw[[x]] <- r
    }
  }
  names(raw) <- var_df$Short.Name
  save(raw, file = './Data/MOM6/raw_MOM6_forecast.RData')

  if(in_par == TRUE){
    cluster <- parallel::makeCluster(10, type='PSOCK')
    doParallel::registerDoParallel(cluster)
    avg <- foreach::foreach(x = 1:nrow(var_df), .packages = c("ncdf4", 'raster', 'jsonlite')) %do% {
      #for(x in 1:nrow(var_df)){
      a <- avg_model_data(raw[[x]])
      avg[[x]] <- a
    }
    parallel::stopCluster(cluster)
  } else {
    for(x in 1:nrow(var_df)){
      a <- avg_model_data(raw[[x]])
      avg[[x]] <- a
    }
  }
  names(avg) <- var_df$Short.Name
  save(avg, file = './Data/MOM6/avg_MOM6_forecast.RData')


  if(in_par == TRUE){
    cluster <- parallel::makeCluster(10, type='PSOCK')
    doParallel::registerDoParallel(cluster)
    sds <- foreach::foreach(y = 1:nrow(var_df), .packages = c("ncdf4", 'raster', 'jsonlite')) %do% {
      #for(y in 1:nrow(var_df)){
      s <- sd_env(raw[[y]])
      sds[[y]] <- s
    }
    parallel::stopCluster(cluster)
  } else {
    for(y in 1:nrow(var_df)){
      s <- sd_env(raw[[y]])
      sds[[y]] <- s
    }
  }
  names(sds) <- var_df$Short.Name
  save(sds, file = './Data/MOM6/sd_MOM6_forecast.RData')

  if(in_par == TRUE){
    cluster <- parallel::makeCluster(10, type='PSOCK')
    doParallel::registerDoParallel(cluster)
    norm <- foreach::foreach(x = 1:nrow(var_df), .packages = c("ncdf4", 'raster', 'jsonlite', 'abind')) %dopar% {
      #for(y in 1:nrow(var_df)){
      n <- normalize_model_data(raw_list = raw[[x]], avg_list = avg[[x]], sd_list = sds[[x]], short_names = var_df$Short.Name[x])
      norm[[x]] <- n
    }
    parallel::stopCluster(cluster)
  } else {
    for(y in 1:nrow(var_df)){
      n <- normalize_model_data(raw_list = raw[[x]], avg_list = avg[[x]], sd_list = sds[[x]], short_names = var_df$Short.Name[x])
      norm[[y]] <- n
    }
  }
  names(norm) <- var_df$Short.Name
  save(norm, file = './Data/MOM6/norm_MOM6_forecast.RData')

  return(norm)
} #end function
