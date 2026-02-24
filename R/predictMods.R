#' @title Predict Models
#' @description The wrapper function for \code{make_predictions}. This function produces log files and has skip functionality to assist in running multiple species in parallel.

#' @param spp Species name to add to log files and save data to correct directory (see vignette for recommended directory set up)
#' @param model Component model name. One of the following: gam, maxent, rf, brt, sdmtmb
#' @param skip TRUE/FALSE indicating whether to skip creating the raster file if file already exists

#' @return returns the file path of the saved predictions. Outputs are saved within specific directories. See the vignette for recommended directory set up.

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
      abund <- make_predictions(mod = mod, model = model, rasts = norm, mask = T, bathyR = bathyR, bathy_max = 1000, se = dfC,  staticVars = staticVars, xy_col = c('x', 'y'), month_col = 'month', year_col = 'year')
      save(abund, file = paste(file.path(getwd(),spp, 'output_rasters'), '/', toupper(model), '.RData', sep = ''))
    }
  } else {
    #predict model
    print(paste0(spp, '- predicting ', model, ' - ', Sys.time()))
    abund <- make_predictions(mod = mod, model = model, rasts = norm, mask = T, bathyR = bathyR, bathy_max = 1000, se = dfC, staticVars = staticVars, xy_col = c('x', 'y'), month_col = 'month', year_col = 'year')
    save(abund, file = paste(file.path(getwd(),spp, 'output_rasters'), '/', toupper(model), '.RData', sep = ''))
  }

  return(paste(file.path(getwd(),spp, 'output_rasters'), '/', toupper(model), '.RData', sep = ''))
}
