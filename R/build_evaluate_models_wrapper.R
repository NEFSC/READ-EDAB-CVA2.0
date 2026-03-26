#' @title Make Component Models
#' @description This is a wrapper function for all of the functions used to build and evaluate the component models: \code{make_sdm}, \code{sdm_cv}, \code{sdm_preds}, \code{sdm_eval}, and \code{sdm_importance}. This function produces log files and has skip functionality to assist in running multiple species in parallel.

#' @param spp Species name to add to log files and save data to correct directory (see vignette for recommended directory set up)
#' @param model Component model name. One of the following: gam, maxent, rf, brt, sdmtmb
#' @param skip TRUE/FALSE indicating whether to skip creating the raster file if file already exists

#' @return returns a summary object of the model built. Outputs from the subsequent functions called within are saved within specific directories. See the vignette for recommended directory set up.

build_evaluate_models_wrapper <- function(spp, model, skip){
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
      mod <- build_sdm(se = dfC, pa_col = 'value', xy_col = c('x', 'y'), month_col = 'month', year_col = 'year', model = model)
      save(mod, file = paste0(file.path(getwd(),spp, 'model_output', 'models'), '/', toupper(model), '.RData'))
    }
  } else {
    print(paste0(spp, '- making model - ', Sys.time()))
    mod <- build_sdm(se = dfC, pa_col = 'value', xy_col = c('x', 'y'), month_col = 'month', year_col = 'year', model = model)
    save(mod, file = paste0(file.path(getwd(),spp, 'model_output', 'models'), '/', toupper(model), '.RData'))
  }

  if(skip){
    if(file.exists(paste0(file.path(getwd(),spp, 'model_output', 'cvs'), '/', toupper(model), '.RData'))){
      print('cv exist and skip == T, so skipping this file!')
      #return(NA)
    } else {
      print(paste0(spp, '- performing CV - ', Sys.time()))
      load(paste0(file.path(getwd(),spp, 'model_output', 'models'), '/', toupper(model), '.RData'))
      cv <- cross_validate_sdm(mod = mod, se = dfC, pa_col = 'value', xy_col = c('x', 'y'), month_col = 'month', year_col = 'year', model = model)
      save(cv, file = paste0(file.path(getwd(),spp, 'model_output', 'cvs'), '/', toupper(model), '.RData'))
    }
  } else {
    print(paste0(spp, '- performing CV - ', Sys.time()))
    load(paste0(file.path(getwd(),spp, 'model_output', 'models'), '/', toupper(model), '.RData'))
    cv <- cross_validate_sdm(mod = mod, se = dfC, pa_col = 'value', xy_col = c('x', 'y'), month_col = 'month', year_col = 'year', model = model)
    save(cv, file = paste0(file.path(getwd(),spp, 'model_output', 'cvs'), '/', toupper(model), '.RData'))
  }

  if(skip){
    if(file.exists(paste0(file.path(getwd(),spp, 'model_output', 'preds'), '/',toupper(model), '.RData'))){
      print('preds exist and skip == T, so skipping this file!')
      #return(NA)
    } else {
      print(paste0(spp, '- Getting Preds - ', Sys.time()))
      load(paste0(file.path(getwd(),spp, 'model_output', 'cvs'), '/', toupper(model), '.RData'))
      preds <- pull_sdm_preds(cv = cv, model = model)
      save(preds, file = paste0(file.path(getwd(),spp, 'model_output', 'preds'), '/',toupper(model), '.RData'))
    }
  } else {
    print(paste0(spp, '- Getting Preds - ', Sys.time()))
    load(paste0(file.path(getwd(),spp, 'model_output', 'cvs'), '/', toupper(model), '.RData'))
    preds <- pull_sdm_preds(cv = cv, model = model)
    save(preds, file = paste0(file.path(getwd(),spp, 'model_output', 'preds'), '/',toupper(model), '.RData'))
  }

  if(skip){
    if(file.exists(paste0(file.path(getwd(),spp, 'model_output', 'eval_metrics'), '/',toupper(model), '.RData'))){
      print('eval metrics exist and skip == T, so skipping this file!')
      #return(NA)
    } else {
      print(paste0(spp, '- Evaluating Model - ', Sys.time()))
      load(paste0(file.path(getwd(),spp, 'model_output', 'preds'), '/',toupper(model), '.RData'))
      ev <- evaluate_sdm(preds = preds, metric = 'auc', model = model)
      save(ev, file = paste0(file.path(getwd(),spp, 'model_output', 'eval_metrics'), '/',toupper(model), '.RData'))
    }
  } else {
    print(paste0(spp, '- Evaluating Model - ', Sys.time()))
    load(paste0(file.path(getwd(),spp, 'model_output', 'preds'), '/',toupper(model), '.RData'))
    ev <- evaluate_sdm(preds = preds, metric = 'auc', model = model)
    save(ev, file = paste0(file.path(getwd(),spp, 'model_output', 'eval_metrics'), '/',toupper(model), '.RData'))
  }

  if(skip){
    if(file.exists(paste0(file.path(getwd(),spp, 'model_output', 'importance'), '/',toupper(model), '.RData'))){
      print('importance exists and skip == T, so skipping this file!')
      return(NA)
    } else {
      print(paste0(spp, '- Getting Variable Importance - ', Sys.time()))
      load(paste0(file.path(getwd(),spp, 'model_output', 'models'), '/', toupper(model), '.RData'))
      imp <- calculate_sdm_variable_importance(mod = mod, se = dfC, pa_col = 'value', xy_col = c('x', 'y'), month_col = 'month', year_col = 'year', model = model)
      save(imp, file = paste0(file.path(getwd(),spp, 'model_output', 'importance'), '/',toupper(model), '.RData'))
    }
  } else {
    print(paste0(spp, '- Getting Variable Importance - ', Sys.time()))
    load(paste0(file.path(getwd(),spp, 'model_output', 'models'), '/', toupper(model), '.RData'))
    imp <- calculate_sdm_variable_importance(mod = mod, se = dfC, pa_col = 'value', xy_col = c('x', 'y'), month_col = 'month', year_col = 'year', model = model)
    save(imp, file = paste0(file.path(getwd(),spp, 'model_output', 'importance'), '/',toupper(model), '.RData'))
  }


  #} #end y
  return(summary(mod))
}

