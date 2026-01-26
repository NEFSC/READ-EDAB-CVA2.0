#' @title Build and Predict Ensemble Model 
#' @description A wrapper function to build and predict the ensemble models using \code{make_sdm} and \code{make_predictions}. Progress is tracked via a log file. NOTE that this function does NOT include skip functionality. This is because it is quicker than the other predictions since it is performing a weighted average of the individual model predictions. 

#' @param spp Species name to add to log files and save data to correct directory (see vignette for recommended directory set up)
#' @param wMetric a character string containing one of two evaluation metrics that can be used to build the model weights: rmse or auc

#' @return returns the file path of the saved predictions. Predicted rasters, model weights, and preds from the ensemble are saved within specific directories. See the vignette for recommended directory set up.

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
