#' @title Build and Predict Ensemble Model
#'
#' @description A wrapper function to build and predict the ensemble models using \code{build_sdm} and \code{make_sdm_predictions}. Progress is tracked via a log file. NOTE that this function does NOT include skip functionality. This is because it is quicker than the other predictions since it is performing a weighted average of the individual model predictions.

#' @param spp Species name to add to log files and save data to correct directory (see vignette for recommended directory set up)
#' @param yr_min start of year range to specify which predictions to use
#' @param yr_max end of year range to specify which predictions to use
#' @param build_ens TRUE/FALSE to switch on/off ensemble model building, set to FALSE to just predict ensemble

#' @return returns the file path of the saved predictions. Predicted rasters, model weights, and preds from the ensemble are saved within specific directories. See the manual for recommended directory set up.

build_ensemble_wrapper <- function(spp, yr_min, yr_max, build_ens = T){
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

  if(build_ens){
    print('generating weights...')
    #load in evaluation metrics
    evalFlist <- dir(file.path(getwd(),spp, 'model_output', 'eval_metrics'), pattern = '.RData', full.names = T)
    eval <- vector(length = length(evalFlist))
    for(y in 1:length(evalFlist)){
      load(evalFlist[y])
      eval[y] <- ev
    }

    #generate weights
    weights <- eval/sum(eval) #we need to make weights like this since AUC bigger = better; whereas RMSE smaller = better
    save(weights, file = paste(file.path(getwd(),spp, 'model_output'), 'ensemble_weights.RData', sep = '/'))

    print('making ensemble...')
    load(paste(file.path(getwd(),spp), 'pa_clean.RData', sep = '/')) #load data - dfC
    #names(dfC)[which(names(dfC) == 'value')] <- 'abund'
    #pull preds
    predFlist <- dir(file.path(getwd(),spp, 'model_output', 'preds'), pattern = '.RData', full.names = T)
    pds <- vector(mode = 'list', length = length(predFlist))
    ind <- NULL
    for(y in 1:length(predFlist)){
      load(predFlist[y])
      if(grepl('BRT', predFlist[y])){
        preds <- merge(preds[!duplicated(preds[c('x', 'y', 'month', 'year')]),c('x', 'y', 'abund', 'month', 'year', 'pred')], dfC, by = c('x', 'y', 'month', 'year'), all.y = T)
      }
      if(grepl('SDMTMB', predFlist[y])){
        names(preds)[which(names(preds) == 'value')] <- 'abund'
        preds <- merge(preds[!duplicated(preds[c('x', 'y', 'month', 'year')]),c('abund', 'month', 'year', 'x', 'y', 'pred')], dfC,  by = c('x', 'y', 'month', 'year'), all.y = T)
      }
      if(grepl('RF', predFlist[y])){
        names(preds)[which(names(preds) == 'X')] <- 'x'
        names(preds)[which(names(preds) == 'Y')] <- 'y'
        preds <- merge(preds[!duplicated(preds[c('x', 'y', 'month', 'year')]),c('abund', 'month', 'year', 'x', 'y', 'pred')], dfC, by = c('x', 'y', 'month', 'year'), all.y = T)
      }
      pds[[y]] <- preds

      ##make index to subset preds to just locations where all models have presence/absence and predicted values
      i <- !is.na(preds$abund) & !is.na(preds$pred)
      i <- replace(i, i == TRUE, 1)
      i <- replace(i, i == FALSE, NA)
      ind <- cbind(ind, i)
    }

    ind2 <- apply(ind, 1, sum, na.rm = T)
    pds <- lapply(pds, FUN = function(x){x[which(ind2 == length(weights)),]})

    #make ensemble
    ens <- build_sdm(model = 'ens', ensemble_weights = weights, ensemble_preds = pds)
    save(ens, file = paste(file.path(getwd(),spp, 'model_output', 'models'), 'ENSEMBLE.RData', sep = '/'))
  } #end buildEns

  load(paste(file.path(getwd(),spp, 'model_output'), 'ensemble_weights.RData', sep = '/')) #weights

  print('predicting ensemble...')

  abundFlist <- dir(file.path(getwd(),spp, 'output_rasters'), pattern = paste0(yr_min, '_', yr_max, '.RData'), full.names = T)
  #test if an ensemble object exists because we don't want to include that
  i <- grep('ENSEMBLE', abundFlist)
  if(length(i) != 0){
    abundFlist <- abundFlist[-i]
  }
  abds <- vector(mode = 'list', length = length(abundFlist))
  for(y in 1:length(abundFlist)){
    load(abundFlist[y])
    abds[[y]] <- raster::stack(abund)
  }

  abund <- make_sdm_predictions(model = 'ens', rasts = abds, weights = weights, static_variables = NULL, mask = F, bathy_raster = NULL, bathy_max = NULL, se = NULL, month_col = NULL, year_col = NULL, xy_col = NULL)
  nms <- expand.grid(month.abb, yr_min:yr_max)
  names(abund) <- paste(nms$Var1, nms$Var2, sep = '.')
  save(abund, file = paste(file.path(getwd(),spp, 'output_rasters'), '/ENSEMBLE', '_', yr_min, '_', yr_max, '.RData', sep = ''))
  print(Sys.time())

  return(paste(file.path(getwd(),spp, 'output_rasters'), '/ENSEMBLE', '_', yr_min, '_', yr_max, '.RData', sep = ''))
}

