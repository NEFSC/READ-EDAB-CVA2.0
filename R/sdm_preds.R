#' @title Extract Predicted Values from Component Model CV
#' @description Extract the predictions from the cross-validation of one of the five component models
#'
#' @param cv output from \code{sdm_cv}
#' @param model one of the following indicating the desired model to extract predicted values from: gam, maxent, brt, rf, or sdmtmb
#'
#' @return a data.frame containing the prediction outputs from the cross-validation necessary to calculate evaluation metric


sdm_preds <- function(cv, model){

  if(model == 'gam' | model == 'maxent'){ #get preds for gam or maxent
    #get evaluation metrics (RMSE & AUC)
    preds <-cv[[1]]
  } #end if gam/maxent

  if(model == 'rf'){
    preds <- cv #the output from cv is also preds here
    colnames(preds)[6:7] <- c('abund', 'pred')
  } #end if RF

  if(model == 'brt'){
    ##put obs and predicted together for ensemble
    preds <- data.frame(abund = cv$value, pred = cv$preds)
    #names(preds) <- c('abund', 'pred')
  } #end if BRT

  if(model == 'sdmtmb'){
    #change names from sdmTMB preds to make it work with EFHSDM
    #preds <- cv$data
    preds <- cv
    colnames(preds)[17] <- c('pred')
  } #end if ensemble

  return(preds)

}
