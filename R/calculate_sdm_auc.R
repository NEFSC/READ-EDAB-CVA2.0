#' @title Calculate Model AUC
#' @description Calculate Area Under the ROC Curve (AUC) for one of the component models or the final ensemble from either cross-validation data or on an external dataset
#'
#' @param model one of the following indicating the desired model to calculate the metric for from: gam, maxent, brt, rf, sdmtmb, or ens
#' @param data_type either 'cv' or 'external'. Determines which methods are used for calculating AUC. If 'cv', then the provided predicted and observed values from a cross-validation are used. If 'external', then predicted values are matched to a provided data.frame of external observations not used in the model construction, which are then used to calculate AUC
#' @param data a data.frame. If \code{data_type} is 'cv', then this is the output from \code{pull_sdm_preds}. If \code{data_type} is 'external', then it is a data.frame of observations not used to train the data
#' @param prediction_rasters spatRaster of predicted model values on desired timeseries. Values will be matched to observations in \code{data}.
#'
#' @return a value representing the AUC for the given model
#'
#'@export

calculate_sdm_auc <- function(model, data_type, data,  prediction_rasters) {

  data <- data[stats::complete.cases(data), ]

  if(data_type == 'cv'){
    #get Pred based on model-specific outputs
    if (model == 'gam' | model == 'maxent') {
      Pred <- ROCR::prediction(data$cvpred, data$abund)
    } #end if gam or maxent

    if (model == 'rf' | model == 'brt' | model == 'ens') {
      Pred <- ROCR::prediction(data$pred, data$abund)
    } #end if rf, brt, or ensemble

    if (model == 'sdmtmb') {
      Pred <- ROCR::prediction(data$pred, preds$value)
    } #end if sdmtmb

    #calculate AUC
    Perf <- ROCR::performance(Pred, 'auc')
    met <- Perf@y.values[[1]]

  } #end if data_type = 'cv'

  if(data_type == 'external'){
    #follow same framework as match_fisheries_environment_data
    #the same methods regardless of model type

    ###building data.frame moved to helper function since this is repeated a bit
    preds <- build_preds_df(observations = data, xy_col = c("grid.lon", "grid.lat"), prediction_rasters = prediction_rasters)

    #get AUC
    Pred <- ROCR::prediction(preds$predicted, preds$pa)
    Perf <- ROCR::performance(Pred, 'auc')
    met <- Perf@y.values[[1]]

  }

  return(met)
}
