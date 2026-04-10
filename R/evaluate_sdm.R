#' @title Calculate Performance Metric
#' @description Calculate RMSE or AUC for one of the component models or the final ensemble
#'
#' @param preds output from \code{pull_sdm_preds}
#' @param model one of the following indicating the desired model to calculate the metric for from: gam, maxent, brt, rf, sdmtmb, or ens
#' @param metric either rmse or auc - determines which metric is calculated
#' @param mod the output from \code{build_sdm} - only used for ensemble
#'
#' @return the desired evaluation metric for the given model

evaluate_sdm <- function(preds, model, metric, mod = NULL) {
  if (model == 'gam' | model == 'maxent') {
    #get metrics for gam or maxent model
    preds2 <- preds[complete.cases(preds), ]
    if (metric == 'rmse') {
      #RMSE
      met <- EFHSDM::RMSE(obs = preds2$abund, pred = preds2$cvpred)
    }

    if (metric == 'auc') {
      #AUC
      Pred <- ROCR::prediction(preds2$cvpred, preds2$abund)
      Perf <- ROCR::performance(Pred, 'auc')
      met <- Perf@y.values[[1]]
    }
  } #end if gam

  if (model == 'rf') {
    if (metric == 'rmse') {
      #get RMSE
      met <- EFHSDM::RMSE(obs = preds$abund, pred = preds$pred)
    }

    if (metric == 'auc') {
      Pred <- ROCR::prediction(preds$pred, preds$abund)
      Perf <- ROCR::performance(Pred, 'auc')
      met <- Perf@y.values[[1]]
    }
  } #end if RF

  if (model == 'brt') {
    if (metric == 'rmse') {
      met <- EFHSDM::RMSE(obs = preds$abund, pred = preds$pred)
    }

    if (metric == 'auc') {
      #calculate AUC
      Pred <- ROCR::prediction(preds$pred, preds$abund)
      Perf <- ROCR::performance(Pred, 'auc')
      met <- Perf@y.values[[1]]
    }
  } #end if BRT

  if (model == 'sdmtmb') {
    if (metric == 'rmse') {
      #calculate RMSE
      met <- EFHSDM::RMSE(obs = preds$value, pred = preds$pred)
    }

    if (metric == 'auc') {
      #calculate AUC
      Pred <- ROCR::prediction(preds$pred, preds$value)
      Perf <- ROCR::performance(Pred, 'auc')
      met <- Perf@y.values[[1]]
    }
  } #end if sdmtmb

  if (model == 'ens') {
    if (metric == 'rmse') {
      #RMSE
      met <- EFHSDM::RMSE(obs = mod$abund, pred = mod$pred)
    }

    if (metric == 'auc') {
      Pred <- ROCR::prediction(mod$pred, mod$abund)
      Perf <- ROCR::performance(Pred, 'auc')
      met <- Perf@y.values[[1]]
    }
  } #end if ensemble

  return(met)
}
