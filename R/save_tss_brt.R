#' @title Calculate True Skill Statistic (TSS) for BRT models
#'
#' @description
#' Calculate TSS for BRT k-folds. Used within \code{eval_brt}. Included with CVA2.0 package to ensure functionality.
#'
#' @param truth observations
#' @param predicted predicted values
#'
#' @return a numeric value
#' @source From Camrin Brawn (WHOI): https://zenodo.org/records/7971532.

save_tss_brt <- function(truth, predicted){
  pred <- ROCR::prediction(as.vector(abs(predicted)), as.vector(truth))
  TSS <- ROCR::performance(pred, "sens", "spec")
  TSSvals <- max(unlist(TSS@y.values) + unlist(TSS@x.values) - 1)
  return(TSSvals)
}
