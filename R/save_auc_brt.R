#' @title Calculate Area under the Curve (AUC) for BRT models
#'
#' @description
#' Calculate AUC for BRT k-folds. Used within \code{eval_brt}. Included with CVA2.0 package to ensure functionality.
#'
#' @param truth observations
#' @param predicted predicted values
#' @param plot_roc TRUE/FALSE to turn on/off plotting
#' @param ... arguments passed to \code{ROCR::plot}
#'
#' @return an objected of class \code{performance}. Same as \code{ROCR::performance}.
#' @source From Camrin Brawn (WHOI): https://zenodo.org/records/7971532.

save_auc_brt <- function(truth, predicted, plot_roc = FALSE, ...) {
  pred <- ROCR::prediction(as.vector(abs(predicted)), as.vector(truth))
  roc <- ROCR::performance(pred, "tpr", "fpr")
  auc <- ROCR::performance(pred, "auc")
  if (plot_roc) {
    plot(roc, ...)
    abline(a = 0, b = 1)
  }
  return(auc@y.values)
}
