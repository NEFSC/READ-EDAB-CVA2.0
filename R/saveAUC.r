saveAUC <- function(truth, predicted, plot_roc = FALSE, ...){
  pred <- ROCR::prediction(as.vector(abs(predicted)), as.vector(truth))    
  roc <- ROCR::performance(pred,"tpr","fpr")
  auc <- ROCR::performance(pred,"auc")
  if (plot_roc){
    plot(roc, ...)
    abline(a=0, b= 1)
  }
  return(auc@y.values)
}