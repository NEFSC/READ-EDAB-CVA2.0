saveTSS <- function(truth, predicted, ...){
  pred <- ROCR::prediction(as.vector(abs(predicted)), as.vector(truth))    
  TSS <- ROCR::performance(pred, "sens", "spec")
  TSSvals <- max(unlist(TSS@y.values) + unlist(TSS@x.values) - 1)
  return(TSSvals)
}