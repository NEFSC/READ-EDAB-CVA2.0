#' @title Calculate BRT model evaluation statistics
#'
#' @description
#' Calculate evaluation statistics for BRT models. From Camrin Brawn (WHOI): https://zenodo.org/records/7971532. Included to ensure functionality.
#'
#' @param model a BRT model object
#' @param test_data is data used to fit the model. Must contain a column name that matches response.
#' @param response is character providing the name of the response variable in test_data
#' @param plot TRUE/FALSE to turn on/off plotting
#'
#' @return model evaluation statistics
#' @source Much of this code is derived from https://github.com/elhazen/PA-paper

eval_brt <- function(model, test_data, response, plot = TRUE) {

  if (!(response %in% names(test_data))) stop('Column name specified in input response object is not available in input data.')

  summarystats <- list()

  pred <-
    gbm::predict.gbm(model,
                     test_data,
                     n.trees = model$gbm.call$best.trees,
                     type = 'response',
                     na.rm = FALSE)

  summarystats$R_squared <- pseudoR2.brt(model)
  summarystats$AUCval <- unlist(saveAUC(test_data[,response], pred))
  summarystats$TSSval <- saveTSS(test_data[,response], pred)
  cm <- caret::confusionMatrix(factor((test_data[,response])), factor(round(pred)))
  ## proportion of true negatives that are correctly predicted
  summarystats$Specificity <- cm$byClass['Specificity']
  ## proportion of true positives that are correctly predicted
  summarystats$Sensitivity <- cm$byClass['Sensitivity']

  summarystats$FalsePos <- mean(pred[test_data[,response] == 0])
  summarystats$FalseNeg <- mean(pred[test_data[,response] == 1])
  summarystats$Accuracy <- unname(cm$overall['Accuracy'])

  summarystats$median_pred_at_pres <- median(gbm::predict.gbm(model, test_data[which(test_data[,response] == 1),], n.trees = model$gbm.call$best.trees, type="response"))
  summarystats$median_pred_at_abs <- median(gbm::predict.gbm(model, test_data[which(test_data[,response] == 0),], n.trees = model$gbm.call$best.trees, type="response"))

  if (plot){
    #pdf(paste("BRT_ROCR_",listnames[i],".pdf",sep=''))
    PredictABEL::plotROC(test_data[,response], pred, colorize = TRUE)
    #dev.off()
  }

  bhvals <- bhattacharyya_stat(test_data, response, model$var.names)
  summarystats$MeanBH <- mean(unlist(bhvals))
  #if (i==1) bhvect<-bhvals else bhvect<-rbind(bhvect,bhvals)

  return(summarystats)
}
