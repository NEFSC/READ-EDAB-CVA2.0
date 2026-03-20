#' @title K-fold BRT fit
#' @description k-fold cross validation code from Camrin Brawn (WHOI): https://zenodo.org/records/7971532. Function originally by B. Abrahms from https://github.com/briana-abrahms/DynamicEnsembleSDM/blob/master/model_evaluation.R. Modified for CVA2.0 to return information necessary to calculate evaluation metrics. Follows many of the same inputs as \code{gbm.step}
#'
#' @param dataInput input data frame
#' @param gbm.x names of predictor variables in \code{dataInput}
#' @param gbm.y name of response variable in \code{dataInput}
#' @param learning.rate sets the weight applied to individual trees input to dismo::gbm.step
#' @param k_folds number of folds to perform
#' @param tree.complexity is tree complexity input to dismo::gbm.step - sets complexity of individual trees
#' @param bag.fraction is bag fraction input to dismo::gbm.step - sets the proportion of observations used in selecting variables
#' @param is_fixed TRUE/FALSE to determine if gbm.step or gbm.fixed should be used.
#' @param max.trees maximum number of trees to fit before stopping
#'
#' @return a list containing 1) the output from \code{eval_brt} and 2) a data.frame of the observed and predicted values used to calculate RMSE/AUC by \code{sdm_cv}
#'
#' @examples
#' \dontrun{
#' eval_kfold_BRT(dataInput_Fit, gbm.x=c("curl","ild", "ssh", "sst","sst_sd"), "presabs")
#' }

eval_kfold_brt <- function(dataInput, gbm.x, gbm.y, learning.rate = 0.05, k_folds = 5, tree.complexity = 3, bag.fraction = 0.6, is_fixed = TRUE, max.trees = 2000){
  dataInput$Kset <- dismo::kfold(dataInput, k_folds) #randomly allocate k groups

  summary_stats <- list()
  kRes <- NULL
  for (i in 1:k_folds){
    print(i)
    train <- dataInput[dataInput$Kset!=i,]
    test <- dataInput[dataInput$Kset==i,]
    if (is_fixed){
      brt.k <- dismo::gbm.fixed(data=train, gbm.x= gbm.x, gbm.y = gbm.y,
                                family="bernoulli", tree.complexity = tree.complexity,
                                learning.rate = learning.rate, bag.fraction = bag.fraction, max.trees = max.trees)
    } else{
      brt.k <- dismo::gbm.step(data=train, gbm.x= gbm.x, gbm.y = gbm.y,
                               family="bernoulli", tree.complexity = tree.complexity,
                               learning.rate = learning.rate, bag.fraction = bag.fraction, max.trees = max.trees)
    }
    preds <- gbm::predict.gbm(brt.k, test,
                              n.trees=brt.k$gbm.call$best.trees, type="response")

    summary_stats[[i]] <- as.data.frame(eval_brt(brt.k, test, response = gbm.y, plot = FALSE))

    ###add training values and predictions to calculate RMSE - written by KLG & added 6/20/2025
    res <- cbind(test, preds)
    kRes <- rbind(kRes, res)

  }
  return(list(summary_stats, kRes))
}
