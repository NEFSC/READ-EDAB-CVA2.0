#' @title K-fold BRT fit
#' @description k-fold cross validation code from Camrin Brawn (WHOI): https://zenodo.org/records/7971532. Function originally by B. Abrahms from https://github.com/briana-abrahms/DynamicEnsembleSDM/blob/master/model_evaluation.R. Modified to return information necessary to calculate evaluation metrics. Follows many of the same inputs as \code{gbm.step}
#'
#' @param data_input input data frame
#' @param gbm_x names of predictor variables in \code{data_input}
#' @param gbm_y name of response variable in \code{data_input}
#' @param learning_rate sets the weight applied to individual trees input to dismo::gbm.step
#' @param k_folds number of folds to perform
#' @param tree_complexity is tree complexity input to dismo::gbm.step - sets complexity of individual trees
#' @param bag_fraction is bag fraction input to dismo::gbm.step - sets the proportion of observations used in selecting variables
#' @param is_fixed TRUE/FALSE to determine if gbm.step or gbm.fixed should be used.
#' @param max_trees maximum number of trees to fit before stopping
#'
#' @return a list containing 1) the output from \code{eval_brt} and 2) a data.frame of the observed and predicted values used to calculate RMSE/AUC by \code{sdm_cv}
#'
#' @examples
#' \dontrun{
#' eval_kfold_brt(data_input_Fit, gbm_x=c("curl","ild", "ssh", "sst","sst_sd"), "presabs")
#' }

eval_kfold_brt <- function(data_input, gbm_x, gbm_y, learning_rate = 0.05, k_folds = 5, tree_complexity = 3, bag_fraction = 0.6, is_fixed = TRUE, max_trees = 2000){
  data_input$Kset <- dismo::kfold(data_input, k_folds) #randomly allocate k groups

  summary_stats <- list()
  kRes <- NULL
  for (i in 1:k_folds){
    print(i)
    train <- data_input[data_input$Kset!=i,]
    test <- data_input[data_input$Kset==i,]
    if (is_fixed){
      brt.k <- dismo::gbm.fixed(data=train, gbm.x= gbm_x, gbm.y = gbm_y,
                                family="bernoulli", tree.complexity = tree_complexity,
                                learning.rate = learning_rate, bag.fraction = bag_fraction, max.trees = max_trees)
    } else{
      brt.k <- dismo::gbm.step(data=train, gbm.x= gbm_x, gbm.y = gbm_y,
                               family="bernoulli", tree.complexity = tree_complexity,
                               learning.rate = learning_rate, bag.fraction = bag_fraction, max.trees = max_trees)
    }
    preds <- gbm::predict.gbm(brt.k, test,
                              n.trees=brt.k$gbm.call$best.trees, type="response")

    summary_stats[[i]] <- as.data.frame(eval_brt(brt.k, test, response = gbm_y, plot = FALSE))

    ###add training values and predictions to calculate RMSE - written by KLG & added 6/20/2025
    res <- cbind(test, preds)
    kRes <- rbind(kRes, res)

  }
  return(list(summary_stats, kRes))
}
