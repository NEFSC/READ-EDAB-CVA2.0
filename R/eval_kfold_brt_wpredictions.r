## K-fold BRT fit
## Fun by B. Abrahms from https://github.com/briana-abrahms/DynamicEnsembleSDM/blob/master/model_evaluation.R
#' same inputs as gbm.step right?
#' @param tc is tree complexity input to dismo::gbm.step
#' @param bf is bag fraction input to dismo::gbm.step
#' @example 
#' eval_kfold_BRT(dataInput_Fit, gbm.x=c("curl","ild", "ssh", "sst","sst_sd", "ssh_sd", "z", "z_sd", "EKE","slope","aspect","BV"), "presabs")

eval_kfold_brt <- function(dataInput, gbm.x, gbm.y, learning.rate = 0.05, k_folds = 5, tree.complexity = 3, bag.fraction = 0.6, is_fixed = TRUE){
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
                             learning.rate = learning.rate, bag.fraction = bag.fraction)
    } else{
      brt.k <- dismo::gbm.step(data=train, gbm.x= gbm.x, gbm.y = gbm.y, 
                               family="bernoulli", tree.complexity = tree.complexity,
                               learning.rate = learning.rate, bag.fraction = bag.fraction)
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
