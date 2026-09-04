#' @title Create a spreadsheet containing all SDM performance metrics
#' @description
#' Compiles sample sizes and performance metrics for the component and ensemble modelsThe final product is a saved csv file containing all of the performance metrics.
#'
#' @param spp_list the data.frame containing species names and alternative names for \code{test_ens}. Must contain the column \code{Name}, which matches the names of the species folders to help pull correct data.
#' @param training_years vector with length equal to 2, indicating the maximum and minimum years that identify the desired training datasets
#' @param pa_col column name for presence/absence column
#' @param release release code for MOM6 data. Helps pull correct training dataset associated with the MOM6 data with the same name
#' @param spatial_temporal TRUE/FALSE to determine method for normalizing. Helps pull correct training dataset associated with the MOM6 data with the same name
#' @param mask_bathy TRUE/FALSE indicating whether or not bathymetry data was used as a mask for raw data before normalization. Helps pull correct training dataset associated with the MOM6 data with the same name
#' @param rm_corr TRUE/FALSE indicating whether or not correlated environmental covariates were removed. Helps to pull correct training dataframe
#' @param add_data TRUE/FALSE indicating whether or not to add additional data to the dataset, for example, AUCs calculated on different time periods
#' @param additional_data a data.frame of additional data to add to the data.frame. column names should be desired column names 
#'
#' @return nothing is returned. The resulting CSV file called 'species_evaluation_metrics.csv' is saved in the working directory.
#'
#' @details
#' The resulting CSV file contains the following columns:
#' \describe{
#'  \item{Common.Name}{name of species}
#'  \item{Managing.Body}{name of group responsible for species management, if exists}
#'  \item{Feeding.Guild}{assigned feeding guild}
#'  \item{Habitat.Guild}{assigned habitat guild}
#'  \item{N.PRESENCE, N.ABSENCE}{number of presences/absences in the data set used to build the model}
#'  \item{BRT, GAM, MAXENT, RF, SDMTMB}{AUC for each component model}
#'  \item{BRT.WT, GAM.WT, MAXENT.WT, RF.WT, SDMTMB.WT}{weights for each component model in the final ensemble}
#'  \item{ENS.AUC}{Ensemble model AUC}
#'  \item{...}{Any additional columns from \code{add_data}}
#' }
#'
#'
#'@export

make_evaluation_csv <- function(spp_list, training_years, pa_col, release, spatial_temporal, mask_bathy, rm_corr, add_data, additional_data) {
  #suffixes to help locate correct data
  suffix <- if(spatial_temporal) "" else "_global"
  bathy_suffix <- if(mask_bathy) "masked" else ""
  corr_suffix <- if(rm_corr) "rmcorr" else ""
  
  #subset spp_list to serve as base for csv
  sppEval <- spp_list[,
    colnames(spp_list) %in%
      c('Common.Name', 'Managing.Body', 'Feeding.Guild', 'Habitat.Guild')
  ]

  message(paste("Pulling Evaluation Metrics..."))
  #create null object to append to sppEval
  sEval <- NULL
  #pull existing metrics calculated and saved in workflow
  for (x in 1:nrow(spp_list)) {
    
    #load in data frame to get the number of presences/absences
    training_name <- file.path(getwd(), spp_list$Name[x], paste0('training_', training_years[1], '_', training_years[2], '_', corr_suffix, '_hindcast_', release, '_', bathy_suffix, suffix, '.csv'))
    
    if (!file.exists(training_name)) {
      stop("Aborting: training dataset not found.")
    }
    dfT <- read.csv(file.path(training_name))
    
    n.pres <- length(which(dfT[,pa_col] == 1))
    n.abs <- length(which(dfT[,pa_col] == 0))

    #pull in other model AUCs **need to check order of these **
    #load in evaluation metrics
    evalFlist <- dir(
      file.path(getwd(), spp_list$Name[x], 'model_output', 'eval_metrics'),
      pattern = '.rds',
      full.names = T
    )
    
    #reorder evalFlist to put ensemble last
    evalFlist <- evalFlist[c(1,3:6,2)]
    
    eval <- vector(length = length(evalFlist))
    for (y in 1:length(evalFlist)) {
      if(!is.na(evalFlist[y])){
        load(evalFlist[y])
        eval[y] <- ev
      } else {
        eval[y] <- NA
      }
    } #eval is a vector of the component model + ensemble AUCs

   load(file.path(getwd(), spp_list$Name[x], 'model_output', 'ensemble_weights.rds')) #weights
   if(length(weights) < 5){
     weights <- c(weights, NA) #if weights only has 4 mondels, append an NA
   }

    #put it all together and add names
    eval <- c(n.pres, n.abs, eval, weights)
    names(eval) <- c(
      'N.PRESENCE',
      'N.ABSENCE',
      'BRT.AUC',
      'GAM.AUC',
      "MAXENT.AUC",
      "RF.AUC",
      'SDMTMB.AUC',
      'ENS.AUC',
      'BRT.WT',
      'GAM.WT',
      "MAXENT.WT",
      "RF.WT",
      'SDMTMB.WT'
    )

    #add to sEval
    sEval <- rbind(sEval, eval)
  } #end x
  
  sppEval <- cbind(sppEval, sEval)

  # an option to add ational data such as other stats
  if (add_data) {
    message(paste("Adding additional data..."))
    sppEval <- cbind(sppEval, additional_data) #combining everything
  }
  utils::write.csv(
    sppEval,
    file = 'species_evaluation_metrics.csv',
    row.names = F
  )
} #end function
