#' @title Create a spreadsheet containing all SDM performance metrics
#' @description
#' Loads existing performance metrics for component models, and calculates Area under the Curve (AUC) for the ensemble model 2 ways. This function also has the option to test the model on an external dataset using \code{testEns}. The final product is a saved csv file containing all of the performance metrics.
#'
#' @param spp.list the data.frame containing species names and alternative names for \code{testEns}. Must contain the column \code{Name}, which matches the names of the species folders to help pull correct data.
#' @param testEns a TRUE/FALSE indicating whether or not to test the ensemble model on an external dataset
#' @param yrMin start of year range to specify which predictions to use. Only used if \code{testEns == TRUE}
#' @param yrMax end of year range to specify which predictions to use. Only used if \code{testEns == TRUE}
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
#'  \item{ENS.AUC}{Ensemble model AUC based on weighted sum of CV predictions from each component model}
#'  \item{WAVG.ENS.AUC}{Ensemble model AUC calculated as a weighted average of AUCs, using the corresponding weights}
#'  \item{AUC.yrMin.yrMax}{Ensemble model AUC calculated on the external data, if \code{testEns == TRUE}. Column name will reflect the timeseries used}
#' }
#'

make_performanceCSV <- function(spp.list, testEns = T, yrMin, yrMax){
#spp.list is the csv of species lists including alternative names to help with matching in testEns

  #subset spp.list to serve as base for csv
  sppEval <- spp.list[,colnames(spp.list) %in% c('Common.Name', 'Managing.Body', 'Feeding.Guild', 'Habitat.Guild')]

  message(paste("Pulling Evaluation Metrics..."))
#create null object to append to sppEval
  sEval <- NULL
  #pull existing metrics calculated and saved in workflow
  for(x in 1:nrow(spp.list)){
    #load in data frame to get the number of presences/absences
    load(paste(file.path(getwd(),spp.list$Name[x]), 'pa_clean.RData', sep = '/')) #load data - dfC
    n.pres <- length(which(dfC$value == 1))
    n.abs <- length(which(dfC$value == 0))

    #pull in other model AUCs
    #load in evaluation metrics
    evalFlist <- dir(file.path(getwd(),spp.list$Name[x], 'model_output', 'eval_metrics'), pattern = '.RData', full.names = T)
    eval <- vector(length = length(evalFlist))
    for(y in 1:length(evalFlist)){
      load(evalFlist[y])
      eval[y] <- ev
    } #eval is a vector of the component model AUCS

    #load in ensemble and calculate raw aucs from weighted sum predictions from component models
    load(paste0(file.path(getwd(),spp.list$Name[x], 'model_output', 'models'), '/ENSEMBLE.RData')) #ens
    #calculate auc
    Pred <- ROCR::prediction(ens$pred, ens$abund)
    Perf <- ROCR::performance(Pred, 'auc')
    auc <- Perf@y.values[[1]]

    #create weighted average auc
    load(paste0(file.path(getwd(),spp.list$Name[x], 'model_output'), '/ensemble_weights.RData')) #weights
    aucW <- weighted.mean(eval, weights)

    #correct length of eval if SDMTMB didn't converge (doing this now because if SDMTMB didn't converge, the length of both weights and eval will be right so we don't need to correct until now)
    if(!any(grepl('SDMTMB', evalFlist))){ #if sdmtmb does not converge
      eval[5] <- NA #add a placeholder to both eval and weights
      weights[5] <- NA
    }

    #put it all together and add names
    eval <- c(n.pres, n.abs, eval, weights, auc, aucW)
    names(eval) <- c('N.PRESENCE', 'N.ABSENCE', 'BRT', 'GAM', "MAXENT", "RF", 'SDMTMB', 'BRT.WT', 'GAM.WT', "MAXENT.WT", "RF.WT", 'SDMTMB.WT', 'ENS.AUC', 'WAVG.ENS.AUC')

    #add to sEval
    sEval <- rbind(sEval, eval)
    print(x)
  } #end x

  #now we test the ensemble on the external dataset, if desired

  if(testEns){
    message(paste("Testing Ensemble on new data..."))
    #create altNames vector for testEns
    altNames <- paste(spp.list$Common.Name, spp.list$COM_NAME, spp.list$Scientific.Name, spp.list$Alternate.Name, spp.list$SCI_NAME, spp.list$SCI_NAME_ALT, spp.list$SCI_NAME_ALT2, sep = ',')

    #predict on external dataset
    evalTest <- vector(length = nrow(sppEval)) #using sppEval as the base
    for(y in 1:nrow(spp.list)){
      ##make source list - we have both rasters and csv files, but for some reason this was easier to do with the csv files. But the raster file list includes all of the data used in the model building (some sources were excluded due to low counts or methods that didn't match the species [i.e. a longline survey would not catch shellfish]). So we use the raster list to subset the csv list appropriately for each species.
      #grab list of standardized csvs
      csvFlist <- dir('./Data/csvs/standardized', pattern = paste(yrMin, yrMax, sep = "_"))
      csvFlist <- gsub('.csv', '', csvFlist) #isolate to names only

      #grab species-specific list of rasters
      rastFlist <- dir(file.path(getwd(),spp.list$Name[y], 'input_rasters'), pattern = paste(yrMin, yrMax, sep = "_"))
      rastFlist <- gsub('.nc', '', rastFlist) #isolate to names only
      i <- grep('combined_rasters', rastFlist) #remove combined raster if present
      if(length(i) != 0){
        rastFlist <- rastFlist[-i]
      }

      csvFlist <- csvFlist[csvFlist %in% rastFlist] #subset csvFlist to just datasets used for species
      csvFlist <- paste0('./Data/csvs/standardized/', csvFlist, '.csv') #recreate full paths
      #test ensemble
      evalTest[y] <- testEns(spp = spp.list$Name[y], sppNames = altNames[y], sources = csvFlist, yrMin = yrMin, yrMax = yrMax)
      print(y)
    }
  } #end testEns

  #combine and save spreadsheet
  if(testEns){ #if testEns is true, add result and rename column
    sppEval <- cbind(sppEval, sEval, evalTest)
    names(sppEval)[ncol(sppEval)] <- paste('AUC', yrMin, yrMax, sep = '.')
  } else { #if testEns is FALSE
    sppEval <- cbind(sppEval, sEval) #just combine everything you made in the first for loop
  }
  utils::write.csv(sppEval, file = 'species_evaluation_metrics.csv', row.names = F)

} #end function
