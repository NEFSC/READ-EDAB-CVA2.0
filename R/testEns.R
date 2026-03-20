#' @title Test the Ensemble Model with Observations
#' @description
#' Calculate Area under the Curve (AUC) for the ensemble species distribution model on an external set of data (not used for building the model). This function builds a dataset of presence/absence with the new data across a variety of sources, extracts the predicted values from the ensemble, and calculates the resulting AUC.
#'
#' @param spp species name. Used to help pull correct ensemble model.
#' @param sppNames a vector containing all possible names for the target species. Must have a length >= 1
#' @param sources a vector containing all of the source csvs containing the observation data to use in the test. These are loaded in directly by \code{read.csv} so they should either be full file paths or in the current working directory.
#' @param yrMin start of year range to specify which predictions to use
#' @param yrMax end of year range to specify which predictions to use
#'
#' @return returns the AUC value for the given model on the new data


testEns <- function(spp, sppNames, sources, yrMin, yrMax){

  nms <- strsplit(sppNames, split = ',')[[1]] #seperate out species names

  #combine csvs
  paDF <- NULL
  for(x in sources){
    csv <- read.csv(x)

    ind <- csv$name %in% nms
    csv$pa <- 0
    csv$pa <- replace(csv$pa, ind == TRUE, 1)

    #csv$source <- sources[x]

    #combining the csvs sometimes results in weird column names which then makes it so the csvs can't get combined, so these checks just remove those
    i <- grep("X", names(csv))
    if(length(i) != 0){
      csv <- csv[,-i]
    }
    i <- grep("...1", names(csv))
    if(length(i) != 0){
      csv <- csv[,-i]
    }
    i <- grep("time", names(csv))
    if(length(i) != 0){
      csv <- csv[,-i]
    }
    i <- grep("date", names(csv))
    if(length(i) != 0){
      csv <- csv[,-i]
    }
    paDF <- rbind(paDF, csv)
  }
  paDF$month <- month.abb[paDF$month]
  paDF$month.year <- paste(paDF$month, paDF$year, sep = '.')

  #load in ensemble predictions
  load(paste0(file.path(getwd(),spp, 'output_rasters'), '/ENSEMBLE_', yrMin, '_', yrMax, '.RData'))
  abund <- raster::stack(abund)

  #extract predicted values at observation locations
  paPred <- NULL
  for(x in 1:raster::nlayers(abund)){ #nlayers(abund) should equal number of unique names in paDF
    sub <- paDF[paDF$month.year == names(abund)[x],]
    sub$pred <- raster::extract(raster::subset(abund,x), sub[,c('lon', 'lat')])
    paPred <- rbind(paPred, sub)
    #print(x)
  }
  paPred <- paPred[complete.cases(paPred),]

  #calculate AUC
  Pred <- ROCR::prediction(paPred$pred, paPred$pa)
  Perf <- ROCR::performance(Pred, 'auc')
  auc <- Perf@y.values[[1]]
  return(auc)

}
