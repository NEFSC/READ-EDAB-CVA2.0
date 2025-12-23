#' @title Exposure Functions
#' @description Functions to calculate raw and species-specific variable and total exposures, using exposure calculation, ranking methods, and logic rule to combine variable-specific exposure from the CVA1.0, while integrating the SDM results

#' \itemize{
#' \item \code{calcExposure} calculates variable exposure as a ratio of the difference between future and present conditions, divided by present standard deviations (a z-score). 
#' \item \code{rankExposure} sets exposure to ranks ranging from 1 - 4 based on thresholds used in CVA1.0. This function also flips variables as needed such that positive exposure is always negative (i.e. both increasing temperature and decreasing oxygen concentrations result in positive exposure)
#' \item \code{makeExposureTimeseries} and \code{makeExposureMaps} combine raw variable exposures and SDM results, using the SDM results as weights to average the exposure across time and space to produce timeseries and maps of variable exposures
#' \item \code{combineWeights} normalizes the component model variable importances, as well as the model weights in the final ensemble, to create a weighted average vector of variable importance for the ensemble model. Options allow for static variables, such as time and space, or non-dynamic environmental variables such as bathymetry, to be excluded from the normalization and weighted average. 
#' \item \code{combineTimeseries} and \code{combineMaps} uses the logic rule from CVA1.0 to combine variable-specific exposures to create a total exposure timeseries and map 
}

#' @param presentRasts,futureRasts lists of rasterStacks, representing the present and future timeseries to calculate exposure with, with each item in the list corresponding to an environmental variable. Both lists need to have the same length (aka number of variables) but the lengths of the timeseries (equal to the number of layers in rasterStacks) can differ. Output of \code{norm_env}
#' @param exposure output from \code{calcExposure}
#' @param flip TRUE/FALSE option to multiply exposure by -1 to ensure that positive exposure values represent negative variable change (increasing temperatures, decreasing oxygen concentrations, etc). 
#' @param noflipList character vector naming which vectors should NOT be flipped - counterintuitive, yes, but this list was actually shorter than the variables that needed to be flipped. Names should match the names in \code{exposure}
#' @param rankExp output from \code{rankExposure}
#' @param sdmRast rasterStack with monthly averages of SDM results to use as weights for weighted average
#' @param vars character vector of variables to generate normalized and averaged ensemble importance for 
#' @param ensWeights vector of component model weights in the ensemble
#' @param impFlist character vector of file paths to variable importance for component models
#' @param matExp,mapExp matrix or rasterStack output from \code{makeExposureTimeseries} and \code{makeExposureMaps}
#' @param countAll TRUE/FALSE to use \code{weights} and \code{wThreshold} to subset variables to only important variables
#' @param weights output from \code{combineWeights} - a vector of variable weights in ensemble SDM
#' @param wThreshold numeric value used to subset weights, variables with weights less to or equal to this value will be excluded from total exposure calculation

#' @return \code{calcExposure} and \code{rankExposure} return a list of rasterStacks with each layer representing monthly raw exposure. \code{calcExposure} returns rasterStacks with the raw exposure (z-score) values and \code{rankExposure} returns ranked values between 1 - 4. The length of the list is equal to the length of the lists supplied as \code{presentRasts} and \code{futureRasts} 
#' @return \code{makeExposureTimeseries} and \code{makeExposureMaps} returns weighted average matrix (\code{makeExposureTimeseries}) and rasterStack (\code{makeExposureMaps}) of exposure, weighted by SDM results. \code{makeExposureTimeseries} returns a matrix with the number of rows equal to the number of variables supplied, and 12 columns (1 for each month). \code{makeExposureMaps} returns a rasterStack with the number of layers equal to the number of variables supplied. 
#' @return \code{combineWeights} returns a vector representing the weighted average of normalized variable importance, representing variable importance in the final ensemble SDM. 
#' @return \code{combineTimeseries} and \code{combineMaps} returns a single vector (\code{combineTimeseries}) or raster (\code{combineMaps}) representing total exposure across time or space

calcExposure <- function(presentRasts, futureRasts){
  EXP <- vector(mode = 'list', length = length(futureRasts))
  for(v in 1:length(futureRasts)){
    ex <- vector(mode = 'list', length = 12)
    for(x in 1:12){
      #take mean of 'present' and 'future'
      mPres <- calc(raster::subset(presentRasts[[v]][[1]], seq(x, nlayers(presentRasts[[v]][[1]]), by = 12)), mean)
      mFut <- calc(raster::subset(futureRasts[[v]][[1]], seq(x, nlayers(futureRasts[[v]][[1]]), by = 12)), mean)
      
      #difference in means
      mDiff <- mFut - mPres
      
      ##calculate SD 
      sdPres <- calc(raster::subset(rawP[[v]][[1]], seq(x, nlayers(rawP[[v]][[1]]), by = 12)), sd)
      sdFut <- calc(raster::subset(futureRasts[[v]][[1]], seq(x, nlayers(futureRasts[[v]][[1]]), by = 12)), sd)
      
      #calculate exposure
      ex[[x]] <- stack(mDiff) / stack(sdPres)
    }
    
    #stack  
    EXP[[v]] <- stack(ex)
    print(v)
  }
  names(EXP) <- names(futureRasts)
  return(EXP)
}

rankExposure <- function(exposure, flip = T, noflipList = c('bottomT', 'surfaceT', 'bottomArg', 'MLD')){
  rankE <- vector(mode = 'list', length = length(exposure)) #list of vectors for 
  for(x in 1:length(exposure)){
    s <- exposure[[x]] 
    
    #change sign - negative = exposure to worse habitat? 
    if(flip){
      if(names(exposure)[x] %in% noflipList){
        s <- s 
      } else {
        s <- s*-1
      }
    }
    
    mapE <- vector(mode = 'list', length = 12) #list for maps
    for(m in 1:12){
      r <- raster::subset(s, m)
      #rank
      QR <- replace(r, r <= 0.5, 1)
      QR <- replace(QR, r > 0.5 & r <= 1.5, 2)
      QR <- replace(QR, r > 1.5 & r <= 2, 3)
      QR <- replace(QR, r > 2, 4)
      
      mapE[[m]] <- QR
    }

    rankE[[x]] <- stack(mapE)
    
  }
  names(rankE) <- names(exposure)
  return(rankE)
}

makeExposureTimeseries <- function(rankExp, sdmRast){
  matExp <- matrix(nrow = length(rankExp), ncol = 12)
  for(x in 1:length(rankExp)){
    s <- rankExp[[x]] 
    
    meanExp <- vector(length = 12) #vector
    for(m in 1:12){
      r <- raster::subset(s, m)
      h <- subset(sdmRast, m)
      r[is.na(h)] <- NA
      
      rDF <- rasterToPoints(r)
      hDF <- rasterToPoints(h)
      
      rh <- merge(hDF, rDF, by = c('x', 'y'), all.x = T)
      meanExp[m] <- weighted.mean(rh[,4], w = rh[,3], na.rm = T)
    }
    matExp[x,] <- meanExp
    #print(x)
  }
  rownames(matExp) <- names(rankExp)
  colnames(matExp) <- month.abb
  
  return(matExp)
}

makeExposureMaps <- function(rankExp, sdmRast){
  mapExp <- vector(mode = 'list', length = length(rankExp))
  for(x in 1:length(rankExp)){
    s <- rankExp[[x]] 
    
    #mapExp[[x]] <- raster::calc(stack(mapE), mean, na.rm = T)
    mapExp[[x]] <- weighted.mean(stack(s), w = sdmRast, na.rm = T)
   # print(x)
  }
  names(mapExp) <- names(rankExp)
  return(stack(mapExp))
}

combineWeights <- function(vars, ensWeights, impFlist){
  #set up data frame
  dfI <- as.data.frame(matrix(nrow = length(impFlist), ncol = length(vars)))
  colnames(dfI) <- vars
  
  #pull importance vectors and add to data frame
  for(x in 1:length(impFlist)){
    load(impFlist[x]) #imp
    if(class(imp) == 'data.frame'){
      imp.vec <- imp$rel.inf ### need to remove spatial variables
      names(imp.vec) <- imp$var
      for(y in 1:length(names(imp.vec))){
        if(names(imp.vec)[y] %in% vars){
          i <- which(vars == names(imp.vec)[y])
          dfI[x,i] <- imp.vec[y]
        }
      }
    } else {
     # imp <- range01(imp)
      for(y in 1:length(names(imp))){
        if(names(imp)[y] %in% vars){
          i <- which(vars == names(imp)[y])
          dfI[x,i] <- imp[y]
        }
      }
    }
  }
  
  #we don't have month and year or the static vars so remove those
  dfV <- replace(dfV, is.na(dfV), 0)
  dfV <- t(apply(dfV, 1, FUN = function(x){x/sum(x,na.rm = T)}))
  dfV <- replace(dfV, is.na(dfV), 0)
  ws <- apply(dfV, MARGIN = 2, FUN = weighted.mean, w = weights, na.rm = T) #weighted average of weights 
  
  return(ws)
}

combineTimeseries <- function(matExp, countAll, weights, wThreshold){
  if(countAll){ #if counting all included factors and not taking weight into account
    matSub <- matExp
  } else {
    wi <- which(weights > wThreshold)
    matSub <- matExp[wi,]
  }
  
  #count each rank in each column 
  hr <- hh <- md <- vector(length = ncol(matSub))
  for(x in 1:ncol(matSub)){
    hr[x] <- length(which(matSub[,x] >= 3.5))
    hh[x] <- length(which(matSub[,x] >= 3))
    md[x] <- length(which(matSub[,x] >= 2.5))
  }
  
  #apply logic rule
  expV <- rep(1, times = ncol(matSub))
  expV <- replace(expV,  md >= 2, 2)
  expV <- replace(expV,  hh >= 2, 3)
  expV <- replace(expV,  hr >= 3, 4)
  
  return(expV)
}

combineMaps <- function(mapExp, countAll, weights, wThreshold){
 if(countAll){
   mapSub <- mapExp
 } else {
  wi <- which(weights >= wThreshold)
  mapSub <- subset(mapExp,wi)
 }

  #multiply mean exposure maps by variable weights to scale, and count how many layers have each rank within each cell 
  hr <- calc(mapSub, function(x){length(x[x >= 3.5])}) 
  hh <- calc(mapSub, function(x){length(x[x >= 3])})
  md <- calc(mapSub, function(x){length(x[x >= 2.5])})
  
  hr[is.na(subset(mapSub,1))] <- NA
  hh[is.na(subset(mapSub,1))] <- NA
  md[is.na(subset(mapSub,1))] <- NA
  
  #apply logic rule from Hare et al 2015
  expL <- subset(mapSub,1)
  expL[] <- 1
  expL <- replace(expL,  md >= 2, 2)
  expL <- replace(expL,  hh >= 2, 3)
  expL <- replace(expL,  hr >= 3, 4)
  expL[is.na(subset(mapSub,1))] <- NA
  
  return(expL)
}

