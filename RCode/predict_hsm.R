###make abundances 

make_predictions <- function(mod, model, rasts, staticData, mask = T, bathy_nm, bathy_max, se = NULL, month_col, year_col, xy_col, weights = NULL){
  #make predictions from model output 
  #mod - model to predict from 
  #model - character string listing the model type
  #rasts - list of rasters to predict to
      #for single model predictions should be a list of raster stacks
      #for ensemble predictions, it should be a length(list) == number of models, where each list contains a list of rasters containing predictions from model; needs to be the same length as weights
  #staticData will be the Rdata file containing all the staticVars - same used to match environmental and species data
  #mask - binary - determines if predictions are limited to a specific bathymetry 
  #bathy_nm - name of bathymetry layer in rasts
  #bathy_max - desired maximum bathymetry 
  #se - spp_env dataset used to build models - only used for rf/sdm models
  #weights - ensemble weights for each model in the ensemble - for ensembles only 
  
  if(mod != 'gam' | mod != 'maxent' | mod != 'rf' | mod != 'brt' | mod != 'sdmtmb' | mod != 'ens'){
    stop('model is not gam, maxent, rf, brt, sdmtmb, or ens')
  }

  #make static vars (month/year) into rasters
  r <- subset(rasts[[1]][[1]], 1)
  rlon<-rlat<-r #copy r to rlon and rlat rasters [1]][1]which will contain the longitude and latitude
  xy<-xyFromCell(r,1:length(r)) #matrix of longitudes (x) and latitudes(y)
  rlon[]<-xy[,1] #raster of longitudes
  rlat[]<-xy[,2] #raster of latitudes
  rMonth <- rYear <- r
  
  hsm <- vector(mode = 'list', length = length(rasts[[1]]))
  
  if(model != 'ens'){
    load(staticData) #staticVars object containing a list of rasters with the same extent as the environmental variables
  }
  
  #slightly different models depending on the model 
  if(model == 'gam'){
    print('Predicting GAM model...')
    for(x in 1:nlayers(rasts[[1]][[1]])){ #all the rasters in rasts have the same number of layers so it doesn't matter which one we call
    
    #may need to change if names aren't always going to be month.year 
    mm.year <- strsplit(names(rasts[[1]][[1]])[x], split = '[.]') #all the rasterbricks in normVars also have the same names so again, doesn't matter which one we call
    
    mm <- as.numeric(gsub('X', '', mm.year[[1]][1])) 
    rMonth[] <- mm
    yr <- as.numeric(mm.year[[1]][2])
    rYear[] <- yr
    
    nStack <- vector(mode = 'list', length = length(rasts))
    for(n in 1:length(rasts)){
      nStack[[n]] <- subset(rasts[[n]][[1]], x)
    }
    nStack <- c(nStack, staticVars)
    nStack <- stack(nStack)
    names(nStack) <- c(names(rasts), names(staticVars))
    crs(nStack) <- crs(rasts[[1]][[1]])
    extent(nStack) <- extent(rasts[[1]][[1]])
    
    sr <- stack(rlon, rlat, rMonth, rYear, nStack)
    names(sr)[1:4] <- c("x", "y", "month", "year")
    
    if(mask){ #if mask == T
      #mask off waters deeper than 1000 m
      i <- which(names(sr) == bathy_nm)
      sr <- replace(sr, abs(subset(sr, i)) > bathy_max, NA) #replace values with an absolute value greater than bathy_max with NA
    }
    
    hsm[[x]] <- MakeGAMAbundance(model = mod, r.stack = sr)
    #print(x)
    }
    names(hsm) <- names(rasts[[1]][[1]])
  } #end if gam
  
  if(model == 'maxent'){
    print('Predicting MAXENT model...')
    for(x in 1:nlayers(rasts[[1]][[1]])){ #all the rasters in rasts have the same number of layers so it doesn't matter which one we call
      
      #may need to change if names aren't always going to be month.year 
      mm.year <- strsplit(names(rasts[[1]][[1]])[x], split = '[.]') #all the rasterbricks in normVars also have the same names so again, doesn't matter which one we call
      
      mm <- as.numeric(gsub('X', '', mm.year[[1]][1])) 
      rMonth[] <- mm
      yr <- as.numeric(mm.year[[1]][2])
      rYear[] <- yr
      
      nStack <- vector(mode = 'list', length = length(rasts))
      for(n in 1:length(rasts)){
        nStack[[n]] <- subset(rasts[[n]][[1]], x)
      }
      nStack <- c(nStack, staticVars)
      nStack <- stack(nStack)
      names(nStack) <- c(names(rasts), names(staticVars))
      crs(nStack) <- crs(rasts[[1]][[1]])
      extent(nStack) <- extent(rasts[[1]][[1]])
      
      sr <- stack(rlon, rlat, rMonth, rYear, nStack)
      names(sr)[1:4] <- c("x", "y", "month", "year")
      
      if(mask){ #if mask == T
        #mask off waters deeper than 1000 m
        i <- which(names(sr) == bathy_nm)
        sr <- replace(sr, abs(subset(sr, i)) > bathy_max, NA) #replace values with an absolute value greater than bathy_max with NA
      }
      
     # require(gbm3)
      hsm[[x]] <- raster(MakeMaxEntAbundance2(model = mod, maxent.stack = sr, type = 'maxnet'))
      #print(x)
    }
    names(hsm) <- names(rasts[[1]][[1]])
  } #end if maxent
  
  if(model == 'rf'){
    print('Predicting RFSI...')
    
    se <- cbind(1:nrow(se), se) #stand in station ids 
    colnames(se)[1] <- "staid"
    se$month.year <- paste(se[,month_col], se[,year_col], sep = '-')
    se$year <- year(my(se$month.year))
    
    #subsample by space-time
    set.seed(2025)
    
    #make regions 
    se$region <- NA
    se$region[which(se[,xy_col[1]] > -70 & se[,xy_col[2]] < 41.5)] <- 'GB' #georges bank
    se$region[which(se[,xy_col[1]] > -71 & se[,xy_col[2]] > 41.5)] <- 'GOM' #gulf of maine
    se$region[which(se[,xy_col[1]] < -70 & se[,xy_col[2]] < 42 & se[,xy_col[2]] > 39.5)] <- 'SNE' #southern new england
    se$region[which(se[,xy_col[2]] < 39.5)] <- 'MAB' #mid-atlantic bight
    
    #make space-time id
    se$sp.tm <- paste(se$month.year, se$region, sep = '-')
    sptm <- unique(se$sp.tm)
    #subsample data 
    seSub <- NULL
    for(x in sptm){
      sub <- se[se$sp.tm == x,]
      
      abs <- sub[sub$value == 0,]
      pres <- sub[sub$value == 1,]
      
      if(nrow(pres) <= 5){ #if there are few presences
        absSub <- abs[sample(x = nrow(abs), size = round(nrow(abs)/4)),] #subsample absences to a minimum number per month and region
        allSub <- rbind(absSub, pres) #combine with presences (if any are absent)
        #this will allow all regions, years, and months to be present in the final time series to help predictions while also making the ratio of presences/absences somewhat more even 
      } else if(nrow(abs) > nrow(pres)){ #if there are enough presences, but absences still outnumber presences
        absSub <- abs[sample(x = nrow(abs), size = nrow(pres)),] #subsample absences 
        allSub <- rbind(absSub, pres)
      } else { #if presences outnumber absences
        allSub <- sub #do nothing and keep it all 
      }
      
      seSub <- rbind(seSub, allSub)
    }
    
    #convert dataframe to spatial object 
    stDF = st_as_sf(seSub, coords = xy_col, crs = 4326, agr = "constant")
    stDF = st_sftime(stDF, time_column_name = year_col)
    
    for(x in 1:nlayers(rasts[[1]][[1]])){ #all the rasters in rasts have the same number of layers so it doesn't matter which one we call
      
      #may need to change if names aren't always going to be month.year 
      mm.year <- strsplit(names(rasts[[1]][[1]])[x], split = '[.]') #all the rasterbricks in normVars also have the same names so again, doesn't matter which one we call
      
      mm <- as.numeric(gsub('X', '', mm.year[[1]][1])) 
      rMonth[] <- mm
      yr <- as.numeric(mm.year[[1]][2])
      rYear[] <- yr
      
      nStack <- vector(mode = 'list', length = length(rasts))
      for(n in 1:length(rasts)){
        nStack[[n]] <- subset(rasts[[n]][[1]], x)
      }
      nStack <- c(nStack, staticVars)
      nStack <- stack(nStack)
      names(nStack) <- c(names(rasts), names(staticVars))
      crs(nStack) <- crs(rasts[[1]][[1]])
      extent(nStack) <- extent(rasts[[1]][[1]])
      
      sr <- stack(rlon, rlat, rMonth, rYear, nStack)
      names(sr)[1:4] <- c("x", "y", "month", "year")
      
      if(mask){ #if mask == T
        #mask off waters deeper than 1000 m
        i <- which(names(sr) == bathy_nm)
        sr <- replace(sr, abs(subset(sr, i)) > bathy_max, NA) #replace values with an absolute value greater than bathy_max with NA
      }
      
      ##convert rasterStack to dataframe to play well with model 
      srDF <- as.data.frame(rasterToPoints(sr))
      #remove NAs
      # srDF <- srDF[-which(srDF$bathy < -1000 | srDF$bathy > 0),]
      #create staid/year 
      srDF$staid <- 1:nrow(srDF) #standin station ids 
      srDF$year <- year(my(paste(srDF$month, srDF$year, sep = '-')))
      
      #convert dataframe to spatial object 
      srDF = st_as_sf(srDF, coords = xy_col, crs = 4326, agr = "constant")
      srDF = st_sftime(srDF, time_column_name = 'year')
      
      predDF <- pred.rfsi(model = mod,
                          data = stDF, 
                          data.staid.x.y.z = c('staid', xy_col, year_col),
                          obs.col = "value",
                          newdata = srDF, 
                          newdata.staid.x.y.z = c('staid', 'X', 'Y', 'year'),
                          output.format = "data.frame", # "sf", # "SpatVector", 
                          cpus = 1, # detectCores()-1,
                          progress = TRUE,
                          classification = F, 
                          no.obs = 'exactly')
      
     raw <- rasterize(x = predDF[,c('X', 'Y')], y = rlon, field = predDF[,'pred'], fun = mean)
     
     #the RF rasters tend to have some mostly small gaps, so using focal to fill those in with nearest neighbor averaging
     w <- matrix(1, 3, 3) 
     hsm[[x]] <- focal(raw, w, mean, na.rm=TRUE, NAonly=TRUE)
      #print(x)
    }
    names(hsm) <- names(rasts[[1]][[1]])
  } #end if rf
  
  if(model == 'brt'){
    print('Predicting Boosted Regression Tree...')
    for(x in 1:nlayers(rasts[[1]][[1]])){ #all the rasters in rasts have the same number of layers so it doesn't matter which one we call
      
      #may need to change if names aren't always going to be month.year 
      mm.year <- strsplit(names(rasts[[1]][[1]])[x], split = '[.]') #all the rasterbricks in normVars also have the same names so again, doesn't matter which one we call
      
      mm <- as.numeric(gsub('X', '', mm.year[[1]][1])) 
      rMonth[] <- mm
      yr <- as.numeric(mm.year[[1]][2])
      rYear[] <- yr
      
      nStack <- vector(mode = 'list', length = length(rasts))
      for(n in 1:length(rasts)){
        nStack[[n]] <- subset(rasts[[n]][[1]], x)
      }
      nStack <- c(nStack, staticVars)
      nStack <- stack(nStack)
      names(nStack) <- c(names(rasts), names(staticVars))
      crs(nStack) <- crs(rasts[[1]][[1]])
      extent(nStack) <- extent(rasts[[1]][[1]])
      
      sr <- stack(rlon, rlat, rMonth, rYear, nStack)
      names(sr)[1:4] <- c("x", "y", "month", "year")
      
      if(mask){ #if mask == T
        #mask off waters deeper than 1000 m
        i <- which(names(sr) == bathy_nm)
        sr <- replace(sr, abs(subset(sr, i)) > bathy_max, NA) #replace values with an absolute value greater than bathy_max with NA
      }
      
     require(gbm)
      srDF <- as.data.frame(rasterToPoints(sr))
      srDF <- srDF[complete.cases(srDF),-c(1:2)]
      hsm[[x]] <- raster::predict(object = sr, model = mod, type="response")
      
     # hsm[[x]] <- raster::rasterize(x = srDF[,1:2], y = rYear, field = p)
      #print(x)
    }
    names(hsm) <- names(rasts[[1]][[1]])
  } #end if brt
  
  if(model == 'sdmtmb'){
    print('Predicting sdmTMB...')
    for(x in 1:nlayers(rasts[[1]][[1]])){ #all the rasters in rasts have the same number of layers so it doesn't matter which one we call
      
      #may need to change if names aren't always going to be month.year 
      mm.year <- strsplit(names(rasts[[1]][[1]])[x], split = '[.]') #all the rasterbricks in normVars also have the same names so again, doesn't matter which one we call
      
      mm <- as.numeric(gsub('X', '', mm.year[[1]][1])) 
      rMonth[] <- mm
      yr <- as.numeric(mm.year[[1]][2])
      rYear[] <- yr
      
      nStack <- vector(mode = 'list', length = length(rasts))
      for(n in 1:length(rasts)){
        nStack[[n]] <- subset(rasts[[n]][[1]], x)
      }
      nStack <- c(nStack, staticVars)
      nStack <- stack(nStack)
      names(nStack) <- c(names(rasts), names(staticVars))
      crs(nStack) <- crs(rasts[[1]][[1]])
      extent(nStack) <- extent(rasts[[1]][[1]])
      
      sr <- stack(rlon, rlat, rMonth, rYear, nStack)
      names(sr)[1:4] <- c("x", "y", "month", "year")
      
      if(mask){ #if mask == T
        #mask off waters deeper than 1000 m
        i <- which(names(sr) == bathy_nm)
        sr <- replace(sr, abs(subset(sr, i)) > bathy_max, NA) #replace values with an absolute value greater than bathy_max with NA
      }
      
      ##convert rasterStack to dataframe to play well with model 
      srDF <- as.data.frame(rasterToPoints(sr)[,-c(1:2)])
      srDF <- srDF[complete.cases(srDF),]

      pred <- predict(mod, newdata = srDF, type = 'response')
      #sdmpred$prob <- exp(sdmpred$est)/(1+exp(sdmpred$est))
      coordinates(pred) <- ~x + y
      proj4string(pred) <- CRS("+proj=longlat +datum=WGS84 +no_defs ")
      hsm[[x]] <- raster::rasterize(x = pred, y = rYear, field = pred$est)
      #print(x)
    }
    names(hsm) <- names(rasts[[1]][[1]])
  } #end if sdmtmb
  
  if(model == 'ens'){
    print('Predicting Ensemble...')
    
    if(length(rasts) == length(weights)){ #make sure they are the same length

          wts <- weights
    
      for(x in 1:length(rasts[[1]])){
        #make list of abundances for timestamp
        abunds <- vector(mode = 'list', length = length(rasts))
        for(m in 1:length(rasts)){
          abunds[[m]] <- terra::rast(rasts[[m]][[x]])
        }
        #make weighted average ensembles 
        hsm[[x]] <- raster(MakeEnsembleAbundance(model.weights = wts, abund.list = abunds))
      } #end x
    } else {
      stop('raster list and weights are not the same length')
    }
  } #end if ensemble

  names(hsm) <- names(rasts[[1]])
  
  return(hsm)
}
