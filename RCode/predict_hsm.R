###make abundances 

make_predictions <- function(mod, model, rasts, mask = T, bathy_nm, bathy_max, se, weights ){
  #make predictions from model output 
  #mod - model to predict from 
  #model - character string listing the model type
  #rasts - list of rasters to predict to
      #for single model predictions should be a list of raster stacks
      #for ensemble predictions, it should be a length(list) == number of models, where each list contains a list of rasters containing predictions from model; needs to be the same length as weights
  #mask - binary - determines if predictions are limited to a specific bathymetry 
  #bathy_nm - name of bathymetry layer in rasts
  #bathy_max - desired maximum bathymetry 
  #se - spp_env dataset used to build models - only used for sdm models
  #weights - ensemble weights for each model in the ensemble - for ensembles only 
  
  #make static vars (month/year) into rasters
  r <- subset(rasts[[1]], 1)
  rlon<-rlat<-r #copy r to rlon and rlat rasters [1]][1]which will contain the longitude and latitude
  xy<-xyFromCell(r,1:length(r)) #matrix of longitudes (x) and latitudes(y)
  rlon[]<-xy[,1] #raster of longitudes
  rlat[]<-xy[,2] #raster of latitides
  rMonth <- rYear <- r
  
  if(model != 'ens'){
    hsm <- vector(mode = 'list', length = nlayers(rasts[[1]]))
  } else {
    hsm <- vector(mode = 'list', length = length(rasts[[1]]))
  }
  
  #slightly different models depending on the model 
  if(model == 'gam'){
    print('Predicting GAM model...')
    for(x in 1:nlayers(rasts[[1]])){ #all the rasters in rasts have the same number of layers so it doesn't matter which one we call
    
    #may need to change if names aren't always going to be month.year 
    mm.year <- strsplit(names(rasts[[1]])[x], split = '[.]') #all the rasterbricks in normVars also have the same names so again, doesn't matter which one we call
    
    mm <- match(mm.year[[1]][1],month.abb) 
    rMonth[] <- mm
    yr <- as.numeric(mm.year[[1]][2])
    rYear[] <- yr
    
    nStack <- vector(mode = 'list', length = length(rasts))
    for(n in 1:length(rasts)){
      nStack[[n]] <- subset(rasts[[n]], x)
    }
    nStack <- stack(nStack)
    crs(nStack) <- crs(rasts[[1]])
    extent(nStack) <- extent(rasts[[1]])
    
    sr <- stack(rlon, rlat, rMonth, rYear, nStack)
    names(sr) <- c("x", "y", "month_num", "year", names(rasts))
    
    if(mask){ #if mask == T
      #mask off waters deeper than 1000 m
      i <- which(names(sr) == bathy_nm)
      sr <- replace(sr, abs(subset(sr, i)) > bathy_max, NA) #replace values with an absolute value greater than bathy_max with NA
    }
    
    hsm[[x]] <- MakeGAMAbundance(model = mod, r.stack = sr)
    #print(x)
    }
  } #end if gam
  
  if(model == 'max'){
    print('Predicting MAXENT model...')
    for(x in 1:nlayers(rasts[[1]])){ #all the rasters in rasts have the same number of layers so it doesn't matter which one we call
      
      #may need to change if names aren't always going to be month.year 
      mm.year <- strsplit(names(rasts[[1]])[x], split = '[.]') #all the rasterbricks in normVars also have the same names so again, doesn't matter which one we call
      
      mm <- match(mm.year[[1]][1],month.abb) 
      rMonth[] <- mm
      yr <- as.numeric(mm.year[[1]][2])
      rYear[] <- yr
      
      nStack <- vector(mode = 'list', length = length(rasts))
      for(n in 1:length(rasts)){
        nStack[[n]] <- subset(rasts[[n]], x)
      }
      nStack <- stack(nStack)
      crs(nStack) <- crs(rasts[[1]])
      extent(nStack) <- extent(rasts[[1]])
      
      sr <- stack(rlon, rlat, rMonth, rYear, nStack)
      names(sr) <- c("x", "y", "month_num", "year", names(rasts))
      
      if(mask){ #if mask == T
        #mask off waters deeper than 1000 m
        i <- which(names(sr) == bathy_nm)
        sr <- replace(sr, abs(subset(sr, i)) > bathy_max, NA) #replace values with an absolute value greater than bathy_max with NA
      }
      
      hsm[[x]] <- raster(MakeMaxEntAbundance(model = mod, maxent.stack = sr, type = 'maxnet'))
      #print(x)
    }
  } #end if maxent
  
  if(model == 'rf'){
    print('Predicting RFSI...')
    
    ##need to make old dataframe for prediction 
    se <- cbind(1:nrow(se), se) #stand in station ids 
    colnames(se)[1] <- "staid"
    se$year <- year(my(paste(se[,month_col], se[,year_col], sep = '-')))
    
    #convert dataframe to spatial object 
    stDF = st_as_sf(se, coords = xy_col, crs = 4326, agr = "constant")
    stDF = st_sftime(stDF, time_column_name = year_col)
    
    
    for(x in 1:nlayers(rasts[[1]])){ #all the rasters in rasts have the same number of layers so it doesn't matter which one we call
      
      #may need to change if names aren't always going to be month.year 
      mm.year <- strsplit(names(rasts[[1]])[x], split = '[.]') #all the rasterbricks in normVars also have the same names so again, doesn't matter which one we call
      
      mm <- match(mm.year[[1]][1],month.abb) 
      rMonth[] <- mm
      yr <- as.numeric(mm.year[[1]][2])
      rYear[] <- yr
      
      nStack <- vector(mode = 'list', length = length(rasts))
      for(n in 1:length(rasts)){
        nStack[[n]] <- subset(rasts[[n]], x)
      }
      nStack <- stack(nStack)
      crs(nStack) <- crs(rasts[[1]])
      extent(nStack) <- extent(rasts[[1]])
      
      sr <- stack(rlon, rlat, rMonth, rYear, nStack)
      names(sr) <- c("x", "y", "month_num", "year", names(rasts))
      
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
      srDF$year <- year(my(paste(srDF$month_num, srDF$year, sep = '-')))
      
      #convert dataframe to spatial object 
      srDF = st_as_sf(srDF, coords = c("x", "y"), crs = 4326, agr = "constant")
      srDF = st_sftime(srDF, time_column_name = 'year')
      
      predDF <- pred.rfsi(model = mod,
                          data = stDF, 
                          data.staid.x.y.z = c('staid', 'X', 'Y', 'year'),
                          obs.col = "value",
                          newdata = srDF, 
                          newdata.staid.x.y.z = c('staid', 'X', 'Y', 'year'),
                          output.format = "data.frame", # "sf", # "SpatVector", 
                          cpus = 1, # detectCores()-1,
                          progress = TRUE,
                          classification = F)
      
      hsm[[x]] <- rasterize(x = predDF[,c('X', 'Y')], y = rlon, field = predDF[,'pred'], fun = mean)
      #print(x)
    }
  } #end if rf
  
  if(model == 'brt'){
    print('Predicting Boosted Regression Tree...')
    for(x in 1:nlayers(rasts[[1]])){ #all the rasters in rasts have the same number of layers so it doesn't matter which one we call
      
      #may need to change if names aren't always going to be month.year 
      mm.year <- strsplit(names(rasts[[1]])[x], split = '[.]') #all the rasterbricks in normVars also have the same names so again, doesn't matter which one we call
      
      mm <- match(mm.year[[1]][1],month.abb) 
      rMonth[] <- mm
      yr <- as.numeric(mm.year[[1]][2])
      rYear[] <- yr
      
      nStack <- vector(mode = 'list', length = length(rasts))
      for(n in 1:length(rasts)){
        nStack[[n]] <- subset(rasts[[n]], x)
      }
      nStack <- stack(nStack)
      crs(nStack) <- crs(rasts[[1]])
      extent(nStack) <- extent(rasts[[1]])
      
      sr <- stack(rlon, rlat, rMonth, rYear, nStack)
      names(sr) <- c("x", "y", "month_num", "year", names(rasts))
      
      if(mask){ #if mask == T
        #mask off waters deeper than 1000 m
        i <- which(names(sr) == bathy_nm)
        sr <- replace(sr, abs(subset(sr, i)) > bathy_max, NA) #replace values with an absolute value greater than bathy_max with NA
      }
      
      p <- dismo::predict(object = mod, x = sr,type="response")
      
      hsm[[x]] <- raster::rasterize(x = xy[,1:2], y = rYear, field = p)
      #print(x)
    }
  } #end if brt
  
  if(model == 'tmb'){
    print('Predicting sdmTMB...')
    for(x in 1:nlayers(rasts[[1]])){ #all the rasters in rasts have the same number of layers so it doesn't matter which one we call
      
      #may need to change if names aren't always going to be month.year 
      mm.year <- strsplit(names(rasts[[1]])[x], split = '[.]') #all the rasterbricks in normVars also have the same names so again, doesn't matter which one we call
      
      mm <- match(mm.year[[1]][1],month.abb) 
      rMonth[] <- mm
      yr <- as.numeric(mm.year[[1]][2])
      rYear[] <- yr
      
      nStack <- vector(mode = 'list', length = length(rasts))
      for(n in 1:length(rasts)){
        nStack[[n]] <- subset(rasts[[n]], x)
      }
      nStack <- stack(nStack)
      crs(nStack) <- crs(rasts[[1]])
      extent(nStack) <- extent(rasts[[1]])
      
      sr <- stack(rlon, rlat, rMonth, rYear, nStack)
      names(sr) <- c("x", "y", "month_num", "year", names(rasts))
      
      if(mask){ #if mask == T
        #mask off waters deeper than 1000 m
        i <- which(names(sr) == bathy_nm)
        sr <- replace(sr, abs(subset(sr, i)) > bathy_max, NA) #replace values with an absolute value greater than bathy_max with NA
      }
      
      ##convert rasterStack to dataframe to play well with model 
      srDF <- as.data.frame(rasterToPoints(sr))
      
      sdmpred<- predict(mod, newdata = srDF)
      sdmpred$prob <- exp(sdmpred$est)/(1+exp(sdmpred$est))
      coordinates(sdmpred) <- ~x + y
      proj4string(sdmpred) <- CRS("+proj=longlat +datum=WGS84 +no_defs ")
      hsm[[x]] <- raster::rasterize(x = xy[,1:2], y = rYear, field = 'prob')
      #print(x)
    }
  } #end if sdmtmb
  
  if(model == 'ens'){
    print('Predicting Ensemble...')
    
    if(length(rasts) == length(weights)){
    
      for(x in 1:nlayer(rasts[[1]])){
        #make list of abundances for timestamp
        abunds <- vector(mode = 'list', length = length(rasts))
        for(m in 1:length(rasts)){
          abunds[[m]] <- terra::rast(rasts[[m]][[x]])
        }
        #make weighted average ensembles 
        hsm[[x]] <- raster(MakeEnsembleAbundance(model.weights = weights, abund.list = abunds))
      } #end x
    } else {
      #print error 
    }
  } #end if ensemble

  names(hsm) <- names(rasts[[1]])
  
  return(hsm)
}