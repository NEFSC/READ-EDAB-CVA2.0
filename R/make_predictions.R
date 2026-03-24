#' @title Predict Component or Ensemble SDM
#' @description Make model predictions for either a component model (GAM/MAXENT/RF/BRT/SDMTMB) or the Ensemble model
#'
#' @param mod the output from \code{make_sdm} - only used for ensemble
#' @param model one of the following indicating the desired model to calculate variable importance for: gam, maxent, brt, rf, sdmtmb, or ens
#' @param rasts list of rasterStacks corresponding to the environmental covariates used to build the models. The number of layers in each rasterStack should be the same and correspond to the length of the timeseries for the models to be predicted on
#' @param staticVars list of rasters containing the static variables used in model. Should be the same object used in \code{merge_spp_env}.
#' @param mask TRUE/FALSE indicating whether to mask off certain depths (i.e. waters deeper than 1000 m)
#' @param bathyR,bathy_max a raster of bathymetry with the same extent and resolution as \code{rasts} and the maximum depth you want included. For example, if you want to mask off waters deeper than 1000 m, \code{bathy_max} would be set to 1000. The value should be positive regardless of the sign of your bathymetry data.
#' @param se data frame containing species presence/absence data and desired environmental covariate data.
#' @param xy_col a vector with a length of 2 indicating the longitude and latitude column names
#' @param month_col,year_col column names for month and year columns respectively
#' @param weights a vector of model weights - used for building the ensemble model
#'
#' @return returns a rasterStack of predicted habitat suitability. The number of layers will be equal to the number of layers in \code{rasts}

make_predictions <- function(mod, model, rasts, staticVars, bathyR, mask = T, bathy_max, se = NULL, month_col, year_col, xy_col, weights = NULL){
  if(model %in% c('gam', 'maxent', 'rf', 'brt', 'sdmtmb', 'ens')){


    #make static vars (month/year) into rasters
    r <- raster::subset(rasts[[1]][[1]], 1)
    rlon<-rlat<-r #copy r to rlon and rlat rasters [1]][1]which will contain the longitude and latitude
    xy<-raster::xyFromCell(r,1:length(r)) #matrix of longitudes (x) and latitudes(y)
    rlon[]<-xy[,1] #raster of longitudes
    rlat[]<-xy[,2] #raster of latitudes
    rMonth <- rYear <- r

    hsm <- vector(mode = 'list', length = length(rasts[[1]]))

    #slightly different models depending on the model
    if(model == 'gam'){
      print('Predicting GAM model...')
      for(x in 1:raster::nlayers(rasts[[1]][[1]])){ #all the rasters in rasts have the same number of layers so it doesn't matter which one we call

        #may need to change if names aren't always going to be month.year
        mm.year <- strsplit(names(rasts[[1]][[1]])[x], split = '[.]') #all the rasterbricks in normVars also have the same names so again, doesn't matter which one we call

        mm <- as.numeric(gsub('X', '', mm.year[[1]][1]))
        rMonth[] <- mm
        yr <- as.numeric(mm.year[[1]][2])
        rYear[] <- yr

        nStack <- vector(mode = 'list', length = length(rasts))
        for(n in 1:length(rasts)){
          nStack[[n]] <- raster::subset(rasts[[n]][[1]], x)
        }
        nStack <- c(nStack, staticVars)
        nStack <- raster::stack(nStack)
        names(nStack) <- c(names(rasts), names(staticVars))
        raster::crs(nStack) <- raster::crs(rasts[[1]][[1]])
        raster::extent(nStack) <- raster::extent(rasts[[1]][[1]])

        sr <- raster::stack(rlon, rlat, rMonth, rYear, nStack)
        names(sr)[1:4] <- c("x", "y", "month", "year")

        if(mask){ #if mask == T
          #mask off waters deeper than 1000 m
          #i <- which(names(sr) == bathy_nm)
          sr <- replace(sr, abs(bathyR) > bathy_max, NA) #replace values with an absolute value greater than bathy_max with NA
        }


        hsm[[x]] <- EFHSDM::MakeGAMAbundance(model = mod, r.stack = sr)
        #print(x)
      }
      names(hsm) <- names(rasts[[1]][[1]])
    } #end if gam

    if(model == 'maxent'){
      print('Predicting MAXENT model...')
      for(x in 1:raster::nlayers(rasts[[1]][[1]])){ #all the rasters in rasts have the same number of layers so it doesn't matter which one we call

        #may need to change if names aren't always going to be month.year
        mm.year <- strsplit(names(rasts[[1]][[1]])[x], split = '[.]') #all the rasterbricks in normVars also have the same names so again, doesn't matter which one we call

        mm <- as.numeric(gsub('X', '', mm.year[[1]][1]))
        rMonth[] <- mm
        yr <- as.numeric(mm.year[[1]][2])
        rYear[] <- yr

        nStack <- vector(mode = 'list', length = length(rasts))
        for(n in 1:length(rasts)){
          nStack[[n]] <- raster::subset(rasts[[n]][[1]], x)
        }
        nStack <- c(nStack, staticVars)
        nStack <- raster::stack(nStack)
        names(nStack) <- c(names(rasts), names(staticVars))
        raster::crs(nStack) <- raster::crs(rasts[[1]][[1]])
        raster::extent(nStack) <- raster::extent(rasts[[1]][[1]])

        sr <- raster::stack(rlon, rlat, rMonth, rYear, nStack)
        names(sr)[1:4] <- c("x", "y", "month", "year")

        if(mask){ #if mask == T
          #mask off waters deeper than 1000 m
          #i <- which(names(sr) == bathy_nm)
          sr <- replace(sr, abs(bathyR) > bathy_max, NA) #replace values with an absolute value greater than bathy_max with NA
        }

        hsm[[x]] <- raster::raster(make_maxent_abundance(model = mod, maxent.stack = sr, type = 'maxnet'))
        #print(x)
      }
      names(hsm) <- names(rasts[[1]][[1]])
    } #end if maxent

    if(model == 'rf'){
      print('Predicting RFSI...')
      se <- se[complete.cases(se),]
      se <- cbind(1:nrow(se), se) #stand in station ids
      colnames(se)[1] <- "staid"
      se$month.year <- paste(se[,month_col], se[,year_col], sep = '-')
      se$year <- lubridate::year(lubridate::my(se$month.year))

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
          absSub <- abs[sample(x = nrow(abs), size = round(nrow(abs)/4)),] #subsample absences to a 1/4 of the absences within month and region
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
      stDF = sf::st_as_sf(seSub, coords = xy_col, crs = 4326, agr = "constant")
      stDF = sftime::st_sftime(stDF, time_column_name = month_col)

      for(x in 1:raster::nlayers(rasts[[1]][[1]])){ #all the rasters in rasts have the same number of layers so it doesn't matter which one we call

        #may need to change if names aren't always going to be month.year
        mm.year <- strsplit(names(rasts[[1]][[1]])[x], split = '[.]') #all the rasterbricks in normVars also have the same names so again, doesn't matter which one we call

        mm <- as.numeric(gsub('X', '', mm.year[[1]][1]))
        rMonth[] <- mm
        yr <- as.numeric(mm.year[[1]][2])
        rYear[] <- yr

        nStack <- vector(mode = 'list', length = length(rasts))
        for(n in 1:length(rasts)){
          nStack[[n]] <- raster::subset(rasts[[n]][[1]], x)
        }
        nStack <- c(nStack, staticVars)
        nStack <- raster::stack(nStack)
        names(nStack) <- c(names(rasts), names(staticVars))
        raster::crs(nStack) <- raster::crs(rasts[[1]][[1]])
        raster::extent(nStack) <- raster::extent(rasts[[1]][[1]])

        sr <- raster::stack(rlon, rlat, rMonth, rYear, nStack)
        names(sr)[1:4] <- c("x", "y", "month", "year")

        if(mask){ #if mask == T
          #mask off waters deeper than 1000 m
          #i <- which(names(sr) == bathy_nm)
          sr <- replace(sr, abs(bathyR) > bathy_max, NA) #replace values with an absolute value greater than bathy_max with NA
        }


        ##convert rasterStack to dataframe to play well with model
        srDF <- as.data.frame(raster::rasterToPoints(sr))
        #remove NAs
        # srDF <- srDF[-which(srDF$bathy < -1000 | srDF$bathy > 0),]
        #create staid/year
        srDF$staid <- 1:nrow(srDF) #standin station ids
        #srDF$month <- month(my(paste(srDF$month, srDF$year, sep = '-')))

        #convert dataframe to spatial object
        srDF = sf::st_as_sf(srDF, coords = xy_col, crs = 4326, agr = "constant")
        srDF = sftime::st_sftime(srDF, time_column_name = 'month')

        predDF <- meteo::pred.rfsi(model = mod,
                            data = stDF,
                            data.staid.x.y.z = c('staid', xy_col),
                            obs.col = "value",
                            newdata = srDF,
                            newdata.staid.x.y.z = c('staid', colnames(srDF)[1:2]),
                            output.format = "data.frame", # "sf", # "SpatVector",
                            cpus = 1, # detectCores()-1,
                            progress = TRUE,
                            classification = F)

        raw <- raster::rasterize(x = predDF[,c('X', 'Y')], y = rlon, field = predDF[,'pred'], fun = mean)


        #the RF rasters tend to have some (mostly small) gaps due to shifts in effort, so performing a quick krige interpolation to smooth gaps
        r <- terra::rast(raw)

        # Convert non-NA cells to points for modeling
        pts <- terra::as.points(r, values=TRUE, na.rm=TRUE)
        pts_df <- as.data.frame(pts, geom="XY")
        colnames(pts_df) <- c("val", "x", "y")

        # 2. Model the spatial correlation (The Variogram)
        # We calculate how variance increases with distance
        v_geo <- gstat::variogram(val ~ 1, loc = ~x+y, data = pts_df)

        # Define the anisotropic variogram
        # psill = 0.06 (from your y-axis)
        # range = 3 (from your x-axis on the 90-degree plot)
        # nugget = 0.01 (where the dots start on the y-axis)

        v_fit_anis <- gstat::fit.variogram(v_geo,
                                    model = gstat::vgm(psill = 0.06,
                                                model = "Sph",
                                                range = 3,
                                                nugget = 0.01,
                                                anis = c(90, 0.66)))

        # 3. Perform Kriging
        # Define the grid where you want predictions (your original raster template)
        grid <- terra::rast(r)
        grid_pts <- as.data.frame(grid, xy=TRUE, na.rm=FALSE)
        sp::coordinates(grid_pts) <- ~x+y

        # Convert training data to SpatialPointsDataFrame for gstat compatibility
        sp::coordinates(pts_df) <- ~x+y

        # Run Ordinary Kriging
        kriged_obj <- gstat::krige(val ~ 1, pts_df, grid_pts, model = v_fit_anis, nmax = 50)

        # 4. Convert back to Raster
        # 1. Convert the kriged output to a standard data frame
        kriged_df <- as.data.frame(kriged_obj)

        # 2. Create the SpatRaster directly from the data frame
        # 'x' and 'y' are the coordinates, 'var1.pred' is the Kriging prediction
        r_filled <- terra::rast(kriged_df[, c("x", "y", "var1.pred")], type="xyz")

        # 3. Set the Coordinate Reference System (CRS)
        # This ensures it matches your original data
        terra::crs(r_filled) <- terra::crs(r)

        hsm[[x]] <- raster::raster(r_filled)
      }

    } #end if rf

    if(model == 'brt'){
      print('Predicting Boosted Regression Tree...')
      for(x in 1:raster::nlayers(rasts[[1]][[1]])){ #all the rasters in rasts have the same number of layers so it doesn't matter which one we call

        #may need to change if names aren't always going to be month.year
        mm.year <- strsplit(names(rasts[[1]][[1]])[x], split = '[.]') #all the rasterbricks in normVars also have the same names so again, doesn't matter which one we call

        mm <- as.numeric(gsub('X', '', mm.year[[1]][1]))
        rMonth[] <- mm
        yr <- as.numeric(mm.year[[1]][2])
        rYear[] <- yr

        nStack <- vector(mode = 'list', length = length(rasts))
        for(n in 1:length(rasts)){
          nStack[[n]] <- raster::subset(rasts[[n]][[1]], x)
        }
        nStack <- c(nStack, staticVars)
        nStack <- raster::stack(nStack)
        names(nStack) <- c(names(rasts), names(staticVars))
        raster::crs(nStack) <- raster::crs(rasts[[1]][[1]])
        raster::extent(nStack) <- raster::extent(rasts[[1]][[1]])

        sr <- raster::stack(rlon, rlat, rMonth, rYear, nStack)
        names(sr)[1:4] <- c("x", "y", "month", "year")

        if(mask){ #if mask == T
          #mask off waters deeper than 1000 m
          #i <- which(names(sr) == bathy_nm)
          sr <- replace(sr, abs(bathyR) > bathy_max, NA) #replace values with an absolute value greater than bathy_max with NA
        }

        srDF <- as.data.frame(raster::rasterToPoints(sr))
        srDF <- srDF[complete.cases(srDF),-c(1:2)]
        hsm[[x]] <- predict(object = sr, model = mod, type="response")

        # hsm[[x]] <- raster::rasterize(x = srDF[,1:2], y = rYear, field = p)
        #print(x)
      }
      names(hsm) <- names(rasts[[1]][[1]])
    } #end if brt

    if(model == 'sdmtmb'){
      print('Predicting sdmTMB...')
      warning('sdmTMB predictions can take a long time in some environments. If this call is taking too long (> 1 hr), try switching to predictSDM and makePredDF')
      for(x in 1:raster::nlayers(rasts[[1]][[1]])){ #all the rasters in rasts have the same number of layers so it doesn't matter which one we call

        #may need to change if names aren't always going to be month.year
        mm.year <- strsplit(names(rasts[[1]][[1]])[x], split = '[.]') #all the rasterbricks in normVars also have the same names so again, doesn't matter which one we call

        mm <- as.numeric(gsub('X', '', mm.year[[1]][1]))
        rMonth[] <- mm
        yr <- as.numeric(mm.year[[1]][2])
        rYear[] <- yr

        nStack <- vector(mode = 'list', length = length(rasts))
        for(n in 1:length(rasts)){
          nStack[[n]] <- raster::subset(rasts[[n]][[1]], x)
        }
        nStack <- c(nStack, staticVars)
        nStack <- raster::stack(nStack)
        names(nStack) <- c(names(rasts), names(staticVars))
        raster::crs(nStack) <- raster::crs(rasts[[1]][[1]])
        raster::extent(nStack) <- raster::extent(rasts[[1]][[1]])

        sr <- raster::stack(rlon, rlat, rMonth, rYear, nStack)
        names(sr)[1:4] <- c("x", "y", "month", "year")

        if(mask){ #if mask == T
          #mask off waters deeper than 1000 m
          #i <- which(names(sr) == bathy_nm)
          sr <- replace(sr, abs(bathyR) > bathy_max, NA) #replace values with an absolute value greater than bathy_max with NA
        }

        ##convert rasterStack to dataframe to play well with model
        srDF <- as.data.frame(raster::rasterToPoints(sr)[,-c(1:2)])
        srDF <- srDF[complete.cases(srDF),]

        pred <- predict(mod, newdata = srDF, type = 'response')
        #sdmpred$prob <- exp(sdmpred$est)/(1+exp(sdmpred$est))
        sp::coordinates(pred) <- ~x + y
        sp::proj4string(pred) <- sp::CRS("+proj=longlat +datum=WGS84 +no_defs ")
        hsm[[x]] <- raster::rasterize(x = pred, y = rYear, field = pred$est)
        print(x)
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
          hsm[[x]] <- raster::raster(EFHSDM::MakeEnsembleAbundance(model.weights = wts, abund.list = abunds))
        } #end x
      } else {
        stop('raster list and weights are not the same length')
      }
    } #end if ensemble

    names(hsm) <- names(rasts[[1]][[1]])

    return(raster::stack(hsm))
  } else {
    stop('model is not gam, maxent, rf, brt, sdmtmb, or ens')
  }
}



