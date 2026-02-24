#' @title Build Component or Ensemble Species Distribution Models
#' @description Build one of the five component models or the final ensemble model
#'
#' @param se data frame containing species presence/absence data and desired environmental covariate data.
#' @param pa_col column name for presence/absence column
#' @param xy_col a vector with a length of 2 indicating the longitude and latitude column names
#' @param month_col,year_col column names for month and year columns respectively
#' @param model one of the following indicating the desired model to build: gam, maxent, brt, rf, sdmtmb, or ens
#' @param ensembleWeights vector of model weights used to build ensemble model. The vector must have the same length and be in the same order as the corresponding predictions list.
#' @param ensemblePreds a list of prediction values from component models used to build ensemble model. The list of predictions must have the same length and be in the same order as the corresponding weight vector.
#'
#' @return the model object from the desired model. Model objects will differ depending on the model type.

make_sdm <- function(se, pa_col, xy_col, month_col, year_col, model, ensembleWeights = NULL, ensemblePreds = NULL){
  if(mod != 'gam' | mod != 'maxent' | mod != 'rf' | mod != 'brt' | mod != 'sdmtmb' | mod != 'ens'){
    stop('model is not gam, maxent, rf, brt, sdmtmb, or ens')
  }

  #now build model of choice
  if(model == 'gam'){ #build gam model
    print('Building GAM...')
    #build formula
    form <- paste0(pa_col, " ~ ")
    if(!is.null(xy_col)){
      form <- paste0(form, 's(', xy_col[1], ',', xy_col[2], ", bs = 'ts', k = 10)")
    }
    if(!is.null(month_col)){
      form <- paste0(form, "+ s(", month_col, ", bs = 'cc', k = 6)")
    }
    #loop through remaining covariates since they will all have the same smoother
    for(x in 1:ncol(se)){
      if(colnames(se)[x] != pa_col & colnames(se)[x] != month_col & colnames(se)[x] != xy_col[1] & colnames(se)[x] != xy_col[2]){ #make sure you don't add the response variable or the variables you've already added
        form <- paste0(form, ' + s(', colnames(se)[x], ", bs = 'ts', k = 6)")
      }
    } #end for x

    #run model
    mod <- EFHSDM::FitGAM(gam.formula = formula(form), data = se, family.gam = "binomial", select = T, reduce = T)

  } #end if gam

  if(model == 'maxent'){
    print('Building MAXENT...')
    #run model - no formula needed
    mod <- EFHSDM::FitMaxnet(data = se, species = pa_col, vars = names(se)[-which(names(se) == pa_col)], reduce = T) #fit maxent with all covariates in dataframe since we've already subset the dataframe to be only relevant covariates

  } #end if maxent

  if(model == 'rf'){
    print('Building Random Forest with Spatial Interpolation...')

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
    stDF = sftime::st_sftime(stDF, time_column_name =  month_col)

    #create formula
    form <- "value ~ "
    for(x in 2:ncol(se)){
      if(colnames(se)[x] != pa_col & colnames(se)[x] != xy_col[1] & colnames(se)[x] != xy_col[2] & colnames(se)[x] != year_col & colnames(se)[x] != month_col
         & colnames(se)[x] != 'month.year' & colnames(se)[x] != 'region' & colnames(se)[x] != 'sp.tm' & colnames(se)[x] != 'staid'){ #make sure you don't add the response variable or the variables you've already added
        form <- paste0(form, ' + ', colnames(se)[x])
      }
    } #end for x

    #tune model
    mod <- meteo::rfsi(formula = formula(form),
                data = stDF,
                data.staid.x.y.z = c('staid', 'X', 'Y'),
                cpus = 1,
                progress = F,
                importance = "impurity",
                seed = 42,
                num.trees = 200,
                s.crs = sf::st_crs(stDF),
                use.idw = T,
                classification = F,
                write.forest = T,
                splitrule = "variance",
                min.node.size = 5,
                sample.fraction = 0.95)

  } #end if RF

  if(model == 'brt'){
    print('Building Boosted Regression Trees...')

    se <- se[complete.cases(se),]

    modAll <- dismo::gbm.step(data = se, gbm.x = names(se)[-which(names(se) == pa_col)],
                       gbm.y = pa_col,
                       family = 'bernoulli', tree.complexity = 5,
                       learning.rate = 0.005, bag.fraction = 0.75, max.trees = 2000, n.folds = 10) #using the same parameters as Braun et al 2023

    ###simplify model before k-fold validations
    simpBRT <- dismo::gbm.simplify(modAll)

    #re-run model with best parameters
    mod <- dismo::gbm.step(data = se, gbm.x = simpBRT$pred.list[[length(simpBRT$pred.list)]],
                    gbm.y = pa_col,
                    family = 'bernoulli', tree.complexity = 5,
                    learning.rate = 0.005, bag.fraction = 0.75, max.trees = 2000, n.folds = 10)


  } #end if BRT

  if(model == 'sdmtmb'){
    print('Building sdmTMB...')

    #build formula
    form <- paste0(pa_col, " ~ ")
    #loop through covariates since they will all have the same smoother
    for(x in 1:ncol(se)){
      if(colnames(se)[x] != pa_col & colnames(se)[x] != xy_col[1] & colnames(se)[x] != xy_col[2] & colnames(se)[x] != year_col){ #make sure you don't add the response variable or the variables you've already added
        form <- paste0(form, ' + s(', colnames(se)[x], ", k = 6)")
      }
    } #end for x


    se <- se[complete.cases(se),]

    #make mesh
    mesh <- sdmTMB::make_mesh(se, xy_cols = xy_col, cutoff = 1) #using lon/lat since this is on the reprojected regular lat/lon grid, and the domain crosses multiple UTM zones
    #MOM6 resolution is 1/12 = ~8 km

    ##add extra years to help with forecasting
    forecast_years <- 2020:2035 #probably more than we need but to be safe - will need to add options to define


    #following general parameters from Andrew Allyn's work with lobster off the coast of Maine
    #see sp/spatiotemporal model code - https://github.com/aallyn/lobSDM/blob/main/Code/5_AdultLob_SDM.R
    tryMod <- tryCatch(
      expr = {
        mod <- sdmTMB::sdmTMB(
          formula = formula(form),
          data = se,
          mesh = mesh,
          family = binomial(link = 'logit'),
          #spatial = "on",
          spatiotemporal = 'iid',
          time = 'year',
          reml = T,
          anisotropy = T,
          share_range = F,
          do_fit = T,
          extra_time = forecast_years
        )
      },
      #end expr
      error = function(e){
        message('model did not converge')
        return(NA)
      }
    ) #end tryCatch

    if(exists('mod')){ #if modS exists == model was simplified; it will not exist if simplifying fails

      try2simp <- tryCatch(
        expr = {
          print('Simplifying model...')
          #remove "bad" covariates using edf following methods from EFHSDM
          sdmEDF <- sdmTMB::cAIC(mod, 'EDF')
          i <- which(round(sdmEDF, digits = 3) < 1)

          while(length(i) > 0){ #if standard errors are large, try reducing model
            #it may not fix the problem but it might help


            cn <- gsub('s\\(', '', names(i))
            badvars <- names(se) %in% cn
            print(paste('Removing', names(se)[badvars]))
            se <- se[,-which(badvars == TRUE)]

            #make new formula
            #build formula
            form <- paste0(pa_col, " ~ ")
            #loop through covariates since they will all have the same smoother
            for(x in 1:ncol(se)){
              if(colnames(se)[x] != pa_col & colnames(se)[x] != xy_col[1] & colnames(se)[x] != xy_col[2] & colnames(se)[x] != year_col & colnames(se)[x] != 'sp.tm' & colnames(se)[x] != 'region'){ #make sure you don't add the responseSub variable or the variables you've already added
                form <- paste0(form, ' + s(', colnames(se)[x], ", k = 6)")
              }
            } #end for x

            #make mesh
            mesh <- sdmTMB::make_mesh(se, xy_cols = xy_col, cutoff = 1)

            #re-run sdm
            modS <- sdmTMB::sdmTMB(
              formula = formula(form),
              data = se,
              mesh = mesh,
              family = binomial(link = 'logit'),
              #spatial = "on",
              spatiotemporal = 'iid',
              time = 'year',
              reml = T,
              anisotropy = T,
              share_range = F,
              do_fit = T,
              extra_time = forecast_years
            )

            #recheck
            sdmEDF <- sdmTMB::cAIC(modS, 'EDF')
            i <- which(round(sdmEDF, digits = 3) < 1)
          }
        }, #end expr
        error = function(e){
          message('model could not be simplified')
          return(NA)
        }
      ) #end tryCatch
    }  else { #end if mod exists
      mod <- NA
    }

    if(exists('modS')){ #if modS exists == model was simplified; it will not exist if simplifying fails
      mod <- modS   #overwrite existing model
    }


  } #end if sdmtmb

  if(model == 'ens'){
    # weights <- MakeEnsemble(rmse = mets) #make weights
    mod <- EFHSDM::ValidateEnsemble(pred.list = ensemblePreds, model.weights = ensembleWeights, make.plots = F, latlon = F) #validate to get preds/obs to metrics

  } #end if ensemble

  return(mod)

}

