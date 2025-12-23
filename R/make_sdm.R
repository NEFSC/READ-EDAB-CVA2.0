#' @title SDM Building, Evaluating, and Prediction Functions
#' @description Functions to build, evaluate, and get variable importance from the different model components, including the ensemble model where appropriate. Also includes prediction functions.

#' \itemize{
#' \item \code{make_sdm} is the function that builds one of the five component models  or the ensemble model
#' \item \code{sdm_cv} performs the cross-validation for one of the five component models
#' \item \code{sdm_preds} isolates the predictions from the cross-validation
#' \item \code{sdm_eval} calculates a performance metric (RMSE or AUC) from the predicted values generated in the cross-validation
#' \item \code{sdm_importance} calculates variable importance for one of the five component models 
#' \item \code{make_predictions} predicts the provided model over the provided timeseries of environmental data
#' \item \code{makePredDF} and \code{predictSDM} are alternatives to \code{make_predictions} specifically for the sdmTMB models. Instead of predicting each timestep in \code{rasts} like \code{make_predictions} does, \code{makePredDF} makes one large data.frame containing all the data in \code{rasts}. A single prediction call can then be applied to the entire timeseries, then \code{predictSDM} will generate the raster layers. 
#' }

#' @param se data frame containing species presence/absence data and desired environmental covariate data.
#' @param pa_col column name for presence/absence column
#' @param xy_col a vector with a length of 2 indicating the longitude and latitude column names
#' @param month_col,year_col column names for month and year columns respectively
#' @param model one of the following indicating the desired model to build: gam, maxent, brt, rf, sdmtmb, or ens
#' @param ensembleWeights,weights vector of model weights used to build ensemble model. The vector must have the same length and be in the same order as the corresponding predictions list.
#' @param ensemblePreds a list of prediction values from component models used to build ensemble model. The list of predictions must have the same length and be in the same order as the corresponding weight vector. 
#' @param mod the output from \code{make_sdm}
#' @param cv output from \code{sdm_cv}
#' @param preds output from \code{sdm_preds}
#' @param rasts list of rasterStacks corresponding to the environmental covariates used to build the models. The number of layers in each rasterStack should be the same and correspond to the length of the timeseries for the models to be predicted on
#' @param staticData local file path to static variables used in model. Should be the same object used in \code{merge_spp_env}. 
#' @param mask TRUE/FALSE indicating whether to mask off certain depths (i.e. waters deeper than 1000 m)
#' @param bathyR,bathy_max a raster of bathymetry with the same extent and resolution as \code{rasts} and the maximum depth you want included. For example, if you want to mask off waters deeper than 1000 m, \code{bathy_max} would be set to 1000. The value should be positive regardless of the sign of your bathymetry data.
#' @param df the output from \code{makePredDF}

#' @return \code{make_sdm} returns the model object from the desired model. Model objects will differ depending on the model type. 
#' @return \code{sdm_cv} returns the cross-validation results from the desired model 
#' @return \code{sdm_preds} returns the prediction outputs from the cross-validation necesary to calculate evaluation metric
#' @return \code{sdm_eval} returns the desired evaluation metric for the given model
#' @return \code{sdm_importance} returns the variable importance for the given model. Each model calculates these differently, so the values should be normalized in order to compare across models. 
#' @return \code{make_predictions} and \code{predictSDM} returns a rasterStack with the same number of layers as in \code{rasts} with probability of occurance predicted for each day with available data. Values will range from 0 to 1. 
#' @return \code{makePredDF} returns a data.frame containing all of the environmental data in \code{rasts} for all timesteps. It will be big depending on the length of your time series. This is meant to be used to help predict sdmTMB models. 

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
    mod <- FitGAM(gam.formula = formula(form), data = se, family.gam = "binomial", select = T, reduce = T)
    
  } #end if gam 
  
  if(model == 'maxent'){
    print('Building MAXENT...')
    #run model - no formula needed
    mod <- FitMaxnet(data = se, species = pa_col, vars = names(se)[-which(names(se) == pa_col)], reduce = T) #fit maxent with all covariates in dataframe since we've already subset the dataframe to be only relevant covariates 
    
  } #end if maxent
  
  if(model == 'rf'){
    print('Building Random Forest with Spatial Interpolation...')
    
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
    stDF = st_as_sf(seSub, coords = xy_col, crs = 4326, agr = "constant")
    stDF = st_sftime(stDF, time_column_name =  month_col)
    
    #create formula 
    form <- "value ~ "
    for(x in 2:ncol(se)){
      if(colnames(se)[x] != pa_col & colnames(se)[x] != xy_col[1] & colnames(se)[x] != xy_col[2] & colnames(se)[x] != year_col & colnames(se)[x] != month_col 
         & colnames(se)[x] != 'month.year' & colnames(se)[x] != 'region' & colnames(se)[x] != 'sp.tm' & colnames(se)[x] != 'staid'){ #make sure you don't add the response variable or the variables you've already added
        form <- paste0(form, ' + ', colnames(se)[x])
      }
    } #end for x 
    
    #tune model 
    mod <- rfsi(formula = formula(form),
                data = stDF,
                data.staid.x.y.z = c('staid', 'X', 'Y'),
                cpus = 1,
                progress = F,
                importance = "impurity",
                seed = 42,
                num.trees = 200,
                s.crs = st_crs(stDF), 
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
    
    modAll <- gbm.step(data = se, gbm.x = names(se)[-which(names(se) == pa_col)], 
                    gbm.y = pa_col, 
                    family = 'bernoulli', tree.complexity = 5,
                    learning.rate = 0.005, bag.fraction = 0.75, max.trees = 2000, n.folds = 10) #using the same parameters as Braun et al 2023
    
    ###simplify model before k-fold validations 
    simpBRT <- gbm.simplify(modAll)
    
    #re-run model with best parameters
    mod <- gbm.step(data = se, gbm.x = simpBRT$pred.list[[length(simpBRT$pred.list)]], 
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
    mesh <- make_mesh(se, xy_cols = xy_col, cutoff = 1) #using lon/lat since this is on the reprojected regular lat/lon grid, and the domain crosses multiple UTM zones  
    #MOM6 resolution is 1/12 = ~8 km 

      ##add extra years to help with forecasting
      forecast_years <- 2020:2035 #probably more than we need but to be safe - will need to add options to define

      
    #following general parameters from Andrew Allyn's work with lobster off the coast of Maine 
    #see sp/spatiotemporal model code - https://github.com/aallyn/lobSDM/blob/main/Code/5_AdultLob_SDM.R
   tryMod <- tryCatch(
       expr = {
      mod <- sdmTMB(
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
                sdmEDF <- cAIC(mod, 'EDF')
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
                  mesh <- make_mesh(se, xy_cols = xy_col, cutoff = 1)
                  
                  #re-run sdm
                  modS <- sdmTMB(
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
                  sdmEDF <- cAIC(modS, 'EDF')
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
    mod <- ValidateEnsemble(pred.list = ensemblePreds, model.weights = weights, make.plots = F, latlon = F) #validate to get preds/obs to metrics 
    
  } #end if ensemble 
  
  return(mod)
  
}

sdm_cv <- function(mod, se, pa_col, xy_col, month_col, year_col, model){

  if(model == 'gam'){ #build gam model 
    #cross validation
    cv <- EFHSDM::CrossValidateModel(model = mod, data = se, folds = 5, model.type = "gam")
  } #end if gam 
  
  if(model == 'maxent'){
    #cross validate
    cv <- EFHSDM::CrossValidateModel(model = mod, data = se, folds = 5, model.type = "maxnet", species = pa_col, scale.preds = F)
  } #end if maxent
  
  if(model == 'rf'){
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
    stDF = st_as_sf(seSub[complete.cases(seSub),], coords = xy_col, crs = 4326, agr = "constant")
    stDF = st_sftime(stDF, time_column_name = month_col)
    
    
    #create formula 
    form <- "value ~ "
    for(x in 2:ncol(se)){
      if(colnames(se)[x] != pa_col & colnames(se)[x] != xy_col[1] & colnames(se)[x] != xy_col[2] & colnames(se)[x] != year_col & colnames(se)[x] != month_col & 
         colnames(se)[x] != 'month.year' & colnames(se)[x] != 'region' & colnames(se)[x] != 'sp.tm' & colnames(se)[x] != 'staid'){ #make sure you don't add the response variable or the variables you've already added
        form <- paste0(form, ' + ', colnames(se)[x])
      }
    } #end for x 

    
    #tuning parameters for model
    n.obs <- 5
    min.node.size <- mod$min.node.size
    sample.fraction <- 0.95 
    splitrule <- "variance"
    ntree <- mod$num.trees # 500
    mtry <- mod$mtry
    tgrid = data.frame(min.node.size=min.node.size, num.trees=ntree,
                       mtry=mtry, n.obs=n.obs, sample.fraction=sample.fraction)
    
    # cross-validate (this will also run tune.rfsi, but it doesn't give you to the output like the raw code suggests! #rude)
    cv <- cv.rfsi(formula = formula(form),
                  data = stDF,
                  data.staid.x.y.z = c('staid', 'X', 'Y', 'year'),
                  cpus = 1,
                  progress = T,
                  importance = "impurity",
                  seed = 42,
                  acc.metric = 'R2',
                  tgrid = tgrid,
                  tgrid.n= 1,
                  tune.type = "LLO", # Leave-Location-Out CV
                  write.forest = T, 
                  s.crs = st_crs(stDF), 
                  k = 5,
                  classification = F)
    
  } #end if RF
  
  if(model == 'brt'){
    ###simplify model before k-fold validations 
    #simpBRT <- mod
    
    se <- se[complete.cases(se),]
    
    #k-fold validation code from Camrin Brawn (WHOI): https://zenodo.org/records/7971532
    cv <- eval_kfold_brt(dataInput = se, 
                                gbm.x = mod$gbm.call$gbm.x, #only use simplified parameters for k-folds 
                                gbm.y = pa_col,  
                                learning.rate = 0.005, 
                                bag.fraction = 0.75, 
                                tree.complexity = 5, k_folds = 5, max.trees = 2000,
                                is_fixed = F)[[2]]
    
  } #end if BRT
  
  if(model == 'sdmtmb'){
    print('Building sdmTMB...')
    #build formula
    #form <- paste0(pa_col, " ~ ")
    #loop through covariates since they will all have the same smoother 
    #for(x in 1:ncol(se)){
     # if(colnames(se)[x] != pa_col & colnames(se)[x] != xy_col[1] & colnames(se)[x] != xy_col[2] & colnames(se)[x] != year_col){ #make sure you don't add the response variable or the variables you've already added
      #  form <- paste0(form, ' + s(', colnames(se)[x], ", k = 6)")
      #}
    #} #end for x 
    
    se2 <- se[complete.cases(se),c(all.vars(formula(mod$formula[[1]])), year_col, xy_col)]
    
    #cross-validation - doing it manually for memory reasons - following format from Braun et al scripts
      se2$Kset <- dismo::kfold(se2, 5) #randomly allocate k groups
    kRes <- NULL 
    for (i in 1:5){
        print(i)
        #seperate training and test sets 
        train <- se2[se2$Kset!=i,]
        test <- se2[se2$Kset==i,]

            #make mesh
    mesh <- make_mesh(train, xy_cols = xy_col, cutoff = 1) #using lon/lat since this is on the reprojected regular lat/lon grid, and the domain crosses multiple UTM zones  
    #MOM6 resolution is 1/12 = ~8 km 
    tryMod <- tryCatch(expr = {
        #create model 
        cv <- sdmTMB(
          formula = formula(mod$formula[[1]]),
          data = train,
          mesh = mesh,
          family = binomial(link = 'logit'),
          #spatial = "on",
          spatiotemporal = 'iid',
          time = 'year',
          reml = T,
          anisotropy = T,
          share_range = F,
          do_fit = T
        )

        #predict model
         predicted <- predict(cv, newdata = test, type = "response")
        rm(cv) #to help with memory
        #return(predicted)
            },  #end expr
        error = function(e){ #to do if fails - this essentially will skip the 
            print('could not build model for test - be sure to test how many iterations failed. may need to be rerun if 2+ fail.')
            test$est <- test$est_non_rf <- test$est_rf <- test$omega_s <- test$epsilon_st <- NA #create NA vectors that would have been added if model was predicted 
            predicted <- test
               # return(predicted)
        },#end error function
            finally = {kRes <- rbind(kRes, predicted)}
                      #return(kRes)} #always put these together
                           ) #end tryCatch
        ###add training values and predictions to calculate RMSE - written by KLG & added 6/20/2025
	
        } #end for i
    cv <- kRes
  } #end if sdmtmb
  
  
  return(cv)
  
}

sdm_preds <- function(cv, model){

  if(model == 'gam' | model == 'maxent'){ #get preds for gam or maxent
    #get evaluation metrics (RMSE & AUC)
    preds <-cv[[1]]
  } #end if gam/maxent
  
  if(model == 'rf'){
    preds <- cv #the output from cv is also preds here
    colnames(preds)[6:7] <- c('abund', 'pred')
  } #end if RF
  
  if(model == 'brt'){
    ##put obs and predicted together for ensemble 
    preds <- data.frame(abund = cv$value, pred = cv$preds)
    #names(preds) <- c('abund', 'pred')
  } #end if BRT
  
  if(model == 'sdmtmb'){
    #change names from sdmTMB preds to make it work with EFHSDM
    #preds <- cv$data
      preds <- cv
    colnames(preds)[17] <- c('pred')
  } #end if ensemble 
  
  return(preds)
  
}

sdm_eval <- function(preds, model, metric, mod = NULL){

  if(model == 'gam' | model == 'maxent'){ #get metrics for gam or maxent model 
    preds2 <- preds[complete.cases(preds),]
    if(metric == 'rmse'){
     #RMSE
      met <- RMSE(obs = preds2$abund, pred = preds2$cvpred)
    }
    
    if(metric == 'auc'){
    #AUC
      Pred <- prediction(preds2$cvpred, preds2$abund)
      Perf <- performance(Pred, 'auc')
      met <- Perf@y.values[[1]]
    }
    
  } #end if gam 
  
  if(model == 'rf'){
    
    if(metric == 'rmse'){
    #get RMSE
        met <- EFHSDM::RMSE(obs = preds$abund, pred = preds$pred)
    }
    
    if(metric == 'auc'){
      Pred <- prediction(preds$pred, preds$abund)
      Perf <- performance(Pred, 'auc')
      met <- Perf@y.values[[1]]
    }
    
  } #end if RF
  
  if(model == 'brt'){

    if(metric == 'rmse'){
      met <- EFHSDM::RMSE(obs = preds$abund, pred = preds$pred)
    }
    
    if(metric == 'auc'){
      #calculate AUC
      Pred <- prediction(preds$pred, preds$abund)
      Perf <- performance(Pred, 'auc')
      met <- Perf@y.values[[1]]
    }

  } #end if BRT
  
  if(model == 'sdmtmb'){

    if(metric == 'rmse'){
    #calculate RMSE
      met <- EFHSDM::RMSE(obs = preds$value, pred = preds$pred)
    }
    
    if(metric == 'auc'){
      #calculate AUC
      Pred <- prediction(preds$pred, preds$value)
      Perf <- performance(Pred, 'auc')
      met <- Perf@y.values[[1]]
    }
  } #end if sdmtmb
  
  if(model == 'ens'){
    if(metric == 'rmse'){
    #RMSE
      met <- RMSE(obs = mod$abund, pred = mod$pred)
    }
    
    if(metric == 'auc'){
      Pred <- prediction(mod$pred, mod$abund)
      Perf <- performance(Pred, 'auc')
      met <- Perf@y.values[[1]]
    }
    
  } #end if ensemble 
  
  return(met)
  
}

sdm_importance <- function(mod, se, pa_col, xy_col, month_col, year_col, model){

  if(model == 'gam'){ #build gam model 
    #extract relative deviance explained
    RDE <- GAMStats(model = mod, data = se)
    
  } #end if gam 
  
  if(model == 'maxent'){
    ##get relative deviance explained
    RDE <- MaxnetStats(model = mod, data = se, species = pa_col)
  } #end if maxent
  
  if(model == 'rf'){
    print('Building Random Forest with Spatial Interpolation...')
    
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
    
    
    #create formula 
    form <- "value ~ "
    for(x in 2:ncol(se)){
      if(colnames(se)[x] != pa_col & colnames(se)[x] != xy_col[1] & colnames(se)[x] != xy_col[2] & colnames(se)[x] != year_col &
         colnames(se)[x] != 'month.year' & colnames(se)[x] != 'region' & colnames(se)[x] != 'sp.tm' & colnames(se)[x] != 'staid'){ #make sure you don't add the response variable or the variables you've already added
        form <- paste0(form, ' + ', colnames(se)[x])
      }
    } #end for x 

    ##get important covariates 
    RDE <- importance(x=mod, method = 'altmann', formula = formula(form), data = stDF)
    
  } #end if RF
  
  if(model == 'brt'){
    #get relative importance
    RDE <- mod$contributions
  } #end if BRT
  
  if(model == 'sdmtmb'){
    se <- se[complete.cases(se),]
      
    #make mesh
    mesh <- make_mesh(se, xy_cols = xy_col, cutoff = 1) #using lon/lat since this is on the reprojected regular lat/lon grid, and the domain crosses multiple UTM zones  
    #MOM6 resolution is 1/12 = ~8 km 
    
    ### get relative importance of model using type 3 anova method
    # model with *only* intercept and no random fields:
    fit_null <- sdmTMB(
      formula(paste0(pa_col, " ~ 1")), 
      spatial = "off", 
      family = binomial(link = 'logit'),
      data = se, mesh = mesh
    )
    
    #loop across variables and get their partial deviance explained 
    RDE <- vector(length = ncol(se))
    for(y in 1:ncol(se)){
      if(colnames(se)[y] != pa_col & colnames(se)[y] != xy_col[1] & colnames(se)[y] != xy_col[2] & colnames(se)[y] != year_col & colnames(se)[y] != 'region' & colnames(se)[y] != 'sp.tm'){ #don't do this for the response variable, xy vars, or year
        #isolate covariate 
        
        #build formula with single covariate
        formSub <- paste0(pa_col, " ~ + s(", colnames(se)[y], ", k = 6)")

          # model with *only* variable of choice and no random fields (can get random fields later):
        fitSub <- sdmTMB(
          formula = formula(formSub),
          spatial = "off", 
          family = binomial(link = 'logit'),
          data = se, 
            mesh = mesh
        )
        
        RDE[y] <- 1 - deviance(fitSub) / deviance(fit_null)
        print(y)
      } #end if
    } #end for
    names(RDE) <- names(se)
    
  } #end if sdmtmb

  return(RDE)
  
}


make_predictions <- function(mod, model, rasts, staticData, bathyR, mask = T, bathy_max, se = NULL, month_col, year_col, xy_col, weights = NULL){
  if(model %in% c('gam', 'maxent', 'rf', 'brt', 'sdmtmb', 'ens')){
  

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
      #i <- which(names(sr) == bathy_nm)
      sr <- replace(sr, abs(bathyR) > bathy_max, NA) #replace values with an absolute value greater than bathy_max with NA
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
        #i <- which(names(sr) == bathy_nm)
        sr <- replace(sr, abs(bathyR) > bathy_max, NA) #replace values with an absolute value greater than bathy_max with NA
      }
      
     # require(gbm3)
      hsm[[x]] <- raster(MakeMaxEntAbundance2(model = mod, maxent.stack = sr, type = 'maxnet'))
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
    stDF = st_as_sf(seSub, coords = xy_col, crs = 4326, agr = "constant")
    stDF = st_sftime(stDF, time_column_name = month_col)
    
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
        #i <- which(names(sr) == bathy_nm)
        sr <- replace(sr, abs(bathyR) > bathy_max, NA) #replace values with an absolute value greater than bathy_max with NA
      }
      
      
      ##convert rasterStack to dataframe to play well with model 
      srDF <- as.data.frame(rasterToPoints(sr))
      #remove NAs
      # srDF <- srDF[-which(srDF$bathy < -1000 | srDF$bathy > 0),]
      #create staid/year 
      srDF$staid <- 1:nrow(srDF) #standin station ids 
      #srDF$month <- month(my(paste(srDF$month, srDF$year, sep = '-')))
      
      #convert dataframe to spatial object 
      srDF = st_as_sf(srDF, coords = xy_col, crs = 4326, agr = "constant")
      srDF = st_sftime(srDF, time_column_name = 'month')
      
      predDF <- pred.rfsi(model = mod,
                          data = stDF, 
                          data.staid.x.y.z = c('staid', xy_col),
                          obs.col = "value",
                          newdata = srDF, 
                          newdata.staid.x.y.z = c('staid', colnames(srDF)[1:2]),
                          output.format = "data.frame", # "sf", # "SpatVector", 
                          cpus = 1, # detectCores()-1,
                          progress = TRUE,
                          classification = F)
      
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
        #i <- which(names(sr) == bathy_nm)
        sr <- replace(sr, abs(bathyR) > bathy_max, NA) #replace values with an absolute value greater than bathy_max with NA
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
    warning('sdmTMB predictions can take a long time in some environments. If this call is taking too long (> 1 hr), try switching to predictSDM and makePredDF')
    require(sdmTMB)
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
        #i <- which(names(sr) == bathy_nm)
        sr <- replace(sr, abs(bathyR) > bathy_max, NA) #replace values with an absolute value greater than bathy_max with NA
      }
      
      ##convert rasterStack to dataframe to play well with model 
      srDF <- as.data.frame(rasterToPoints(sr)[,-c(1:2)])
      srDF <- srDF[complete.cases(srDF),]

      require(sdmTMB)
      pred <- predict(mod, newdata = srDF, type = 'response')
      #sdmpred$prob <- exp(sdmpred$est)/(1+exp(sdmpred$est))
      coordinates(pred) <- ~x + y
      proj4string(pred) <- CRS("+proj=longlat +datum=WGS84 +no_defs ")
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
        hsm[[x]] <- raster(MakeEnsembleAbundance(model.weights = wts, abund.list = abunds))
      } #end x
    } else {
      stop('raster list and weights are not the same length')
    }
  } #end if ensemble

  names(hsm) <- names(rasts[[1]][[1]])
  
  return(hsm)
  } else {
    stop('model is not gam, maxent, rf, brt, sdmtmb, or ens')
  }
}



makePredDF <- function(rasts, staticData, bathyR, bathy_max, mask){

      #make static vars (month/year) into rasters
  r <- subset(rasts[[1]][[1]], 1)
  rlon<-rlat<-r #copy r to rlon and rlat rasters [1]][1]which will contain the longitude and latitude
  xy<-xyFromCell(r,1:length(r)) #matrix of longitudes (x) and latitudes(y)
  rlon[]<-xy[,1] #raster of longitudes
  rlat[]<-xy[,2] #raster of latitides
  rMonth <- rYear <- r
  
    load(staticData) #staticVars object containing a list of rasters with the same extent as the environmental variables
    allDF <- NULL
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
            #i <- which(names(sr) == bathy_nm)
            sr <- replace(sr, abs(bathyR) > bathy_max, NA) #replace values with an absolute value greater than bathy_max with NA
          }
          
          ##convert rasterStack to dataframe to play well with model 
          srDF <- as.data.frame(rasterToPoints(sr)[,-c(1:2)])
          srDF <- srDF[complete.cases(srDF),]
            allDF <- rbind(allDF, srDF)
        }
          save(allDF, file = paste0('./Data/prediction_dataframe_', min(allDF$year, na.rm = T), '_', max(allDF$year, na.rm = T), '.RData'))
          print('allDF created and saved in Data directory')
            return(allDF)
    } #end function


predictSDM <- function(mod, df, staticData){
        load(staticData) #staticVars object containing a list of rasters with the same extent as the environmental variables
          Sys.time()
            require(sdmTMB)
          
          hsm <- vector(mode = 'list', length = length(unique(df$my)))
          for(x in 1:length(unique(df$my))){
            sub <- df[df$my == unique(df$my)[x],]
            coordinates(sub) <- ~x + y
            proj4string(sub) <- CRS("+proj=longlat +datum=WGS84 +no_defs ")
            hsm[[x]] <- raster::rasterize(x = sub, y = staticVars[[1]], field = sub$est)
            print(x)
          }
          
          return(hsm)
    }

