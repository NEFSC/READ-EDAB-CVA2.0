#code for paper Camrin et al 2023: https://esajournals.onlinelibrary.wiley.com/doi/10.1002/eap.2893
#source('~/ClimateVulnerabilityAssessment2.0)/Code/eval_kfold_brt_wpredictions.r')
#source('~/ClimateVulnerabilityAssessment2.0)/Code/eval_brt.r')
#source('~/ClimateVulnerabilityAssessment2.0)/Code/pseudoR2.brt.r')
#source('~/ClimateVulnerabilityAssessment2.0)/Code/saveTSS.r')
#source('~/ClimateVulnerabilityAssessment2.0)/Code/saveAUC.r')
#source('~/ClimateVulnerabilityAssessment2.0)/Code/bhattacharyya.stat.r')
#source('~/ClimateVulnerabilityAssessment2.0)/Code/bhatt.coef.r')


make_sdm <- function(se, pa_col, xy_col, month_col, year_col, model, ensembleWeights = NULL, ensemblePreds = NULL){
  #spp-env data frame subset to important covariates with match_guilds 
  #model - one of the following to build models appropriately: gam, maxent, rf, brt, sdmtmb
  #pa_col - name of columns containing response data 
  #xy_col - names of columns containing positional data, where the first is the x(lon) and second is the y(lat)
  #month/year_col - names of columns containing months and years associated with observations 
  #ensembleMods/Mets/Preds - vector or list of evaluation metric to use to generate weights and prediction dataframes from CVs to use to validate model -- ALL NEED TO BE THE SAME LENGTH AND ORDER 
  #metric = name of metric to use for model weights - currently accepts auc or rmse
  
  
  #now build model of choice
  #each produces the following objects: mod, cv, preds, rmse, auc, RDE
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
      
      if(nrow(pres) <= 5){
        allSub <- NULL
      } else if(nrow(abs) > nrow(pres)){
        absSub <- abs[sample(x = nrow(abs), size = nrow(pres)),]
        allSub <- rbind(absSub, pres)
      } else {
        allSub <- sub
      }
      
      seSub <- rbind(seSub, allSub)
    }
    
    #convert dataframe to spatial object 
    stDF = st_as_sf(seSub, coords = xy_col, crs = 4326, agr = "constant")
    stDF = st_sftime(stDF, time_column_name = year_col)
    
    #create formula 
    form <- "value ~ "
    for(x in 2:ncol(se)){
      if(colnames(se)[x] != pa_col & colnames(se)[x] != xy_col[1] & colnames(se)[x] != xy_col[2] & colnames(se)[x] != year_col
         & colnames(se)[x] != 'month.year' & colnames(se)[x] != 'region' & colnames(se)[x] != 'sp.tm' & colnames(se)[x] != 'staid'){ #make sure you don't add the response variable or the variables you've already added
        form <- paste0(form, ' + ', colnames(se)[x])
      }
    } #end for x 
    
    #tune model 
    mod <- rfsi(formula = formula(form),
                data = stDF,
                data.staid.x.y.z = c('staid', 'X', 'Y', 'year'),
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
    
    #make mesh
    mesh <- make_mesh(se, xy_cols = xy_col, cutoff = 6/12) #using lon/lat since this is on the reprojected regular lat/lon grid, and the domain crosses multiple UTM zones  
    #MOM6 resolution is 1/12 = ~8 km 
    
    #build formula
    form <- paste0(pa_col, " ~ ")
    #loop through covariates since they will all have the same smoother 
    for(x in 1:ncol(se)){
      if(colnames(se)[x] != pa_col & colnames(se)[x] != xy_col[1] & colnames(se)[x] != xy_col[2] & colnames(se)[x] != year_col){ #make sure you don't add the response variable or the variables you've already added
        form <- paste0(form, ' + s(', colnames(se)[x], ", bs = 'ts', k = 6)")
      }
    } #end for x 
    
    
    #following general parameters from Andrew Allyn's work with lobster off the coast of Maine 
    #see sp/spatiotemporal model code - https://github.com/aallyn/lobSDM/blob/main/Code/5_AdultLob_SDM.R
    mod <- sdmTMB(
      formula = formula(form),
      data = se,
      mesh = mesh,
      family = binomial(link = 'logit'),
      spatial = "on",
      spatiotemporal = 'ar1',
      time = 'year',
      reml = T,
      anisotropy = F,
      share_range = T,
      do_fit = T
    )
    
    print('Simplifying model...')
    #remove "bad" covariates using edf following methods from EFHSDM 
    sdmEDF <- cAIC(mod, 'EDF')
    i <- which(round(sdmEDF, digits = 3) < 1)
    
    while(length(i) > 0){ #if standard errors are large, try reducing model 
      #it may not fix the problem but it might help
      
      
      cn <- gsub("[\\(\\)]", "", regmatches(names(i), gregexpr("\\(.*?\\)", names(i)))[[1]])
      badvars <- names(se) %in% cn
      print(paste('Removing', names(se)[badvars]))
      se <- se[,-which(names(se) == cn)]
      
      #make new formula 
      #build formula
      form <- paste0(pa_col, " ~ ")
      #loop through covariates since they will all have the same smoother 
      for(x in 1:ncol(se)){
        if(colnames(se)[x] != pa_col & colnames(se)[x] != xy_col[1] & colnames(se)[x] != xy_col[2] & colnames(se)[x] != year_col){ #make sure you don't add the response variable or the variables you've already added
          form <- paste0(form, ' + s(', colnames(se)[x], ", bs = 'ts', k = 6)")
        }
      } #end for x 
      
      #re-run sdm
      mod <- sdmTMB(
        formula = formula(form),
        sea = se,
        mesh = mesh,
        family = binomial(link = 'logit'),
        spatial = "on",
        spatiotemporal = 'ar1',
        time = 'year',
        reml = T,
        anisotropy = F,
        share_range = T,
        do_fit = T
      )
      
      #recheck 
      sdmEDF <- cAIC(mod, 'EDF')
      i <- which(round(sdmEDF, digits = 3) < 1)
    }
    
   
  } #end if sdmtmb
  
  if(model == 'ens'){
   # weights <- MakeEnsemble(rmse = mets) #make weights
    mod <- ValidateEnsemble(pred.list = ensemblePreds, model.weights = weights, make.plots = F) #validate to get preds/obs to metrics 
    
  } #end if ensemble 
  
  return(mod)
  
}

sdm_cv <- function(mod, se, pa_col, xy_col, month_col, year_col, model){
 #mod - model output from make_sdm
   #se - spp-env data frame subset to important covariates with match_guilds 
  #model - one of the following to build models appropriately: gam, maxent, rf, brt, sdmtmb
  #pa_col - name of columns containing response data 
  #xy_col - names of columns containing positional data, where the first is the x(lon) and second is the y(lat)
  #month/year_col - names of columns containing months and years associated with observations 


  if(model == 'gam'){ #build gam model 
    #cross validation
    cv <- EFHSDM::CrossValidateModel(model = mod, data = se, folds = 10, model.type = "gam")
  } #end if gam 
  
  if(model == 'maxent'){
    #cross validate
    cv <- EFHSDM::CrossValidateModel(model = mod, data = se, folds = 10, model.type = "maxnet", species = pa_col, scale.preds = F)
  } #end if maxent
  
  if(model == 'rf'){
    se <- cbind(1:nrow(se), se) #stand in station ids 
    colnames(se)[1] <- "staid"
    se$year <- year(my(paste(se[,month_col], se[,year_col], sep = '-')))
    
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
      
      if(nrow(pres) <= 5){
        allSub <- NULL
      } else if(nrow(abs) > nrow(pres)){
        absSub <- abs[sample(x = nrow(abs), size = nrow(pres)),]
        allSub <- rbind(absSub, pres)
      } else {
        allSub <- sub
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
                  cpus = 3,
                  progress = T,
                  importance = "impurity",
                  seed = 42,
                  acc.metric = 'R2',
                  tgrid = tgrid,
                  tgrid.n= 1,
                  tune.type = "LLO", # Leave-Location-Out CV
                  write.forest = T, 
                  s.crs = st_crs(stDF), 
                  k = 10,
                  classification = F)
    
  } #end if RF
  
  if(model == 'brt'){
    ###simplify model before k-fold validations 
    #simpBRT <- mod
    
    #k-fold validation code from Camrin Brawn (WHOI): https://zenodo.org/records/7971532
    cv <- eval_kfold_brt(dataInput = se, 
                                gbm.x = mod$gbm.call$gbm.x, #only use simplified parameters for k-folds 
                                gbm.y = pa_col,  
                                learning.rate = 0.005, 
                                bag.fraction = 0.75, 
                                tree.complexity = 5, k_folds = 10, max.trees = 2000,
                                is_fixed = F)[[2]]
    
  } #end if BRT
  
  if(model == 'sdmtmb'){
    print('Building sdmTMB...')
    
    #make mesh
    mesh <- make_mesh(se, xy_cols = xy_col, cutoff = 6/12) #using lon/lat since this is on the reprojected regular lat/lon grid, and the domain crosses multiple UTM zones  
    #MOM6 resolution is 1/12 = ~8 km 
    
    #build formula
    form <- paste0(pa_col, " ~ ")
    #loop through covariates since they will all have the same smoother 
    for(x in 1:ncol(se)){
      if(colnames(se)[x] != pa_col & colnames(se)[x] != xy_col[1] & colnames(se)[x] != xy_col[2] & colnames(se)[x] != year_col){ #make sure you don't add the response variable or the variables you've already added
        form <- paste0(form, ' + s(', colnames(se)[x], ", bs = 'ts', k = 6)")
      }
    } #end for x 

    #cross-validation
    cv <- sdmTMB_cv(
      formula = formula(form),
      data = se,
      mesh = mesh,
      family = binomial(link = 'logit'),
      spatial = "on",
      spatiotemporal = 'ar1',
      time = 'year',
      reml = T,
      anisotropy = F,
      share_range = T,
      do_fit = T,
      k_folds = 10,
      parallel = F
    )
    
  } #end if sdmtmb
  
  
  return(cv)
  
}

sdm_preds <- function(cv, model){
  #cv - output from sdm_cv
  #model - one of the following to build models appropriately: gam, maxent, rf, brt, sdmtmb

  
  #now build model of choice
  #each produces the following objects: mod, cv, preds, rmse, auc, RDE
  if(model == 'gam' | model == 'maxent'){ #get preds for gam or maxent
    #get evaluation metrics (RMSE & AUC)
    preds <-cv[[1]]
  } #end if gam/maxent
  
  if(model == 'rf'){
    preds <- cv #the output from cv is also preds here
  } #end if RF
  
  if(model == 'brt'){
    ##put obs and predicted together for ensemble 
    preds <- data.frame(obs = cv$value, preds = cv$preds)
    names(preds) <- c('value', 'preds')
  } #end if BRT
  
  if(model == 'sdmtmb'){
    #change names from sdmTMB preds to make it work with EFHSDM
    preds <- cv$data
    colnames(preds)[9:11] <- c('fold', 'pred', 'loglik')
  } #end if ensemble 
  
  return(preds)
  
}

sdm_eval <- function(preds, model, metric, mod = NULL){
  #preds - results from sdm_preds
  #model - one of the following to build models appropriately: gam, maxent, rf, brt, sdmtmb, ens
  #metric = name of metric to calculate - rmse or auc
  #mod - ensemble model object from make_sdm

  #now evaluate model of choice
  #each produces the following objects: mod, cv, preds, rmse, auc, RDE
  if(model == 'gam' | model == 'maxent'){ #get metrics for gam or maxent model 

    if(metric == 'rmse'){
     #RMSE
      met <- RMSE(obs = preds$abund, pred = preds$cvpred)
    }
    
    if(metric == 'auc'){
    #AUC
      Pred <- prediction(preds$cvpred, preds$abund)
      Perf <- performance(Pred, 'auc')
      met <- Perf@y.values[[1]]
    }
    
  } #end if gam 
  
  if(model == 'rf'){
    
    if(metric == 'rmse'){
    #get RMSE
        met <- EFHSDM::RMSE(obs = preds$obs, pred = preds$pred)
    }
    
    if(metric == 'auc'){
      Pred <- prediction(preds$pred, preds$obs)
      Perf <- performance(Pred, 'auc')
      met <- Perf@y.values[[1]]
    }
    
  } #end if RF
  
  if(model == 'brt'){

    if(metric == 'rmse'){
      met <- EFHSDM::RMSE(obs = preds$value, pred = preds$preds)
    }
    
    if(metric == 'auc'){
      #calculate AUC
      Pred <- prediction(preds$preds, preds$value)
      Perf <- performance(Pred, 'auc')
      met <- Perf@y.values[[1]]
    }

  } #end if BRT
  
  if(model == 'sdmtmb'){

    if(metric == 'rmse'){
    #calculate RMSE
      met <- EFHSDM::RMSE(obs = preds$value, pred = preds$cv_predicted)
    }
    
    if(metric == 'auc'){
      #calculate AUC
      Pred <- prediction(preds$cv_predicted, preds$value)
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
  #se - spp-env data frame subset to important covariates with match_guilds 
  #mod - output from make_sdm
  #model - one of the following to build models appropriately: gam, maxent, rf, brt, sdmtmb
  #pa_col - name of columns containing response data 
  #xy_col - names of columns containing positional data, where the first is the x(lon) and second is the y(lat)
  
  #now build model of choice
  #each produces the following objects: mod, cv, preds, rmse, auc, RDE
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
    se$year <- year(my(paste(se[,month_col], se[,year_col], sep = '-')))
    
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
      
      if(nrow(pres) <= 5){
        allSub <- NULL
      } else if(nrow(abs) > nrow(pres)){
        absSub <- abs[sample(x = nrow(abs), size = nrow(pres)),]
        allSub <- rbind(absSub, pres)
      } else {
        allSub <- sub
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
    
    #make mesh
    mesh <- make_mesh(se, xy_cols = xy_col, cutoff = 6/12) #using lon/lat since this is on the reprojected regular lat/lon grid, and the domain crosses multiple UTM zones  
    #MOM6 resolution is 1/12 = ~8 km 
    
    ### get relative importance of model using type 3 anova method
    # model with *only* intercept and no random fields:
    fit_null <- sdmTMB(
      value ~ 1, 
      spatial = "off", 
      family = binomial(link = 'logit'),
      data = se,
    )
    
    #loop across variables and get their partial deviance explained 
    RDE <- vector(length = ncol(se))
    for(x in 1:ncol(se)){
      if(colnames(se)[x] != pa_col & colnames(se)[x] != xy_col[1] & colnames(se)[x] != xy_col[2] & colnames(se)[x] != year_col){ #don't do this for the response variable, xy vars, or year
        #isolate covariate 
        
        #build formula with single covariate
        formSub <- paste0(pa_col, " ~ + s(", colnames(se)[x], ", bs = 'ts', k = 6)")
        
        fitSub <- sdmTMB(
          formula = formula(formSub),
          data = se,
          mesh = mesh,
          family = binomial(link = 'logit'),
          spatial = "on",
          spatiotemporal = 'ar1',
          time = 'year',
          reml = T,
          anisotropy = F,
          share_range = T,
          do_fit = T
        )
        
        RDE[x] <- 1 - deviance(fitSub) / deviance(fit_null)
        print(x)
      }
    } #end if
    names(RDE) <- names(se)
    
  } #end if sdmtmb

  return(RDE)
  
}
