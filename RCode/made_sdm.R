#code for paper Camrin et al 2023: https://esajournals.onlinelibrary.wiley.com/doi/10.1002/eap.2893
source('~/ClimateVulnerabilityAssessment2.0)/Code/eval_kfold_brt_wpredictions.r')
source('~/ClimateVulnerabilityAssessment2.0)/Code/eval_brt.r')
source('~/ClimateVulnerabilityAssessment2.0)/Code/pseudoR2.brt.r')
source('~/ClimateVulnerabilityAssessment2.0)/Code/saveTSS.r')
source('~/ClimateVulnerabilityAssessment2.0)/Code/saveAUC.r')
source('~/ClimateVulnerabilityAssessment2.0)/Code/bhattacharyya.stat.r')
source('~/ClimateVulnerabilityAssessment2.0)/Code/bhatt.coef.r')


make_sdm <- function(se, pa_col, xy_col, month_col, year_col, model, metricVec, predList){
  #spp-env data frame subset to important covariates with match_guilds 
  #model - one of the following to build models appropriately: gam, maxent, rf, brt, sdmtmb
  #pa_col - name of columns containing response data 
  #xy_col - names of columns containing positional data, where the first is the x(lon) and second is the y(lat)
  #month/year_col - names of columns containing months and years associated with observations 
  #metricVec/predList - vector or list of evaluation metric to use to generate weights and prediction dataframes from CVs to use to validate model
  
  #but first, remove correlated variables 
  corInd <- findCorrelation(cor(se[,-c(pa_col, xy_col)]), names = T) #find correlated variables 
  se <- se[,-which(colnames(se) == corInd)]
  print('Removing ', corInd, 'due to high colinearity...')
  
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
      form <- paste0(form, "+ s(", month_col, ", bs = 'cc', k = 6")
    }
    #loop through remaining covariates since they will all have the same smoother 
    for(x in 1:ncol(se)){
      if(colnames(se)[x] != pa_col | colnames(se)[x] != month_col | colnames(se)[x] != xy_col[1] | colnames(se)[x] != xy_col[2]){ #make sure you don't add the response variable or the variables you've already added
        form <- paste0(form, ' + s(', colnames(se)[x], ", bs = 'ts', k = 6)")
      }
    } #end for x 
    
    #run model
    mod <- FitGAM(gam.formula = formula(form), data = se, family.gam = "binomial", select = T, reduce = T)
    
    #cross validation
    cv <- CrossValidateModel(model = mod, data = se, folds = 10, model.type = "gam")
    
    #get evaluation metrics (RMSE & AUC)
    preds <-cv[[1]]
    rmse <- RMSE(obs = preds$abund, pred = preds$cvpred)
    Pred <- prediction(preds$cvpred, preds$abund)
    Perf <- performance(Pred, 'auc')
    auc <- Perf@y.values[[1]]
    
    #extract relative deviance explained
    RDE <- GAMStats(model = bi.model, data = se)
    
  } #end if gam 
  
  if(model == 'maxent'){
    print('Building MAXENT...')
    #run model - no formula needed
    mod <- FitMaxnet(data = se, species = pa_col, vars = names(se)[-which(names(se) == pa_col)], reduce = T) #fit maxent with all covariates in dataframe since we've already subset the dataframe to be only relevant covariates 
    
    #cross validate and get RMSE
    cv <- CrossValidateModel(model = mod, data = se, folds = 10, model.type = "maxnet", species = pa_col, scale.preds = F)
    preds <- cv[[1]]
    rmse <- RMSE(obs = preds$abund, pred = preds$cvpred)
    Pred <- prediction(preds$cvpred, preds$abund)
    Perf <- performance(Pred, 'auc')
    auc <- Perf@y.values[[1]]
    
    ##get relative deviance explained
    RDE <- MaxnetStats(model = mod, data = se, species = pa_col)
  } #end if maxent
  
  if(model == 'rf'){
    print('Building Random Forest with Spatial Interpolation...')
    
    se <- cbind(1:nrow(se), se) #stand in station ids 
    colnames(se)[1] <- "staid"
    se$year <- year(my(paste(se[,month_col], se[,year_col], sep = '-')))
    
    #convert dataframe to spatial object 
    stDF = st_as_sf(se, coords = xy_col, crs = 4326, agr = "constant")
    stDF = st_sftime(stDF, time_column_name = year_col)
    
    #create formula 
    form <- "value ~ "
    for(x in 1:ncol(se)){
      if(colnames(se)[x] != pa_col | colnames(se)[x] != xy_col[1] | colnames(se)[x] != xy_col[2]){ #make sure you don't add the response variable or the variables you've already added
        form <- paste0(form, ' + ', colnames(se)[x])
      }
    } #end for x 
    
    #tuning parameters for model
    n.obs <- 5:10
    min.node.size <- 2:10
    sample.fraction <- seq(1, 0.632, -0.05) # 0.632 without / 1 with replacement
    splitrule <- "variance"
    ntree <- 200 # 500
    mtry <- 3:(2+2*max(n.obs))
    tgrid = expand.grid(min.node.size=min.node.size, num.trees=ntree,
                        mtry=mtry, n.obs=n.obs, sample.fraction=sample.fraction)
    
    #tune model 
    rfsi_tune <- tune.rfsi(formula = formula(form),
                           data = stDF,
                           data.staid.x.y.z = c('staid', 'X', 'Y', 'year'),
                           cpus = 10,
                           progress = T,
                           importance = "permutation",
                           seed = 42,
                           acc.metric = 'R2',
                           tgrid = tgrid,
                           tgrid.n= length(tgrid),
                           tune.type = "LLO", # Leave-Location-Out CV
                           write.forest = T, 
                           s.crs = st_crs(stDF), 
                           k = 10, 
                           use.idw = T,
                           classification = F)
    
    #save best model 
    mod <- rfsi_tune$final.model
    
    # cross-validate (this will also run tune.rfsi, but it doesn't give you to the output like the raw code suggests! #rude)
    cv <- cv.rfsi(formula = formula(form),
                      data = stDF,
                      data.staid.x.y.z = c('staid', 'X', 'Y', 'year'),
                      cpus = 10,
                      progress = 1,
                      importance = "permutation",
                      seed = 42,
                      acc.metric = 'R2',
                      tgrid = tgrid,
                      tgrid.n= length(tgrid),
                      tune.type = "LLO", # Leave-Location-Out CV
                      write.forest = T, 
                      s.crs = st_crs(stDF), 
                      k = 10,
                      classification = F)
    
    
    #get RMSE/AUC
    rmse <- EFHSDM::RMSE(obs = cv$obs, pred = cv$pred)
    
    Pred <- prediction(cv$pred, cv$obs)
    Perf <- performance(Pred, 'auc')
    auc <- Perf@y.values[[1]]
    
    ##get important covariates 
    RDE <- importance(x=rfsi_model, method = 'altmann', formula = formula(form), data = stDF)
    
  } #end if RF
  
  if(model == 'brt'){
    print('Building Boosted Regression Trees...')
    
    modAll <- gbm.step(data = se, gbm.x = names(se)[-which(names(se) == pa_col)], 
                    gbm.y = pa_col, 
                    family = 'bernoulli', tree.complexity = 5,
                    learning.rate = 0.005, bag.fraction = 0.75, n.folds = 10) #using the same parameters as Braun et al 2023
    
    ###simplify model before k-fold validations 
    simpBRT <- gbm.simplify(modAll)
    
    #re-run model with best parameters
    mod <- gbm.step(data = se, gbm.x = simpBRT$pred.list[[length(simpBRT$pred.list)]], 
                     gbm.y = pa_col, 
                     family = 'bernoulli', tree.complexity = 5,
                     learning.rate = 0.005, bag.fraction = 0.75, n.folds = 10)
    
    #get relative importance
    RDE <- modS$contributions
    
    #k-fold validation code from Camrin Brawn (WHOI): https://zenodo.org/records/7971532
    brt_kfold <- eval_kfold_brt(dataInput = spDF2, 
                                gbm.x = simpBRT$pred.list[[length(simpBRT$pred.list)]], #only use simplified parameters for k-folds 
                                gbm.y = pa_col,  
                                learning.rate = 0.005, 
                                bag.fraction = 0.75, 
                                tree.complexity = 5, 
                                k_folds = 10,
                                is_fixed = F)
    
    ##edited eval_kfold_brt to provide all training data subsets with predictions/obs to calculate RMSE 
    # eval_kfold_brt will now return the summary_stats table as item 1 in the list, 
    # and a datatable with the observations and predictions in item 2 of the list
    #calculate RMSE
    rmse <- EFHSDM::RMSE(obs = brt_kfold[[2]]$value, pred = brt_kfold[[2]]$preds)
    #calculate AUC
    Pred <- prediction(brt_kfold[[2]]$preds, brt_kfold[[2]]$value)
    Perf <- performance(Pred, 'auc')
    auc <- Perf@y.values[[1]]
    
    ##put obs and predicted together for ensemble 
    preds <- data.frame(obs = brt_kfold[[2]]$value, preds = brt_kfold[[2]]$preds)
    
  } #end if BRT
  
  if(model == 'sdmtmb'){
    print('Building sdmTMB...')
    
    #make mesh
    mesh <- make_mesh(se, xy_cols = xy_cols, cutoff = 6/12) #using lon/lat since this is on the reprojected regular lat/lon grid, and the domain crosses multiple UTM zones  
    #MOM6 resolution is 1/12 = ~8 km 
    
    #build formula
    form <- paste0(pa_col, " ~ ")
    #loop through covariates since they will all have the same smoother 
    for(x in 1:ncol(se)){
      if(colnames(se)[x] != pa_col | colnames(se)[x] != xy_col[1] | colnames(se)[x] != xy_col[2] | colnames(se)[x] != year_col){ #make sure you don't add the response variable or the variables you've already added
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
      time = year_col,
      time_varying_type = 'ar1',
      reml = T,
      anisotropy = F,
      share_range = F,
      do_fit = T
    )
    
    #cross-validation
    cv <- sdmTMB_cv(
      formula = formula(form),
      data = se,
      mesh = mesh,
      family = binomial(link = 'logit'),
      spatial = "on",
      spatiotemporal = 'ar1',
      time = year_col,
      time_varying_type = 'ar1',
      reml = T,
      anisotropy = F,
      share_range = F,
      do_fit = T,
      k_folds = 10,
      parallel = F
    )
    
    #calculate RMSE
    rmse <- EFHSDM::RMSE(obs = cv$data$value, pred = cv$data$cv_predicted)
    
    #calculate AUC
    Pred <- prediction(cv$data$cv_predicted, cv$data$value)
    Perf <- performance(Pred, 'auc')
    auc <- Perf@y.values[[1]]
    
    #change names from sdmTMB preds to make it work with EFHSDM
    preds <- fitCV$data
    colnames(preds)[9:11] <- c('fold', 'pred', 'loglik')
    
    ### get relative importance of model using type 3 anova method
    # model with *only* intercept and no random fields:
    fit_null <- sdmTMB(
      value ~ 1, 
      spatial = "off", 
      family = binomial(link = 'logit'),
      data = spDF2,
    )
    
    #loop across variables and get their partial deviance explained 
    RDE <- vector(length = ncol(se))
    for(x in 1:ncol(se)){
      if(colnames(se)[x] != pa_col | colnames(se)[x] != xy_col[1] | colnames(se)[x] != xy_col[2] | colnames(se)[x] != year_col){ #don't do this for the response variable, xy vars, or year
        #remove one covariate 
        seSub <- se[,-x]
        
      #build formula
        formSub <- paste0(pa_col, " ~ ")
        #loop through covariates since they will all have the same smoother 
        for(x in 1:ncol(seSub)){
          if(colnames(seSub)[x] != pa_col | colnames(seSub)[x] != xy_col[1] | colnames(seSub)[x] != xy_col[2] | colnames(seSub)[x] != year_col){ #make sure you don't add the response variable or the variables you've already added
            formSub <- paste0(formSub, ' + s(', colnames(seSub)[x], ", bs = 'ts', k = 6)")
          }
        } #end for x 
      
      fitSub <- sdmTMB(
        formula = formula(formSub),
        data = seSub,
        mesh = mesh,
        family = binomial(link = 'logit'),
        spatial = "on",
        spatiotemporal = 'ar1',
        time = year_col,
        time_varying_type = 'ar1',
        reml = T,
        anisotropy = F,
        share_range = F,
        do_fit = T
      )
      
      RDE[x] <- 1 - deviance(fitSub) / deviance(fit_null)
      print(x)
      }
    } #end if
    names(RDE) <- names(spDF2)
    
  } #end if sdmtmb
  
  if(model == 'ens'){
    weights <- MakeEnsemble(rmse = metricVec) #make weights
    mod <- ValidateEnsemble(pred.list = predList, model.weights = weights, make.plots = F) #validate to get preds/obs to metrics 
    
    #evaluation metrics 
    rmse <- RMSE(obs = mod$abund, pred = mod$pred)
    Pred <- prediction(mod$pred, mod$abund)
    Perf <- performance(Pred, 'auc')
    auc <- Perf@y.values[[1]]
    
  } #end if ensemble 
  
  #change output slightly for ensemble 
  if(model != 'ens'){
     return(list(mod, cv, preds, rmse, auc, RDE))
  } else {
    return(list(mod, weights, rmse, auc))
  }
  
}