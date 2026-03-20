#' @title Component Model Cross-Validation
#' @description Perform cross-validation on one of the five component model types. k is 5 for all models.
#'
#' @param mod the output from \code{make_sdm}
#' @param se data frame containing species presence/absence data and desired environmental covariate data.
#' @param pa_col column name for presence/absence column
#' @param xy_col a vector with a length of 2 indicating the longitude and latitude column names
#' @param month_col,year_col column names for month and year columns respectively
#' @param model one of the following indicating the desired model to perform cross validation on: gam, maxent, brt, rf, or sdmtmb
#'
#' @return the cross-validation results from the desired model. Objects will differ slightly depending on model type.


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
    stDF = sf::st_as_sf(seSub[complete.cases(seSub),], coords = xy_col, crs = 4326, agr = "constant")
    stDF = sftime::st_sftime(stDF, time_column_name = month_col)


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
    cv <- meteo::cv.rfsi(formula = formula(form),
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
                  s.crs = sf::st_crs(stDF),
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
      mesh <- sdmTMB::make_mesh(train, xy_cols = xy_col, cutoff = 1) #using lon/lat since this is on the reprojected regular lat/lon grid, and the domain crosses multiple UTM zones
      #MOM6 resolution is 1/12 = ~8 km
      tryMod <- tryCatch(expr = {
        #create model
        cv <- sdmTMB::sdmTMB(
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

