#' @title Calculate Variable Importance
#' @description Calculate variable importance for one of the component models
#'
#' @param mod the output from \code{make_sdm} - only used for ensemble
#' @param se data frame containing species presence/absence data and desired environmental covariate data.
#' @param pa_col column name for presence/absence column
#' @param xy_col a vector with a length of 2 indicating the longitude and latitude column names
#' @param month_col,year_col column names for month and year columns respectively. Defaults to 'month' and 'year' respectively. 
#' @param var_names a vector of covariate names to use in the desired model. Should match some or all of the column names in \code{se}. 
#' @param model one of the following indicating the desired model to calculate variable importance for: gam, maxent, brt, rf, or sdmtmb
#'
#' @return a vector of the variable importance for the given model. Each model calculates these differently, so the values should be normalized in order to compare across models.

calculate_sdm_variable_importance <- function(mod, 
                                              se, 
                                              pa_col, 
                                              xy_col, 
                                              month_col = 'month', 
                                              year_col = 'year', 
                                              var_names,
                                              model){
  
  if(model == 'gam'){ #build gam model
    #extract relative deviance explained
    RDE <- EFHSDM::GAMStats(model = mod, data = se)
    
  } #end if gam
  
  if(model == 'maxent'){
    ##get relative deviance explained
    RDE <- EFHSDM::MaxnetStats(model = mod, data = se, species = pa_col)
  } #end if maxent
  
  if(model == 'rf'){
    print('Building Random Forest with Spatial Interpolation...')
    
    se <- cbind(1:nrow(se), se) #stand in station ids
    colnames(se)[1] <- "staid"
    
    #subsample by space-time
    set.seed(2025)
    
    #make regions
    se$region <- NA
    se$region[which(se[, xy_col[1]] > -70 & se[, xy_col[2]] < 41.5)] <- 'GB' #georges bank
    se$region[which(se[, xy_col[1]] > -71 & se[, xy_col[2]] > 41.5)] <- 'GOM' #gulf of maine
    se$region[which(
      se[, xy_col[1]] < -70 & se[, xy_col[2]] < 42 & se[, xy_col[2]] > 39.5
    )] <- 'SNE' #southern new england
    se$region[which(se[, xy_col[2]] < 39.5)] <- 'MAB' #mid-atlantic bight
    
    #make space-time id
    se$sp.tm <- paste(se$month.year, se$region, sep = '-')
    sptm <- unique(se$sp.tm)
    #subsample data
    seSub <- NULL
    for (x in sptm) {
      sub <- se[se$sp.tm == x, ]
      
      abs <- sub[sub[,pa_col] == 0, ]
      pres <- sub[sub[,pa_col] == 1, ]
      
      if (nrow(pres) <= 5) {
        #if there are few presences
        absSub <- abs[sample(x = nrow(abs), size = round(nrow(abs) / 4)), ] #subsample absences to a 1/4 of the absences within month and region
        allSub <- rbind(absSub, pres) #combine with presences (if any are absent)
        #this will allow all regions, years, and months to be present in the final time series to help predictions while also making the ratio of presences/absences somewhat more even
      } else if (nrow(abs) > nrow(pres)) {
        #if there are enough presences, but absences still outnumber presences
        absSub <- abs[sample(x = nrow(abs), size = nrow(pres)), ] #subsample absences
        allSub <- rbind(absSub, pres)
      } else {
        #if presences outnumber absences
        allSub <- sub #do nothing and keep it all
      }
      
      seSub <- rbind(seSub, allSub)
    }
    
    #convert dataframe to spatial object
    stDF = sf::st_as_sf(seSub, coords = xy_col, crs = 4326, agr = "constant")
    stDF = sftime::st_sftime(stDF, time_column_name = month_col)
    
    #create formula
    form <- paste0(pa_col, " ~ ")
    #add covariates - don't need to add space/time since they are already accounted for in spatial object
    for (x in var_names) {
      form <- paste0(form, ' + ', x)
    } #end for x
    
    ##get important covariates
    RDE <- ranger::importance(x=mod, method = 'altmann', formula = formula(form), data = stDF)
    
  } #end if RF
  
  if(model == 'brt'){
    #get relative importance
    RDE <- mod$contributions
  } #end if BRT
  
  if(model == 'sdmtmb'){
    # 1. Predict the fixed-effects component only (setting spatial fields to 0)
    # This isolates environmental signals from spatial absorption
    requireNamespace("sdmTMB", quietly = TRUE) # Forces R to load maxnet and register all its S3 methods (like predict.maxnet)
    pred_fixed <- predict(mod, re_form = NA)
    
    # Total variance explained by all environmental variables combined
    total_fixed_var <- var(pred_fixed$est)
    
    importance_df <- data.frame(Variable = var_names, Var_Contribution = NA, Pct_Importance = NA)
    
    # 2. Drop each variable's prediction contribution to see what is lost
    fe_coefs <- broom::tidy(mod, effects = "fixed")
    
    for (i in seq_along(var_names)) {
      v <- var_names[i]
      coef_val <- fe_coefs$estimate[fe_coefs$term == v]
      
      # Calculate what the fixed prediction would look like WITHOUT this variable
      # (Subtracting its linear effect: Beta * X)
      isolated_pred <- pred_fixed$est - (coef_val * pred_fixed[[v]])
      
      # Importance = Total environmental variance minus variance without this variable
      dropped_var <- total_fixed_var - var(isolated_pred)
      importance_df$Var_Contribution[i] <- max(0, dropped_var) # Bound at 0
    }
    
    # Normalize to percentages
    importance_df$Pct_Importance <- (importance_df$Var_Contribution / sum(importance_df$Var_Contribution)) * 100
    RDE <- importance_df
    
  } #end if sdmtmb
  
  return(RDE)
  
}


