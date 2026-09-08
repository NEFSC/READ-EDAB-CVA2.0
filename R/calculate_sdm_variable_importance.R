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
    
    # 1. Create a proper Date column by appending "01." (the 1st day of the month)
    seSub$true_date <- as.Date(paste0("01.", seSub$month.year), format = "%d.%m.%Y")
    
    # 2. Convert dataframe to spatial object
    stDF = sf::st_as_sf(seSub, coords = xy_col, crs = 4326, agr = "constant")
    
    # 3. Use the new true_date column for your sftime temporal dimension
    stDF = sftime::st_sftime(stDF, time_column_name = "true_date")
    
    ##get important covariates
    # 1. Strip the spatial geometry for ranger compatibility
    se_df <- sf::st_drop_geometry(stDF)
    
    # 2. Get baseline predictions and calculate a performance metric 
    # (Assuming probability predictions for Presence/Absence)
    base_preds <- stats::predict(mod, data = se_df)$predictions
    
    # If your model outputs probabilities for classes, make sure to select the "Presence" column
    if(is.matrix(base_preds)) base_preds <- base_preds[, "1"] 
    
    # Calculate baseline performance (e.g., using a simple metric like Log Loss or Brier Score)
    # Here we use Brier Score (Mean Squared Error for probabilities, lower is better)
    base_brier <- mean((base_preds - se_df[[pa_col]])^2)
    
    importance_df <- data.frame(Variable = var_names, Importance = NA)
    
    # 3. Perform Block Permutation Importance
    for (i in seq_along(var_names)) {
      v <- var_names[i]
      perm_data <- se_df
      
      # SHUFFLE WITHIN BLOCKS: ave() applies the sample function within each sp.tm group
      perm_data[[v]] <- stats::ave(perm_data[[v]], perm_data$sp.tm, FUN = sample)
      
      # Predict on the spatially shuffled data
      perm_preds <- stats::predict(mod, data = perm_data)$predictions
      if(is.matrix(perm_preds)) perm_preds <- perm_preds[, "1"]
      
      # Calculate degraded performance
      perm_brier <- mean((perm_preds - perm_data[[pa_col]])^2)
      
      # Importance is the INCREASE in error (larger increase = more important)
      # Bounded at 0 in case random noise slightly improves the model by chance
      importance_df$Importance[i] <- max(0, perm_brier - base_brier)
    }
    
    # Normalize to percentages
    importance_df$Pct_Importance <- (importance_df$Importance / sum(importance_df$Importance)) * 100
    
    # Sort from most to least important
    RDE <- importance_df[order(-importance_df$Pct_Importance), ]
    
  } #end if RF
  
  if(model == 'brt'){
    #get relative importance
    RDE <- mod$contributions
  } #end if BRT
  
  if(model == 'sdmtmb'){
    requireNamespace("sdmTMB", quietly = TRUE)
    
    # 1. Calculate a baseline performance metric (e.g., Pseudo R-squared of the fixed effects)
    pred_baseline <- stats::predict(mod, re_form = NA)$est
    
    # Replace 'response_var' with the actual name of your dependent variable
    base_r2 <- stats::cor(pred_baseline, mod$data[,pa_col], use = "complete.obs")^2 
    
    dyn_names <- c(var_names, month_col, year_col)
    
    importance_df <- data.frame(Variable = dyn_names, Importance = NA)
    
    # 2. Shuffle each variable to see how much performance drops
    for (i in seq_along(dyn_names)) {
      v <- dyn_names[i]
      
      # Create a copy of the data and randomly shuffle the target variable
      perm_data <- mod$data
      perm_data[[v]] <- sample(perm_data[[v]])
      
      # Predict using the shuffled data
      pred_perm <- stats::predict(mod, newdata = perm_data, re_form = NA)$est
      
      # Calculate performance with the shuffled variable
      perm_r2 <- stats::cor(pred_perm, mod$data[,pa_col], use = "complete.obs")^2
      
      # Importance is the drop in R-squared (larger drop = more important)
      importance_df$Importance[i] <- max(0, base_r2 - perm_r2)
    }
    
    # Normalize to percentages
    importance_df$Pct_Importance <- (importance_df$Importance / sum(importance_df$Importance)) * 100
    RDE <- importance_df
    
  } #end if sdmtmb
  
  return(RDE)
  
}


