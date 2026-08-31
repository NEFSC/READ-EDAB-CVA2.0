#' @title Build Component or Ensemble Species Distribution Models
#' @description Build one of the five component models or the final ensemble model
#'
#' @param se data frame containing species presence/absence data and desired environmental covariate data.
#' @param pa_col column name for presence/absence column
#' @param xy_col a vector with a length of 2 indicating the longitude and latitude column names
#' @param month_col,year_col column names for month and year columns respectively. Defaults to 'month' and 'year' respectively.
#' @param var_names a vector of covariate names to use in the desired model. Should match some or all of the column names in \code{se}.
#' @param model one of the following indicating the desired model to build: gam, maxent, brt, rf, sdmtmb, or ens
#' @param year_range vector with length of two including the maximum and minimum years to include in the model. Only used if \code{model = 'sdmtmb'}. Should include maximum desired year if predicting as well; i.e. if training data range from 1993-2019, but you want to forecast out to 2035, the \code{year_range} should be c(1993,2035). Defaults to the range of the \code{year_col} in \code{se}
#' @param ensemble_weights vector of model weights used to build ensemble model. The vector must have the same length and be in the same order as the corresponding predictions list.
#' @param ensemble_preds a list of prediction values from component models used to build ensemble model. The list of predictions must have the same length and be in the same order as the corresponding weight vector.
#'
#' @return the model object from the desired model. Model objects will differ depending on the model type.
#'
#'@export

build_sdm <- function(
    se,
    pa_col,
    xy_col,
    month_col = 'month',
    year_col = 'year',
    var_names,
    model,
    year_range = range(se[,year_col], na.rm = T),
    ensemble_weights = NULL,
    ensemble_preds = NULL
) {
  # Fail fast check
  valid_models <- c("gam", "maxent", "rf", "brt", "sdmtmb", "ensemble")
  if (!model %in% valid_models) stop("Model must be one of: ", paste(valid_models, collapse = ", "))


  #now build model of choice
  if (model == "gam") {
    #build gam model
    print('Building GAM...')
    #build formula
    form <- paste0(pa_col, " ~ ")
    if (!is.null(xy_col)) {
      form <- paste0(
        form,
        's(',
        xy_col[1],
        ',',
        xy_col[2],
        ", bs = 'ts', k = 10)"
      )
    }
    if (!is.null(month_col)) {
      form <- paste0(form, "+ s(", month_col, ", bs = 'cc', k = 6)")
    }
    if (!is.null(year_col)) {
      form <- paste0(form, "+ s(", year_col, ", bs = 'ts', k = 6)")
    }
    #loop through covariates since they will all have the same smoother
    for (x in var_names) {
      form <- paste0(form, ' + s(', x, ", bs = 'ts', k = 6)")
    } #end for x

    #run model
    mod <- EFHSDM::FitGAM(
      gam.formula = stats::formula(form),
      data = se,
      family.gam = "binomial",
      select = T,
      reduce = T
    )
  } #end if gam

  if (model == "maxent") {
    print('Building MAXENT...')
    #run model - no formula needed
    #need to remove extraneous variables
    seSub <- se[,names(se) %in% c(pa_col, xy_col, month_col, year_col, var_names)]
    mod <- EFHSDM::FitMaxnet(
      data = seSub,
      species = pa_col,
      vars = names(seSub)[-which(names(seSub) == pa_col)],
      reduce = T
    ) #fit maxent with all covariates in dataframe since we've already subset the dataframe to be only relevant covariates
  } #end if maxent

  if (model == "rf") {
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

    #make clean coordinate columns to help with cv
    coords <- sf::st_coordinates(stDF)
    stDF$X <- coords[, "X"]
    stDF$Y <- coords[, "Y"]

    #create formula
    # Build formula string
    form <- paste(pa_col, "~", paste(var_names, collapse = " + "))

    #build model
    mod <- meteo::rfsi(
      formula = stats::formula(form),
      data = stDF,
      data.staid.x.y.z = c('staid', 'X', 'Y'),
      cpus = 1,
      progress = F,
      importance = "impurity",
      seed = 42,
      num.trees = 200,
      s.crs = sf::st_crs(stDF),
      use.idw = T,
      write.forest = T,
      splitrule = "extratrees",
      min.node.size = 5,
      sample.fraction = 0.95
    )

    #reduce parameters
    history <- list()
    current_predictors <- var_names

    # Create spatial LLO (Leave-Location-Out) folds manually based on your unique stations
    unique_stations <- unique(stDF$staid)
    set.seed(42)
    num_folds <- 5
    station_folds <- split(sample(unique_stations), rep(1:num_folds, length.out = length(unique_stations)))

    repeat {
      # Build formula string
      form_string <- paste(pa_col, "~", paste(current_predictors, collapse = " + "))
      current_formula <- stats::formula(form_string)

      # 2. Spatial Cross-Validation with Probability output
      # 'probability = TRUE' forces ranger to output probability vectors instead of hard 0/1s
      # ... Inside your repeat/loop block ...

      # Vectors to pool predictions across folds
      all_observed <- c()
      all_predicted_probs <- c()

      # --- 2. Manual Cross-Validation Loop ---
      for (f in 1:num_folds) {
        val_stations <- station_folds[[f]]

        # Split data based on spatial station IDs
        train_df <- stDF[!stDF$staid %in% val_stations, ]
        val_df   <- stDF[stDF$staid %in% val_stations, ]

        # Fit the RFSI model using the exact parameters that work for you
        # Note: probability = TRUE must be supplied via ranger arguments (...)
        model_fit <- tryCatch({
          meteo::rfsi(
            formula = current_formula,
            data = train_df,
            data.staid.x.y.z = c('staid', 'X', 'Y'),
            cpus = 1,
            progress = FALSE,
            num.trees = 200,
            s.crs = sf::st_crs(stDF),
            use.idw = TRUE,
            splitrule = "extratrees",        # Gini for classification/probability
            min.node.size = 5,
            sample.fraction = 0.95,
            probability = TRUE,        # Forces probability forest execution
            importance = "impurity"
          )
        }, error = function(e) NULL)

        if (is.null(model_fit)) next

        # Generate predictions on the validation fold
        # pred.rfsi returns probability vectors when given a probability model
        predictions <- meteo::pred.rfsi(
          model = model_fit,
          data = train_df,             # Conditioning data points
          obs.col = pa_col,
          data.staid.x.y.z = c('staid', 'X', 'Y'),
          newdata = val_df,            # Locations to predict
          newdata.staid.x.y.z = c('staid', 'X', 'Y'),
          s.crs = sf::st_crs(stDF),
          newdata.s.crs = sf::st_crs(stDF),
          progress = FALSE
        )

        # 1. Convert val_df to a plain data frame to strip spatial geometry for a clean merge
        val_plain <- as.data.frame(val_df)

        # This matches 'staid' and 'time', plus coordinates if named exactly the same
        merge_keys <- intersect(c("staid", "time", "X", "Y"), names(predictions))

        # 3. Merge them together. This keeps ONLY rows successfully predicted by pred.rfsi
        aligned_results <- merge(val_plain, predictions, by.x = c("staid", month_col, "X", "Y"), by.y = merge_keys)

        # 4. Safely append the perfectly paired observed outcomes and predicted probabilities
        all_observed        <- c(all_observed, aligned_results[[pa_col]])
        all_predicted_probs <- c(all_predicted_probs, aligned_results[["2"]])
      }

      # --- 3. Compute AUC Score ---
      roc_obj <- pROC::roc(all_observed, all_predicted_probs, quiet = TRUE)
      current_auc <- as.numeric(pROC::auc(roc_obj))

      cat(sprintf("Variables (%d): %s | CV AUC: %.4f\n",
                  length(current_predictors),
                  paste(current_predictors, collapse = ", "),
                  current_auc))

      # Save iteration details
      history[[length(current_predictors)]] <- list(vars = current_predictors, auc = current_auc)

      # Base Case: Stop if only 1 environmental predictor remains
      if (length(current_predictors) == 1) { break }

      # --- 4. Identify Weakest Variable to Drop ---
      # Fit once on full dataset to get final importance
      global_fit <- meteo::rfsi(
        formula = current_formula, data = stDF,
        data.staid.x.y.z = c('staid', 'X', 'Y'), cpus = 1, progress = FALSE,
        num.trees = 200, s.crs = sf::st_crs(stDF), use.idw = TRUE,
        splitrule = "extratrees", min.node.size = 5, sample.fraction = 0.95,
        probability = TRUE, importance = "impurity"
      )

      imp_scores <- global_fit$variable.importance[current_predictors]
      weakest_var <- names(which.min(imp_scores))

      # Remove the weakest link
      current_predictors <- setdiff(current_predictors, weakest_var)
    }

    # --- 5. Print out the Optimal Parsimonious Model Matrix ---
    history_df <- do.call(rbind, lapply(history, function(x) {
      if(is.null(x)) return(NULL)
      data.frame(num_vars = length(x$vars), auc = x$auc, vars = paste(x$vars, collapse=", "))
    }))

    # Parse out the 1% parsimony option we implemented earlier
    max_auc <- max(history_df$auc)
    parsimonious_row <- history_df[history_df$auc >= (max_auc * 0.99), ]
    best_row <- parsimonious_row[which.min(parsimonious_row$num_vars), ]

    cat("\n--- Manual Optimization Complete ---\n")
    cat("Selected Parsimonious Model:\n")
    cat("Number of Predictors:", best_row$num_vars, "\n")
    cat("Parsimonious CV AUC :", round(best_row$auc, 4), "\n")
    cat("Selected Variables  :", best_row$vars, "\n")

    ##build final model
    form_string <- stats::formula(paste(pa_col, "~", paste(strsplit(best_row$vars, ', ')[[1]], collapse = " + ")))

    mod <- meteo::rfsi(
      formula = stats::formula(form_string),
      data = stDF,
      data.staid.x.y.z = c('staid', 'X', 'Y'),
      cpus = 1,
      progress = F,
      importance = "impurity",
      seed = 42,
      num.trees = 200,
      s.crs = sf::st_crs(stDF),
      use.idw = T,
      write.forest = T,
      splitrule = "extratrees",
      min.node.size = 5,
      sample.fraction = 0.95
    )

  } #end if RF

  if (model == "brt") {
    print('Building Boosted Regression Trees...')

    #need to remove extraneous variables
    seSub <- se[,names(se) %in% c(pa_col, xy_col, month_col, year_col, var_names)]

    modAll <- dismo::gbm.step(
      data = seSub,
      gbm.x = names(seSub)[-which(names(seSub) == pa_col)],
      gbm.y = pa_col,
      family = 'bernoulli',
      tree.complexity = 5,
      learning.rate = 0.005,
      bag.fraction = 0.75,
      max.trees = 2000,
      n.folds = 10
    ) #using the same parameters as Braun et al 2023

    ###simplify model before k-fold validations
    simpBRT <- dismo::gbm.simplify(modAll)

    #re-run model with best parameters
    mod <- dismo::gbm.step(
      data = seSub,
      gbm.x = simpBRT$pred.list[[length(simpBRT$pred.list)]],
      gbm.y = pa_col,
      family = 'bernoulli',
      tree.complexity = 5,
      learning.rate = 0.005,
      bag.fraction = 0.75,
      max.trees = 2000,
      n.folds = 10
    )
  } #end if BRT

  if (model == "sdmtmb") {
    print('Building sdmTMB...')

    se <- se[stats::complete.cases(se),]

    #build formula
    form <- paste0(pa_col, " ~ ")
    # Keep seasonal anchors as smoothers if desired (common practice)
    if (!is.null(month_col)) {
      form <- paste0(form, " s(", month_col, ", k = 4)") # Lowered k slightly for stability
    }
    if (!is.null(year_col)) {
      form <- paste0(form, " + s(", year_col, ", k = 4)")
    }

    # LOOP UPDATE: Add environmental covariates as strictly LINEAR effects
    for (x in var_names) {
      form <- paste0(form, " + ", x) # <-- No more s() or k = 6!
    }


    #make mesh
    mesh <- sdmTMB::make_mesh(se, xy_cols = xy_col, cutoff = 1) #using lon/lat since this is on the reprojected regular lat/lon grid, and the domain crosses multiple UTM zones
    #MOM6 resolution is 1/12 = ~8 km

    # --- SECTION 1: Try Initial Global Model ---
    mod <- tryCatch(
      expr = {
        mod <- sdmTMB::sdmTMB(
          formula        = stats::formula(form),
          data           = se,
          mesh           = mesh,
          family         = stats::binomial(link = 'logit'),
          spatiotemporal = 'ar1',
          time           = year_col,
          reml           = FALSE, # ML for fixed-effect AIC comparison
          anisotropy     = TRUE,
          share_range    = TRUE,
          do_fit         = TRUE,
          extra_time     = year_range[1]:year_range[2]
        )
        # Check gradient and presence of NA standard errors
        max_grad  <- max(abs(mod$gradients), na.rm = TRUE)
        fe_tidy   <- tryCatch(broom::tidy(mod, effects = "fixed"), error = function(e) NULL)
        has_na_se <- if (!is.null(fe_tidy)) any(is.na(fe_tidy$std.error)) else TRUE
        
        if (max_grad > 0.001 || has_na_se) {
            sdmTMB::run_extra_optimization(mod, nlminb_loops = 1, newton_steps = 1)
        }
        
      },
      error = function(e) {
        message('Initial model did not converge')
        return(NA)
      }
    )

    # --- SECTION 2: Automated Fast AIC Reduction ---
    if (exists('mod') && inherits(mod, 'sdmTMB')) {

      try2simp <- tryCatch(
        expr = {
          print('Simplifying model using fast AIC evaluations...')

          # Step 2a: Get baseline AIC from the initial global model
          # Using stats::AIC() is standard and instantaneous
          best_aic <- stats::AIC(mod)

          simplifying <- TRUE

          while (simplifying && length(var_names) > 1) {
            # Step 2b: Identify the weakest fixed-effect or smooth term link
            fe_summary <- broom::tidy(mod, effects = "fixed")
            env_summary <- fe_summary[fe_summary$term %in% var_names, ]


            # Standard Z-score selection for valid/mixed SEs
            env_summary$z_stat <- abs(env_summary$estimate / env_summary$std.error)
              
              weakest_var <- env_summary$term[which.min(env_summary$z_stat)]

            test_vars <- setdiff(var_names, weakest_var)

            # Step 2c: Build the candidate test formula string with linear effects
            test_form_str <- paste0(pa_col, " ~ ")
            if (!is.null(month_col)) {
              test_form_str <- paste0(test_form_str, " + s(", month_col, ", k = 4)")
            }
            if (!is.null(year_col)) {
              test_form_str <- paste0(test_form_str, " + s(", year_col, ", k = 4)")
            }

            # Append the remaining test covariates linearly
            for (x in test_vars) {
              test_form_str <- paste0(test_form_str, " + ", x) # <-- Strict linear integration
            }
            test_formula <- stats::formula(test_form_str)

            print(paste('Testing removal of:', weakest_var))

            # Step 2d: Run the candidate model ONCE (no cross-validation loops!)
            test_fit <- sdmTMB::sdmTMB(
              formula = test_formula,
              data = se,
              mesh = mesh,          # Keep original high-resolution mesh since speed is no longer an issue
              family = stats::binomial(link = 'logit'),
              spatiotemporal = 'ar1',
              time = year_col,
              reml = FALSE,
              anisotropy = TRUE,
              share_range = TRUE,
              do_fit = TRUE,
              extra_time = year_range[1]:year_range[2]
            )

            # Calculate the candidate model's AIC
            test_aic <- stats::AIC(test_fit)

            # Information Theory Rule: Accept the drop if the AIC stays lower,
            # flat, or increases by less than 2 points (AIC tolerance rule-of-thumb).
            if (test_aic <= (best_aic + 2)) {
              print(paste('Successfully removed:', weakest_var, "| New AIC:", round(test_aic, 2)))

              # Update tracking metrics and permanent variables
              best_aic <- test_aic
              var_names <- test_vars
              form <- test_form_str

              # Save this candidate as our current champion
              mod <- test_fit
            } else {
              print(paste('Drop rejected. AIC spiked too high for:', weakest_var))
              simplifying <- FALSE # Stop reducing if dropping this variable hurts the model
            }
          }
        },
        error = function(e) {
          message('Model could not be simplified further due to an internal optimization error.')
          return(NA)
        }
      )
    } else {
      mod <- NA
    }
    
    # --- SECTION 3: Final Champion Re-Fit (REML = TRUE) ---
    if (!is.null(mod) && inherits(mod, "sdmTMB")) {
      
      print('Re-fitting final champion model with REML = TRUE for optimal spatial variance estimation...')
      
      # Option A: Using stats::update (Fastest & standard in R)
      final_mod <- tryCatch(
        expr = {
          mod <- stats::update(mod, reml = TRUE)
          # Check gradient and presence of NA standard errors
          max_grad  <- max(abs(mod$gradients), na.rm = TRUE)
          fe_tidy   <- tryCatch(broom::tidy(mod, effects = "fixed"), error = function(e) NULL)
          has_na_se <- if (!is.null(fe_tidy)) any(is.na(fe_tidy$std.error)) else TRUE
          
          if (max_grad > 0.001 || has_na_se) {
            sdmTMB::run_extra_optimization(mod, nlminb_loops = 1, newton_steps = 1)
          }
        },
        error = function(e) {
          message('REML update failed. Falling back to explicit sdmTMB fit...')
          NULL
        }
      )
      
      # Option B: Fallback explicit call if update() fails
      if (is.null(final_mod)) {
        final_mod <- tryCatch(
          expr = {
            #build formula
            form <- paste0(pa_col, " ~ ")
            # Keep seasonal anchors as smoothers if desired (common practice)
            if (!is.null(month_col)) {
              form <- paste0(form, " s(", month_col, ", k = 4)") # Lowered k slightly for stability
            }
            if (!is.null(year_col)) {
              form <- paste0(form, " + s(", year_col, ", k = 4)")
            }
            
            # LOOP UPDATE: Add environmental covariates as strictly LINEAR effects
            for (x in var_names) {
              form <- paste0(form, " + ", x) # <-- No more s() or k = 6!
            }
            
            mod <- sdmTMB::sdmTMB(
              formula        = stats::formula(form),
              data           = se,
              mesh           = mesh,
              family         = stats::binomial(link = 'logit'),
              spatiotemporal = 'ar1',
              time           = year_col,
              reml           = TRUE, # Final model fitted with REML
              anisotropy     = TRUE,
              share_range    = TRUE,
              do_fit         = TRUE,
              extra_time     = year_range[1]:year_range[2]
            )
            # Check gradient and presence of NA standard errors
            max_grad  <- max(abs(mod$gradients), na.rm = TRUE)
            fe_tidy   <- tryCatch(broom::tidy(mod, effects = "fixed"), error = function(e) NULL)
            has_na_se <- if (!is.null(fe_tidy)) any(is.na(fe_tidy$std.error)) else TRUE
            
            if (max_grad > 0.001 || has_na_se) {
              sdmTMB::run_extra_optimization(mod, nlminb_loops = 1, newton_steps = 1)
            }
          },
          error = function(e) {
            message('Final REML fit failed completely. Returning ML champion model instead.')
            return(mod)
          }
        )
      }
      
      # Assign final model object
      mod <- final_mod
      print('Final model fit complete.')
      
    } else {
      message('No valid model was produced.')
      mod <- NA
    }

  } #end if sdmtmb

  if (model == "ensemble") {
    # weights <- MakeEnsemble(rmse = mets) #make weights
    mod <- EFHSDM::ValidateEnsemble(
      pred.list = ensemble_preds,
      model.weights = ensemble_weights,
      make.plots = F,
      latlon = T
    ) #validate to get preds/obs to metrics
  } #end if ensemble

  return(mod)
}
