#' @title Predict Component or Ensemble SDM
#' @description Make model predictions for either a component model (GAM/MAXENT/RF/BRT/SDMTMB) or the Ensemble model
#'
#' @param mod the output from \code{make_sdm} - only used for ensemble
#' @param model one of the following indicating the desired model to calculate variable importance for: gam, maxent, brt, rf, sdmtmb, or ens
#' @param rasts for component models (GAM/MAXENT/RF/BRT/SDMTMB) list of spatRasters corresponding to the environmental covariates used to build the models. The number of layers in each spatRaster should be the same and correspond to the length of the timeseries for the models to be predicted on. For the Ensemble, a list of spatRasters containing the predicted values for each of the component models. The length of the list should be equal to the length of \code{weights}
#' @param static_variables spatRaster containing the static variables used in model. 
#' @param se data frame containing species presence/absence data and desired environmental covariate data.
#' @param xy_col a vector with a length of 2 indicating the longitude and latitude column names
#' @param month_col,year_col column names for month and year columns respectively
#' @param pa_col column name for presence/absence column
#' @param weights a vector of model weights - used for building the ensemble model
#'
#' @return returns a rasterStack of predicted habitat suitability. The number of layers will be equal to the number of layers in \code{rasts}
#'
#'@export

# --- MAIN SDM PREDICTION FUNCTION ---
make_sdm_predictions <- function(
    mod, model, rasts, static_variables, se = NULL, 
    month_col, pa_col, year_col, xy_col, weights = NULL
) {
  
  # Fail fast check
  valid_models <- c('gam', 'maxent', 'rf', 'brt', 'sdmtmb', 'ensemble')
  if (!model %in% valid_models) stop("Model must be one of: ", paste(valid_models, collapse = ", "))
  
  # Handle Ensemble separately as it doesn't need covariate stacks
  if (model == 'ensemble') {
    if (length(rasts) != length(weights)) stop('raster list and weights are not the same length')
    print('Predicting Ensemble...')
    hsm <- vector(mode = 'list', length = length(rasts[[1]]))
    
    for (x in seq_along(rasts[[1]])) {
      abunds <- lapply(rasts, function(r_list) terra::rast(r_list[[x]]))
      hsm[[x]] <- terra::rast(EFHSDM::MakeEnsembleAbundance(model.weights = weights, abund.list = abunds))
    }
    names(hsm) <- names(rasts[[1]])
    return(terra::rast(hsm))
  }
  
  # --- COMMON SETUP FOR SPATIAL MODELS USING TERRA ---
  # Convert inputs to terra pointers for fast coordinate generation
  template_r <- rasts[[1]][[1]] 
  
  # Get exact coordinates per cell efficiently
  xy_mat <- terra::xyFromCell(template_r, 1:terra::ncell(template_r))
  rlon   <- terra::rast(template_r, vals = xy_mat[, 1])
  rlat   <- terra::rast(template_r, vals = xy_mat[, 2])
  
  num_layers <- terra::nlyr(rasts[[1]])
  hsm <- vector(mode = 'list', length = num_layers)
  
  # --- SPECIAL PRE-PROCESSING: RANDOM FOREST DATA SUBSAMPLING ---
  if (model == 'rf') {
    print('Processing RFSI space-time data framing...')
    se <- cbind(1:nrow(se), se)
    colnames(se)[1] <- "staid"
    
    # Stratified Spatial Regions Setup
    se$region <- NA
    se$region[which(se[, xy_col[1]] > -70 & se[, xy_col[2]] < 41.5)] <- 'GB'
    se$region[which(se[, xy_col[1]] > -71 & se[, xy_col[2]] > 41.5)] <- 'GOM'
    se$region[which(se[, xy_col[1]] < -70 & se[, xy_col[2]] < 42 & se[, xy_col[2]] > 39.5)] <- 'SNE'
    se$region[which(se[, xy_col[2]] < 39.5)] <- 'MAB'
    
    se$sp.tm <- paste(se$month.year, se$region, sep = '-')
    sptm <- unique(se$sp.tm)
    
    set.seed(2025)
    seSub <- do.call(rbind, lapply(sptm, function(x) {
      sub  <- se[se$sp.tm == x, ]
      abs  <- sub[sub[,pa_col] == 0, ]
      pres <- sub[sub[,pa_col] == 1, ]
      
      if (nrow(pres) <= 5) {
        return(rbind(abs[sample(nrow(abs), round(nrow(abs) / 4)), ], pres))
      } else if (nrow(abs) > nrow(pres)) {
        return(rbind(abs[sample(nrow(abs), nrow(pres)), ], pres))
      } else {
        return(sub)
      }
    }))
    
    stDF <- sf::st_as_sf(seSub, coords = xy_col, crs = 4326, agr = "constant")
    stDF <- sftime::st_sftime(stDF, time_column_name = month_col)
  }
  
  if (model == 'sdmtmb') {
    warning('sdmTMB predictions can take a long time (approximately 1 min per timestamp).')
  }
  
  # --- MAIN TIME-SERIES LOOP (CONSOLIDATED) ---
  message("Predicting ", toupper(model), " model across time steps...")
  
  for (x in 1:num_layers) {
    # Run our helper function to build the stack
    sr_terra <- prep_time_step_stack(x, rasts, template_r, rlon, rlat, static_variables, month_col, year_col, xy_col)
    
    # Model-specific executions
    if (model == 'gam') {
      # Convert each layer to a standard legacy RasterLayer first, then stack them
      sr_legacy <- raster::stack(lapply(1:terra::nlyr(sr_terra), function(i) {
        raster::raster(sr_terra[[i]])
      }))
      hsm[[x]]  <- terra::rast(EFHSDM::MakeGAMAbundance(model = mod, r.stack = sr_legacy))
      
    } else if (model == 'maxent') {
      requireNamespace("maxnet", quietly = TRUE) # Forces R to load maxnet and register all its S3 methods (like predict.maxnet)
      # Convert each layer to a standard legacy RasterLayer first, then stack them
      sr_legacy <- raster::stack(lapply(1:terra::nlyr(sr_terra), function(i) {
        raster::raster(sr_terra[[i]])
      }))
      hsm[[x]]  <- EFHSDM::MakeMaxEntAbundance(model = mod, maxent.stack = sr_legacy, type = 'maxnet')
      
    } else if (model == 'brt') {
      sr_df <- prep_time_step_df(x, rasts, static_variables)
      colnames(sr_df)[1:2] <- xy_col
      p <- gbm::predict.gbm(newdata = sr_df, object = mod, type = "response")
      p_df <- cbind(sr_df[,1:2], p)
      # hsm[[x]] <- terra::rasterize(as.matrix(p_df[, 1:2]), template_r, field = p_df$p) 
      
      r_pred <- terra::rast(p_df, type = "xyz")
      terra::crs(r_pred) <- terra::crs(template_r)
      r_pred_extended    <- terra::extend(r_pred, template_r) # Ensures alignment
      hsm[[x]] <- r_pred_extended
      
    } else if (model == 'sdmtmb') {
      requireNamespace("sdmTMB", quietly = TRUE) # Forces R to load maxnet and register all its S3 methods (like predict.maxnet)
      
      sr_df <- prep_time_step_df(x, rasts, static_variables)
      colnames(sr_df)[1:2] <- xy_col
      
      p  <- stats::predict(mod, newdata = sr_df, type = 'response')
      
      p_df <- cbind(sr_df[,1:2], p$est)
      r_pred <- terra::rast(p_df, type = "xyz")
      terra::crs(r_pred) <- terra::crs(template_r)
      hsm[[x]]   <- terra::extend(r_pred, template_r) # Ensures alignment
      
    } else if (model == 'rf') {
      sr_df <- prep_time_step_df(x, rasts, static_variables)
      sr_df$staid <- 1:nrow(sr_df)
      
      sr_sf <- sf::st_as_sf(sr_df, coords = c('x', 'y'), crs = 4326, agr = "constant") |> 
        sftime::st_sftime(time_column_name = month_col)
      
      predDF <- meteo::pred.rfsi(
        model = mod, data = stDF, data.staid.x.y.z = c('staid', 'x', 'y'), obs.col = pa_col,
        newdata = sr_sf, newdata.staid.x.y.z = c('staid', 'x', 'y'),
        output.format = "data.frame", cpus = 1, progress = FALSE, classification = FALSE
      )
      
      # Fast spatial interpolation (Kriging step) using clean data frame mapping
      pts_df <- predDF[, c('pred', 'X', 'Y')]
      colnames(pts_df) <- c("val", "x", "y")
      
      v_geo <- gstat::variogram(val ~ 1, loc = ~ x + y, data = pts_df)
      v_fit_anis <- gstat::fit.variogram(v_geo, model = gstat::vgm(psill=0.06, model="Sph", range=3, nugget=0.01, anis=c(90, 0.66)))
      
      grid_pts <- as.data.frame(template_r, xy = TRUE, na.rm = FALSE)
      sp::coordinates(grid_pts) <- ~ x + y
      sp::coordinates(pts_df)  <- ~ x + y
      
      kriged_df <- as.data.frame(gstat::krige(val ~ 1, pts_df, grid_pts, model = v_fit_anis, nmax = 50))
      
      r_filled <- terra::rast(kriged_df[, c("x", "y", "var1.pred")], type = "xyz")
      terra::crs(r_filled) <- terra::crs(template_r)
      hsm[[x]] <- r_filled
    }
  }
  
  hsm <- terra::rast(hsm)
  names(hsm) <- names(rasts[[1]])
  return(hsm)
}