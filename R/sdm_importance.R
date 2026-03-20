#' @title Calculate Variable Importance
#' @description Calculate variable importance for one of the component models
#'
#' @param mod the output from \code{make_sdm} - only used for ensemble
#' @param se data frame containing species presence/absence data and desired environmental covariate data.
#' @param pa_col column name for presence/absence column
#' @param xy_col a vector with a length of 2 indicating the longitude and latitude column names
#' @param month_col,year_col column names for month and year columns respectively
#' @param model one of the following indicating the desired model to calculate variable importance for: gam, maxent, brt, rf, or sdmtmb
#'
#' @return a vector of the variable importance for the given model. Each model calculates these differently, so the values should be normalized in order to compare across models.

sdm_importance <- function(mod, se, pa_col, xy_col, month_col, year_col, model){

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
    stDF = sf::st_as_sf(seSub, coords = xy_col, crs = 4326, agr = "constant")
    stDF = sftime::st_sftime(stDF, time_column_name = year_col)


    #create formula
    form <- "value ~ "
    for(x in 2:ncol(se)){
      if(colnames(se)[x] != pa_col & colnames(se)[x] != xy_col[1] & colnames(se)[x] != xy_col[2] & colnames(se)[x] != year_col &
         colnames(se)[x] != 'month.year' & colnames(se)[x] != 'region' & colnames(se)[x] != 'sp.tm' & colnames(se)[x] != 'staid'){ #make sure you don't add the response variable or the variables you've already added
        form <- paste0(form, ' + ', colnames(se)[x])
      }
    } #end for x

    ##get important covariates
    RDE <- ranger::importance(x=mod, method = 'altmann', formula = formula(form), data = stDF)

  } #end if RF

  if(model == 'brt'){
    #get relative importance
    RDE <- mod$contributions
  } #end if BRT

  if(model == 'sdmtmb'){
    se <- se[complete.cases(se),]

    #make mesh
    mesh <- sdmTMB::make_mesh(se, xy_cols = xy_col, cutoff = 1) #using lon/lat since this is on the reprojected regular lat/lon grid, and the domain crosses multiple UTM zones
    #MOM6 resolution is 1/12 = ~8 km

    ### get relative importance of model using type 3 anova method
    # model with *only* intercept and no random fields:
    fit_null <- sdmTMB::sdmTMB(
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
        fitSub <- sdmTMB::sdmTMB(
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


