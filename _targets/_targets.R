# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline

# Load packages required to define the pipeline:
library(targets)
library(crew)
library(tarchetypes)
# library(tarchetypes) # Load other packages as needed.

# Set target options:
tar_option_set(
  packages = c('ncdf4',
               'caret',
               'DescTools',
               'fields',
               'abind',
               'sf',
               'survdat',
               'dbutils',
               'measurements',
               'lubridate',
               'raster',
               'reshape2',
               'Matrix',
               'TMB',
               'sdmTMB',
               'sdmTMBextra',
               'future',
               'ranger',
               'sp',
               'akgfmaps',
               'EFHSDM',
               'terra',
               'meteo',
               'dismo',
               'gbm3',
               'gamm4',
               'ROCR',
               'tibble','sftime'),
   # Packages that your targets need for their tasks.
  error = 'null',

  
  # format = "qs", # Optionally set the default storage format. qs is fast.
  #
  # Pipelines that take a long time to run may benefit from
  # optional distributed computing. To use this capability
  # in tar_make(), supply a {crew} controller
  # as discussed at https://books.ropensci.org/targets/crew.html.
  # Choose a controller that suits your needs. For example, the following
  # sets a controller that scales up to a maximum of two workers
  # which run as local R processes. Each worker launches when there is work
  # to do and exits if 60 seconds pass with no tasks to run.
  #
  #   controller = crew::crew_controller_local(workers = 4)
  #
  # Alternatively, if you want workers to run on a high-performance computing
  # cluster, select a controller from the {crew.cluster} package.
  # For the cloud, see plugin packages like {crew.aws.batch}.
  # The following example is a controller for Sun Grid Engine (SGE).
  # 
  #   controller = crew.cluster::crew_controller_sge(
  #     # Number of workers that the pipeline can scale up to:
  #     workers = 10,
  #     # It is recommended to set an idle time so workers can shut themselves
  #     # down if they are not running tasks.
  #     seconds_idle = 120,
  #     # Many clusters install R as an environment module, and you can load it
  #     # with the script_lines argument. To select a specific verison of R,
  #     # you may need to include a version string, e.g. "module load R/4.3.2".
  #     # Check with your system administrator if you are unsure.
  #     script_lines = "module load R"
  #   )
  #
  # Set other options as needed.
)

# Run the R scripts in the R/ folder with your custom functions:
tar_source()
# tar_source("other_functions.R") # Source other scripts as needed.

#modify species list before target list because this is the only way it works for some reason
var.list <- create_var_list(varLong = c('Bottom Temperature'), 
                                        #'Bottom Oxygen', 
                                       # 'Sea Water Salinity at Sea Floor', 
                                        #'Bottom Aragonite Solubility'), 
                                        #'Sea Surface Temperature', 
                                        #'Sea Surface Salinity', 
                                        #'Surface pH', 
                                        #'Mixed layer depth (delta rho = 0.03)',
                                        #'Diazotroph new (NO3-based) prim. prod. integral in upper 100m', 
                                        #'Small phyto. new (NO3-based) prim. prod. integral in upper 100m', 
                                        #'Medium phyto. new (NO3-based) prim. prod. integral in upper 100m', 
                                        #'Large phyto. new (NO3-based) prim. prod. integral in upper 100m',
                                        #'Small zooplankton nitrogen biomass in upper 100m',
                                        #'Medium zooplankton nitrogen biomass in upper 100m',
                                        #'Large zooplankton nitrogen biomass in upper 100m',
                                        #'Water column net primary production vertical integral', 
                                        #'Downward Flux of Particulate Organic Carbon'), 
                            varShort = c('bottomT'), #'bottomO2', 'bottomS', 'bottomArg'), 
                                         #'surfaceT', 'surfaceS', 'surfacepH', 'MLD', 
                                         #'diazPP', 'smallPP', 'mediumPP', 'largePP', 
                                         #'smallZoo', 'mediumZoo', 'largeZoo', 'intNPP', 'POC'),
                            names = c("Long.Name", "Short.Name"))


# Replace the target list below with your own:
list(
  #get data 
  # tar_target(
  #   name = get_survey,
  #  command = standardize_data(dataType = 'Surveys', channel = dbutils::connect_to_database(server="NEFSC_pw_oraprod",uid="KGALLAGHER"))
  # ),
  # tar_target(
  #   name = get_observer,
  #  command = standardize_data(dataType = 'Observer', channel = dbutils::connect_to_database(server="NEFSC_pw_oraprod",uid="KGALLAGHER"))
  # ),
  tar_target(name = surv, command = 'nefsc_survey_data.csv', format = 'file'), #load in surv - already run standardize_data on NEFSC survey data
  tar_target(name = obs, command = 'observer_data.csv', format = 'file'), #load in observer data - already run standardize_data on observer dataset 
  tar_target(
    name = get_MENH,
    command = standardize_data(dataType = 'CSV', csv = "./csvs/MaineDMR_Trawl_Survey_Tow_Catch_2025-06-30.csv", csvCols = c('Start_Longitude', 'Start_Latitude', 'Start_Date', 'Number_Caught', 'Common_Name'))
  ),
  tar_target(
    name = get_MA,
    command = standardize_data(dataType = 'CSV', csv = "./csvs/MABottom_Trawl_2025-07-01.csv", csvCols = c('Lon', 'Lat', 'Date', 'Num', 'SCI_NAME'))
  ),
  tar_target(
    name = get_NJ,
    command = standardize_data(dataType = 'CSV', csv = "./csvs/NJOT_Tow_Catch_2025-07-01.csv", csvCols = c('START_LON', 'START_LAT', 'DATE.FORMAT', 'NUMBER', 'LATIN_NAME'))
  ),
  #normalize environmental data - now branched to save my sanity!! 
  #dynamic? branching
    getMOM6 <- tar_map(
      values = var.list,
      tar_target(
        name = pull_MOM6,
        command = pull_env(varURL = "https://psl.noaa.gov/cefi_portal/data_index/cefi_data_indexing.Projects.CEFI.regional_mom6.cefi_portal.northwest_atlantic.full_domain.hindcast.json", reqVars = Long.Name, shortNames = Short.Name)
      ),
      tar_target(
        name = avg_MOM6,
        command = avg_env(rastList = pull_MOM6)
    ), 
      tar_target(
        name = norm_MOM6,
        command = norm_env(rawList = pull_MOM6, avgList = avg_MOM6, shortNames = Short.Name)
      )
    ),
  tar_combine(
    name = combine_env, 
    getMOM6[['norm_MOM6']],
    command = list(!!!.x)
 ),
  #make rasters - will add static branching from here onwards in full runs
  tar_target(
    name = make_survey_rast,
    command = create_rast(data = surv, dataType = 'Surveys', grid = "http://psl.noaa.gov/thredds/dodsC/Projects/CEFI/regional_mom6/cefi_portal/northwest_atlantic/full_domain/hindcast/monthly/regrid/r20230520/tob.nwa.full.hcast.monthly.regrid.r20230520.199301-201912.nc", targetVec = c('Summer Flounder', 'Fluke', "PARALICHTHYS DENTATUS"))
  ), 
  tar_target(
    name = make_observer_rast,
    command = create_rast(data = obs, dataType = 'Observer', grid = "http://psl.noaa.gov/thredds/dodsC/Projects/CEFI/regional_mom6/cefi_portal/northwest_atlantic/full_domain/hindcast/monthly/regrid/r20230520/tob.nwa.full.hcast.monthly.regrid.r20230520.199301-201912.nc", targetVec = c('Summer Flounder', 'Fluke', "PARALICHTHYS DENTATUS"))
  ),
  tar_target(
    name = make_MENH_rast,
    command = create_rast(get_MENH, dataType = 'csv', grid = "http://psl.noaa.gov/thredds/dodsC/Projects/CEFI/regional_mom6/cefi_portal/northwest_atlantic/full_domain/hindcast/monthly/regrid/r20230520/tob.nwa.full.hcast.monthly.regrid.r20230520.199301-201912.nc", targetVec = c('Summer Flounder', 'Fluke', "PARALICHTHYS DENTATUS"))
  ),
  tar_target(
    name = make_MA_rast,
    command = create_rast(get_MA, dataType = 'csv', grid = "http://psl.noaa.gov/thredds/dodsC/Projects/CEFI/regional_mom6/cefi_portal/northwest_atlantic/full_domain/hindcast/monthly/regrid/r20230520/tob.nwa.full.hcast.monthly.regrid.r20230520.199301-201912.nc", targetVec = c('Summer Flounder', 'Fluke', "PARALICHTHYS DENTATUS"))
  ),
  tar_target(
    name = make_NJ_rast,
    command = create_rast(get_NJ, dataType = 'csv', grid = "http://psl.noaa.gov/thredds/dodsC/Projects/CEFI/regional_mom6/cefi_portal/northwest_atlantic/full_domain/hindcast/monthly/regrid/r20230520/tob.nwa.full.hcast.monthly.regrid.r20230520.199301-201912.nc", targetVec = c('Summer Flounder', 'Fluke', "PARALICHTHYS DENTATUS"))
  ),
  #merge rasters
  #tar_combine(
   # name = raster_list, 
  #  make_survey_rast, 
   # make_observer_rast, 
   # make_MENH_rast, 
   # make_MA_rast, 
   # make_NJ_rast,
   # command = list(!!!.x)
 # ),
  tar_target(
    name = merge_rasters,
    command = merge_rasts(list(make_survey_rast, 
                                make_observer_rast, 
                                make_MENH_rast, 
                                make_MA_rast, 
                                make_NJ_rast))
  ),
  #merge environmental and species data
  tar_target(
    name = merge_env_spp,
    command = merge_spp_env(rastStack = merge_rasters, envData = combine_env, addStatic = TRUE, staticData = 'staticVariables.RData')
  )
  #match guild
 # tar_target(
  #  name = match_guild,
  #  command = match_guilds(spp_env = merge_env_spp, spp = c('Summer Flounder', 'Fluke', "PARALICHTHYS DENTATUS"), spp_col = 'name', spp_guild = 'spp_list.csv', feeding_key = 'feeding_guilds.csv', feeding_col = 'Feeding.Guild', habitat_key = 'habitat_guilds.csv',  habitat_col = 'Habitat.Guild', static_vars = c('x', 'y', 'month_num', 'year', 'bathy', 'rugosity', 'coast_dist'), pa_col = 'value')
 # ), 
 #remove correlated variables 
# tar_target(
 #  name = rm_corr, 
 #  command = remove_corr(match_guild, pa_col = 'value', xy_col = c('x', 'y'), month_col = 'month', year_col = 'year')
# )
  #make models
  #tar_target(
   # name = gam.model,
   # command = make_sdm(se = rm_corr, pa_col = 'value', xy_col = c('x', 'y'), month_col = 'month', year_col = 'year', model = 'gam')
 # ),
 # tar_target(
  #  name = maxent.model,
  #  command = make_sdm(se = rm_corr, pa_col = 'value', xy_col = c('x', 'y'), month_col = 'month', year_col = 'year', model = 'maxent')
 # ),
 # tar_target(
  #  name = rf.model,
   # command = make_sdm(se = rm_corr, pa_col = 'value', xy_col = c('x', 'y'), month_col = 'month', year_col = 'year', model = 'rf')
 # ),
 # tar_target(
  #  name = brt.model,
  #  command = make_sdm(se = rm_corr, pa_col = 'value', xy_col = c('x', 'y'), month_col = 'month', year_col = 'year', model = 'brt')
 # ),
 # tar_target(
  #  name = tmb.model,
   # command = make_sdm(se = rm_corr, pa_col = 'value', xy_col = c('x', 'y'), month_col = 'month', year_col = 'year', model = 'sdmtmb')
 # ),
 #cross-validate models
 #tar_target(
 # name = gam.cv,
 # command = sdm_cv(mod = gam.model, se = rm_corr, pa_col = 'value', xy_col = c('x', 'y'), month_col = 'month', year_col = 'year', model = 'gam')
 # ),
 # tar_target(
 #  name = maxent.cv,
 #  command = sdm_cv(mod = maxent.model, se = rm_corr, pa_col = 'value', xy_col = c('x', 'y'), month_col = 'month', year_col = 'year', model = 'maxent')
 # ),
 # tar_target(
 #  name = rf.cv,
 # command = sdm_cv(mod = rf.model, se = rm_corr, pa_col = 'value', xy_col = c('x', 'y'), month_col = 'month', year_col = 'year', model = 'rf')
 # ),
 # tar_target(
 #  name = brt.cv,
 #  command = sdm_cv(mod = brt.model, se = rm_corr, pa_col = 'value', xy_col = c('x', 'y'), month_col = 'month', year_col = 'year', model = 'brt')
 # ),
 # tar_target(
 #  name = tmb.cv,
 # command = sdm_cv(mod = tmb.model, se = rm_corr, pa_col = 'value', xy_col = c('x', 'y'), month_col = 'month', year_col = 'year', model = 'sdmtmb')
 # ),
 #get predictions from models
 #tar_target(
 # name = gam.preds,
 # command = sdm_preds(cv = gam.cv, model = 'gam')
 # ),
 # tar_target(
 #  name = maxent.preds,
 #  command = sdm_preds(cv = maxent.cv, model = 'maxent')
 # ),
 # tar_target(
 #  name = rf.preds,
 # command = sdm_preds(cv = rf.cv, model = 'rf')
 # ),
 # tar_target(
 #  name = brt.preds,
 #  command = sdm_preds(cv = brt.cv, model = 'brt')
 # ),
 # tar_target(
 #  name = tmb.preds,
 # command = sdm_preds(cv = tmb.cv, model = 'sdmtmb')
 # ),
 #get evaluation metric from models
 #tar_target(
 # name = gam.eval,
 # command = sdm_eval(preds = gam.preds, metric = 'auc', model = 'gam')
 # ),
 # tar_target(
 #  name = maxent.eval,
 #  command = sdm_eval(preds = maxent.preds, metric = 'auc', model = 'maxent')
 # ),
 # tar_target(
 #  name = rf.eval,
 # command = sdm_eval(preds = rf.preds, metric = 'auc', model = 'rf')
 # ),
 # tar_target(
 #  name = brt.eval,
 #  command = sdm_eval(preds = brt.preds, metric = 'auc', model = 'brt')
 # ),
 # tar_target(
 #  name = tmb.eval,
 # command = sdm_eval(preds = tmb.preds, metric = 'auc', model = 'sdmtmb')
 # ),
 #get deviance explained/importance from models
 #tar_target(
 # name = gam.importance,
 # command = sdm_importance(mod = gam.model, se = rm_corr, pa_col = 'value', xy_col = c('x', 'y'), month_col = 'month', year_col = 'year', model = 'gam')
 # ),
 # tar_target(
 #  name = maxent.importance,
 #  command = sdm_importance(mod = maxent.model, se = rm_corr, pa_col = 'value', xy_col = c('x', 'y'), month_col = 'month', year_col = 'year', model = 'maxent')
 # ),
 # tar_target(
 #  name = rf.importance,
 # command = sdm_importance(mod = rf.model, se = rm_corr, pa_col = 'value', xy_col = c('x', 'y'), month_col = 'month', year_col = 'year', model = 'rf')
 # ),
 # tar_target(
 #  name = brt.importance,
 #  command = sdm_importance(mod = brt.model, se = rm_corr, pa_col = 'value', xy_col = c('x', 'y'), month_col = 'month', year_col = 'year', model = 'brt')
 # ),
 # tar_target(
 #  name = tmb.importance,
 # command = sdm_importance(mod = tmb.model, se = rm_corr, pa_col = 'value', xy_col = c('x', 'y'), month_col = 'month', year_col = 'year', model = 'sdmtmb')
 # ),
  #predict models
  #tar_target(
  #  name = predict.gam,
  #  command = make_predictions(mod = gam.model, model = 'gam', rasts = combine_env, mask = T, bathy_nm = 'bathy', bathy_max = 1000, se = rm_corr)
 # ),
 # tar_target(
  #  name = predict.maxent,
  #  command = make_predictions(mod = maxent.model, model = 'maxent', rasts = combine_env, mask = T, bathy_nm = 'bathy', bathy_max = 1000, se = rm_corr)
# ),
 # tar_target(
  #  name = predict.rf,
  #  command = make_predictions(mod = rf.model, model = 'rf', rasts = combine_env, mask = T, bathy_nm = 'bathy', bathy_max = 1000, se = rm_corr)
 # ),
 # tar_target(
  #  name = predict.brt,
   # command = make_predictions(mod = brt.model, model = 'brt', rasts = combine_env, mask = T, bathy_nm = 'bathy', bathy_max = 1000, se = rm_corr)
 # ),
 # tar_target(
  #  name = predict.tmb,
  #  command = make_predictions(mod = tmb.model, model = 'sdmtmb', rasts = combine_env, mask = T, bathy_nm = 'bathy', bathy_max = 1000, se = rm_corr)
 # ),
  #build ensemble
#tar_target(
# name = ensemble.weights,
#command = MakeEnsemble(rmse = c(gam.eval, maxent.eval, rf.eval, brt.eval, tmb.eval)) #make weights
#)
 # tar_target(
  #  name = build.ensemble,
  #  command = make_sdm(model = 'ens', ensembleWeights = ensemble.weights, ensemblePreds = list(gam.preds, maxent.preds, brt.preds, rf.preds, tmb.preds))
 # ),
  #predict ensemble
 # tar_target(
  #  name = predict.ensemble,
  #  command = make_predictions(model = 'ens', rasts = list(predict.gam, predict.maxent, predict.brt, predict.rf, predict.tmb), weights = ensemble.weights)
 # )
)

