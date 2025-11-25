###building out targets workflow independently of targets
#uses functions built for targets and saves everything in species specific folders

##################################
#####SET UP - LOAD EVERY TIME ####
##################################

###load libraries
library(ncdf4)
library(caret)
library(DescTools)
library(fields)
library(abind)
library(sf)
library(sftime)
library(survdat)
library(dbutils)
library(measurements)
library(lubridate)
library(raster)
library(reshape2)
library(Matrix)
library(TMB)
library(sdmTMB)
library(sdmTMBextra)
library(future)
library(ranger)
library(sp)
library(akgfmaps)
library(EFHSDM)
library(terra)
library(meteo)
library(dismo)
library(gbm3)
library(gamm4)
library(ROCR)
library(sftime)
library(furrr)

### set working directory
setwd('/home/kgallagher/ClimateVulnerabilityAssessment2.0/SDMs')


### source functions
targets::tar_source("/home/kgallagher/ClimateVulnerabilityAssessment2.0/functions")

#create species folders and appropriate subfolders 
#spp.list <- create_spp_list('spp_list.csv')

spp.list <- read.csv('spp_list.csv')
#spp.list <- spp.list[,c(1:6)]
spp.list$Name <- gsub(' ', '', spp.list$Common.Name)


#make directory for each species if it doesn't exist; if directory exists, it is not changed
for(x in 1:nrow(spp.list)){
  dir.create(file.path(getwd(),spp.list$Name[x]), showWarnings = T) #main folder 
  dir.create(file.path(getwd(),spp.list$Name[x], 'input_rasters'), showWarnings = T) #input raster folder
  dir.create(file.path(getwd(),spp.list$Name[x], 'output_rasters'), showWarnings = T) #output data folder
  dir.create(file.path(getwd(),spp.list$Name[x], 'model_output'), showWarnings = T) #model_output folder
  dir.create(file.path(getwd(),spp.list$Name[x], 'model_output', 'models'), showWarnings = T) #model_output/models folder
  dir.create(file.path(getwd(),spp.list$Name[x], 'model_output', 'cvs'), showWarnings = T) #model_output/cvs folder
  dir.create(file.path(getwd(),spp.list$Name[x], 'model_output', 'preds'), showWarnings = T) #model_output/preds folder
  dir.create(file.path(getwd(),spp.list$Name[x], 'model_output', 'eval_metrics'), showWarnings = T) #model_output/eval_metrics folder
  dir.create(file.path(getwd(),spp.list$Name[x], 'model_output', 'importance'), showWarnings = T) #model_output/importance folder
}

##################################

##############################
#####GET FISHERIES DATA ######
##############################

#survey data - needs VPN 
surv <- standardize_data(dataType = 'Surveys', channel = dbutils::connect_to_database(server="NEFSC_pw_oraprod",uid="KGALLAGHER"))
write.csv(surv, './Data/csvs/standardized/survey.csv')

#observer data - needs VPN
obs <- standardize_data(dataType = 'Observer', channel = dbutils::connect_to_database(server="NEFSC_pw_oraprod",uid="KGALLAGHER"))
write.csv(obs, './Data/csvs/standardized/observer.csv')

## State run surveys 
#maine/new hampshire
menh <-  standardize_data(dataType = 'CSV', csv = "./Data/csvs/raw/MaineDMR_Trawl_Survey_Tow_Catch_2025-07-17.csv", csvCols = c('towID', 'Start_Longitude', 'Start_Latitude', 'Start_Date', 'Number_Caught', 'Common_Name'))
write.csv(menh, './Data/csvs/standardized/MENH.csv')

#mass
mass <- standardize_data(dataType = 'CSV', csv = "./Data/csvs/raw/MABottom_Trawl_2025-08-6.csv", csvCols = c('towID', 'Lon', 'Lat', 'Date', 'Num', 'SCI_NAME'))
write.csv(mass, './Data/csvs/standardized/MA.csv')

#new jersey
nj <- standardize_data(dataType = 'CSV', csv = "./Data/csvs/raw/NJOT_Tow_Catch_2025-07-01.csv", csvCols = c('TOW_ID', 'START_LON', 'START_LAT', 'DATE.FORMAT', 'NUMBER', 'LATIN_NAME'))
write.csv(nj, './Data/csvs/standardized/NJ.csv')

#ct
ct <- standardize_data(dataType = 'CSV', csv = "./Data/csvs/raw/CT_Tow_Catch_2025-08-06.csv", csvCols = c('Sample.Number', 'Longitude', 'Latitude', 'Date', 'TotalCount', 'SCI_NAME'))
write.csv(ct, './Data/csvs/standardized/CT.csv')

#delaware
de <- standardize_data(dataType = 'CSV', csv = "./Data/csvs/raw/DE_Tow_Catch_2025-07-18.csv", csvCols = c('towID', 'LONDD', 'LATDD', 'date', 'number', 'SCI_NAME'))
write.csv(de, './Data/csvs/standardized/DE.csv')

#neamap 
neamap <- standardize_data(dataType = 'CSV', csv = "./Data/csvs/raw/NEAMAP_Tow_Catch_2025-08-14.csv", csvCols = c('station', 'lon', 'lat', 'date', 'present_absent', 'SCI_NAME'))
write.csv(neamap, './Data/csvs/standardized/NEAMAP.csv')


##############################

##############################
#####GET MOM6 DATA ###########
##############################
var.list <- create_var_list(varLong = c('Bottom Temperature', 
                             'Bottom Oxygen', 
                             'Sea Water Salinity at Sea Floor', 
                            'Bottom Aragonite Solubility', 
                             'Sea Surface Temperature', 
                             'Sea Surface Salinity', 
                             'Surface pH', 
                             'Mixed layer depth (delta rho = 0.03)',
                             'Diazotroph new (NO3-based) prim. prod. integral in upper 100m', 
                             'Small phyto. new (NO3-based) prim. prod. integral in upper 100m', 
                             'Medium phyto. new (NO3-based) prim. prod. integral in upper 100m', 
                            'Large phyto. new (NO3-based) prim. prod. integral in upper 100m',
                             'Small zooplankton nitrogen biomass in upper 100m',
                             'Medium zooplankton nitrogen biomass in upper 100m',
                             'Large zooplankton nitrogen biomass in upper 100m',
                             'Water column net primary production vertical integral', 
                             'Downward Flux of Particulate Organic Carbon'), 
                            varShort = c('bottomT', 'bottomO2', 'bottomS', 'bottomArg', 
                             'surfaceT', 'surfaceS', 'surfacepH', 'MLD', 
                             'diazPP', 'smallPP', 'mediumPP', 'largePP', 
                             'smallZoo', 'mediumZoo', 'largeZoo', 'intNPP', 'POC'),
                            names = c("Long.Name", "Short.Name"))

raw <- avg <- norm <- vector(mode = 'list', length = nrow(var.list))

cluster <- makeCluster(10, type='PSOCK')
registerDoParallel(cluster)
raw <- foreach(x = 1:nrow(var.list), .packages = c("ncdf4", 'raster', 'jsonlite')) %dopar% {
#for(x in 1:nrow(var.list)){
  r <- pull_env(varURL = "https://psl.noaa.gov/cefi_portal/data_index/cefi_data_indexing.Projects.CEFI.regional_mom6.cefi_portal.northwest_atlantic.full_domain.hindcast.json", reqVars = var.list$Long.Name[x], shortNames = var.list$Short.Name[x], release = 'r20230520')
  
  #a <- avg_env(r)
  
  #n <- norm_env(rawList = r, avgList = a, shortNames = var.list$Short.Name[x])
  
  raw[[x]] <- r
  #avg[[x]] <- a
  #norm[[x]] <- n
  
  #print(x)
}
stopCluster(cluster)
names(raw) <- var.list$Short.Name
save(raw, file = '/home/kgallagher/ClimateVulnerabilityAssessment2.0/SDMs/Data/MOM6/raw_MOM6_082025.RData')

cluster <- makeCluster(10, type='PSOCK')
registerDoParallel(cluster)
avg <- foreach(x = 1:nrow(var.list), .packages = c("ncdf4", 'raster', 'jsonlite')) %do% {
  #for(x in 1:nrow(var.list)){
 # r <- pull_env(varURL = "https://psl.noaa.gov/cefi_portal/data_index/cefi_data_indexing.Projects.CEFI.regional_mom6.cefi_portal.northwest_atlantic.full_domain.hindcast.json", reqVars = var.list$Long.Name[x], shortNames = var.list$Short.Name[x], release = 'r20230520')
  
  a <- avg_env(raw[[x]])
  
  #n <- norm_env(rawList = r, avgList = a, shortNames = var.list$Short.Name[x])
  
 # raw[[x]] <- r
  avg[[x]] <- a
  #norm[[x]] <- n
  
  #print(x)
}
stopCluster(cluster)
names(avg) <- var.list$Short.Name
save(avg, file = '/home/kgallagher/ClimateVulnerabilityAssessment2.0/SDMs/Data/MOM6/avg_MOM6_082025.RData')

cluster <- makeCluster(10, type='PSOCK')
registerDoParallel(cluster)
norm <- foreach(x = 1:nrow(var.list), .packages = c("ncdf4", 'raster', 'jsonlite', 'abind')) %dopar% {
  #for(x in 1:nrow(var.list)){
  #r <- pull_env(varURL = "https://psl.noaa.gov/cefi_portal/data_index/cefi_data_indexing.Projects.CEFI.regional_mom6.cefi_portal.northwest_atlantic.full_domain.hindcast.json", reqVars = var.list$Long.Name[x], shortNames = var.list$Short.Name[x], release = 'r20230520')
  
  #a <- avg_env(r)
  
  n <- norm_env(rawList = raw[[x]], avgList = avg[[x]], shortNames = var.list$Short.Name[x])
  
 # raw[[x]] <- r
  #avg[[x]] <- a
  norm[[x]] <- n
  
  #print(x)
}
stopCluster(cluster)
names(norm) <- var.list$Short.Name
save(norm, file = '/home/kgallagher/ClimateVulnerabilityAssessment2.0/SDMs/Data/MOM6/norm_MOM6_082025.RData')

##############################

##############################
##### BUILD FISHERIES RASTERS#
##############################

#build rasters seperately for each source - takes about 24 hours in parallel
sources <- c('survey',  
             'MENH', 
             'MA', 
             'NJ', 
             'CT', 
             "DE", 
             'NEAMAP',
             'observer')

#skip = F #overwrite existing files
saveRast <- function(csvName, spp, sppNames, skip){

  #require(logr)
  
  data <- read.csv(paste0('./Data/csvs/standardized/', csvName, '.csv'))
  
  #open log file
  sink(file = file.path(getwd(), 'logs', paste0(csvName, '.log')), append = T)
  #sink(file = file.path(getwd(), 'logs', paste0(csvName, '.log')), append = T, type = 'message')
  
  # Ensure the sinks are closed when the function exits, regardless of how it exits.
  on.exit({
    #sink(type = "message")
    sink()
  })
  
  print(paste(csvName , spp, sep = '-'))
  
 # for(s in 1:length(spp.list$Common.Name)){
    print(Sys.time())
    #print(spp)
    
      if(skip){
        if(file.exists(paste0(file.path(getwd(), spp, 'input_rasters'), '/', csvName, '.nc'))){
            print('file exists and skip == T, so skipping this file!')
            return(NA)
        } else {
    
          nms <- strsplit(sppNames, split = ',')[[1]]
          
          if(csvName != 'observer'){
            rast <- create_rast(data = data, dataType = 'Surveys', grid = "http://psl.noaa.gov/thredds/dodsC/Projects/CEFI/regional_mom6/cefi_portal/northwest_atlantic/full_domain/hindcast/monthly/regrid/r20230520/tob.nwa.full.hcast.monthly.regrid.r20230520.199301-201912.nc", targetVec = nms)
            print(range(rast[]))
            raster::writeRaster(rast, filename = paste0(file.path(getwd(),spp, 'input_rasters'), '/', csvName, '.nc'), bylayer = F,overwrite = T)
          } else {
            rast <- create_rast(data = data, dataType = 'Observer', grid = "http://psl.noaa.gov/thredds/dodsC/Projects/CEFI/regional_mom6/cefi_portal/northwest_atlantic/full_domain/hindcast/monthly/regrid/r20230520/tob.nwa.full.hcast.monthly.regrid.r20230520.199301-201912.nc", targetVec = nms)
            if(is.null(rast)){
              print('rast is NULL - minimum conditions not met')
            } else {
              print(range(rast[]))
              raster::writeRaster(rast, filename = paste0(file.path(getwd(),spp, 'input_rasters'), '/', csvName, '.nc'), bylayer = F,overwrite = T)
            }
          }
      
          return(range(rast[]))
    
      }#end skip && if file is present
      } else { #if skip = F, just run it without checking
        nms <- strsplit(sppNames, split = ',')[[1]]
        
        if(csvName != 'observer'){
          rast <- create_rast(data = data, dataType = 'Surveys', grid = "http://psl.noaa.gov/thredds/dodsC/Projects/CEFI/regional_mom6/cefi_portal/northwest_atlantic/full_domain/hindcast/monthly/regrid/r20230520/tob.nwa.full.hcast.monthly.regrid.r20230520.199301-201912.nc", targetVec = nms)
          print(range(rast[]))
          raster::writeRaster(rast, filename = paste0(file.path(getwd(),spp, 'input_rasters'), '/', csvName, '.nc'), bylayer = F,overwrite = T)
        } else {
          rast <- create_rast(data = data, dataType = 'Observer', grid = "http://psl.noaa.gov/thredds/dodsC/Projects/CEFI/regional_mom6/cefi_portal/northwest_atlantic/full_domain/hindcast/monthly/regrid/r20230520/tob.nwa.full.hcast.monthly.regrid.r20230520.199301-201912.nc", targetVec = nms)
          if(is.null(rast)){
            print('rast is NULL - minimum conditions not met')
          } else {
            print(range(rast[]))
            raster::writeRaster(rast, filename = paste0(file.path(getwd(),spp, 'input_rasters'), '/', csvName, '.nc'), bylayer = F,overwrite = T)
          }
        }
        
        return(range(rast[]))
  }
  
  sink() #close log file
} #end function

#for(y in 1:length(sources)){
#saveRast(csvName = sources[y], spp.list = spp.list, skip = T)
#}

args <- tidyr::expand_grid(csvName = sources, spp = spp.list$Name, skip = T)

args2 <- merge(x = args, y = spp.list[,c(1:7,11)], by.x = 'spp', by.y = 'Name')

altNames <- paste(args2$Common.Name, args2$COM_NAME, args2$Scientific.Name, args2$Alternate.Name, args2$SCI_NAME, args2$SCI_NAME_ALT, args2$SCI_NAME_ALT2, sep = ',')   

plan(multisession, workers = 2)
#sink(file = 'rasters.log', append = T)
checks <- future_pmap(list(..1 = args2$csvName, ..2 = args2$spp, ..3 = altNames, ..4 = args2$skip), ~ saveRast(csvName = ..1, spp = ..2, sppNames = ..3, skip = ..4), .progress = T)
#sink()
plan(sequential)

#cl <- makeCluster(3)
#clusterExport(cl, c('saveRast', 'create_rast', 'spp.list', 'sources'))
#clusterApplyLB(cl, sources, saveRast, spp.list = spp.list, skip = T)
#stopCluster(cl)

### put them all together - takes about an hour and a half
args <- tidyr::expand_grid(name = spp.list$Name, skip = T)

combineSave <- function(name, skip){
  sink(file.path(getwd(), 'logs', 'combineRasters.log'), append = T)
  # Ensure the sinks are closed when the function exits, regardless of how it exits.
  on.exit({
    #sink(type = "message")
    sink()
  })
  print(Sys.time())
  print(name)
  
  if(skip){
    if(file.exists(paste(file.path(getwd(),name, 'input_rasters'), 'combined_rasters.nc', sep = '/'))){
      print('file exists and skip == T, so skipping this file!')
      return(NA)
    } else {
      flist <-  dir(file.path(getwd(),name, 'input_rasters'), full.names = T) 
      rasts <- vector(mode = 'list', length = length(flist))
      for(r in 1:length(flist)){
        rast <- brick(flist[r]) #rast
        rasts[[r]] <- rast
        #print(r)
      }
      
      combinedRasts <- merge_rasts(rasts)
      print(range(combinedRasts[]))
      writeRaster(combinedRasts, filename = paste(file.path(getwd(),name, 'input_rasters'), 'combined_rasters.nc', sep = '/'), bylayer = F, overwrite = T)
    } #e3nd else
  } else { #if skip = F, do it anyway
    
    flist <-  dir(file.path(getwd(),name, 'input_rasters'), full.names = T) 
    rasts <- vector(mode = 'list', length = length(flist))
    for(r in 1:length(flist)){
      rast <- brick(flist[r]) #rast
      rasts[[r]] <- rast
      #print(r)
    }
    
    combinedRasts <- merge_rasts(rasts)
    print(range(combinedRasts[]))
    writeRaster(combinedRasts, filename = paste(file.path(getwd(),name, 'input_rasters'), 'combined_rasters.nc', sep = '/'), bylayer = F, overwrite = T)
    #sink()
  } #end else 
}
#skip = F #overwrite existing files

#for(x in 1:length(spp.list$Name)){
plan(multisession, workers = 2)
combs <- future_pmap(list(..1 = args$name, ..2 = args$skip), ~ combineSave(name = ..1, skip = ..2), .progress = T)
plan(sequential)
#sink()

###CHECKS BEFORE MOVING ON  
flist <- dir(path = getwd(), pattern = 'combined_rasters.nc', recursive = T, full.names = T)
for(x in 24:length(flist)){
  r <- brick(flist[x])
  print(range(r[]))
}

##############################

##############################
##### MAKE DATA FRAMES #######
##############################

load('./Data/MOM6/norm_MOM6_082025.RData') #norm

#for(x in 1:nrow(spp.list)){
makeDF <- function(name, skip, mMin, mMax, yMin, yMax){
  sink(file.path(getwd(), 'logs', 'build_data_frames.log'), append = T)
  # Ensure the sinks are closed when the function exits, regardless of how it exits.
  on.exit({
    #sink(type = "message")
    sink()
  })
  
  print(Sys.time())
  print(name)
  
  if(skip){
    if(file.exists(paste(file.path(getwd(),name), 'pa_clean.RData', sep = '/'))){
      print('file exists and skip == T, so skipping this file!')
      return(NA)
    } else {
      combinedRasts <- brick(paste(file.path(getwd(),name, 'input_rasters'), 'combined_rasters.nc', sep = '/')) #combinedRasts
      
      my <- expand.grid(mMin:mMax, yMin:yMax)
      names(combinedRasts) <- paste(sprintf("%02d", my$Var1), my$Var2, sep = '.')
      
      print(paste0(name, '- merging spp and env data -', Sys.time()))
      df <- merge_spp_env(rastStack = combinedRasts, envData = norm, addStatic = TRUE, staticData = './Data/staticVariables_cropped.RData')
      save(df, file = paste(file.path(getwd(),name), 'presence_absence_environment.RData', sep = '/'))
      
      print(paste0(name, '- matching guild -', Sys.time()))
      
      dfG <- match_guilds(spp_env = df, spp = c('Summer Flounder', 'Fluke', "PARALICHTHYS DENTATUS"), spp_col = 'SCI_NAME', spp_guild = 'spp_list.csv', feeding_key = 'feeding_guilds.csv', feeding_col = 'Feeding.Guild', habitat_key = 'habitat_guilds.csv',  habitat_col = 'Habitat.Guild', static_vars = c('x', 'y', 'month', 'year', 'bathy', 'rugosity', 'dist2coast'), pa_col = 'value')
      save(dfG, file = paste(file.path(getwd(),name), 'pa_guild.RData', sep = '/'))
      
      print(paste0(name, '- converting to binary -', Sys.time()))
      dfB <- clean_data(dfG, pa_col = 'value')
      save(dfB, file = paste(file.path(getwd(),name), 'pa_binary.RData', sep = '/'))
      
      print(paste0(name, ' - removing correlated variables -', Sys.time()))
      dfC <- remove_corr(dfB, pa_col = 'value', xy_col = c('x', 'y'), month_col = 'month', year_col = 'year')
      save(dfC, file = paste(file.path(getwd(),name), 'pa_clean.RData', sep = '/'))
    }
  } else {
    combinedRasts <- brick(paste(file.path(getwd(),name, 'input_rasters'), 'combined_rasters.nc', sep = '/')) #combinedRasts
    
    my <- expand.grid(mMin:mMax, yMin:yMax)
    names(combinedRasts) <- paste(sprintf("%02d", my$Var1), my$Var2, sep = '.')
    
    print(paste0(name, '- merging spp and env data -', Sys.time()))
    df <- merge_spp_env(rastStack = combinedRasts, envData = norm, addStatic = TRUE, staticData = './Data/staticVariables_cropped.RData')
    save(df, file = paste(file.path(getwd(),name), 'presence_absence_environment.RData', sep = '/'))
    
    print(paste0(name, '- matching guild -', Sys.time()))
    
    dfG <- match_guilds(spp_env = df, spp = c('Summer Flounder', 'Fluke', "PARALICHTHYS DENTATUS"), spp_col = 'SCI_NAME', spp_guild = 'spp_list.csv', feeding_key = 'feeding_guilds.csv', feeding_col = 'Feeding.Guild', habitat_key = 'habitat_guilds.csv',  habitat_col = 'Habitat.Guild', static_vars = c('x', 'y', 'month', 'year', 'bathy', 'rugosity', 'dist2coast'), pa_col = 'value')
    save(dfG, file = paste(file.path(getwd(),name), 'pa_guild.RData', sep = '/'))
    
    print(paste0(name, '- converting to binary -', Sys.time()))
    dfB <- clean_data(dfG, pa_col = 'value')
    save(dfB, file = paste(file.path(getwd(),name), 'pa_binary.RData', sep = '/'))
    
    print(paste0(name, ' - removing correlated variables -', Sys.time()))
    dfC <- remove_corr(dfB, pa_col = 'value', xy_col = c('x', 'y'), month_col = 'month', year_col = 'year')
    save(dfC, file = paste(file.path(getwd(),name), 'pa_clean.RData', sep = '/'))
  }
}

#args <- tidyr::expand_grid(name = spp.list$Name, skip = T) #create list of arguments for loop

#plan(multisession, workers = 2) ###cannot do in parallel because norm data is too big to share across nodes, but this should only take about 2.5 hours 
for(x in 1:nrow(spp.list)){
  makeDF(name = spp.list$Name[x], skip = T, mMin = 1, mMax = 12, yMin = 1993, yMax = 2019)
}
#plan(sequential)

##############################

##############################
##### MAKE MODELS  ###########
##############################

models <- c(#'gam', 
           # 'maxent', 
           # 'rf', 
            'brt') 
          #'sdmtmb')

#load('norm_MOM6_082025.RData') #norm

args <- tidyr::expand_grid(model = models, spp = spp.list$Name, skip = T)

#for(x in 1:nrow(spp.list)){
makeMods <- function(spp, model, skip){
  #open log file
  sink(file = file.path(getwd(), 'logs', paste0(model, '.log')), append = T)
  #sink(file = file.path(getwd(), 'logs', paste0(csvName, '.log')), append = T, type = 'message')
  
  # Ensure the sinks are closed when the function exits, regardless of how it exits.
  on.exit({
    #sink(type = "message")
    sink()
  })
  
  print(Sys.time())
  print(spp)
  
  load(paste(file.path(getwd(),spp), 'pa_clean.RData', sep = '/')) #load data - dfC
  
  #for(y in models){ 
    #print(y)
    
  if(skip){
    if(file.exists(paste0(file.path(getwd(),spp, 'model_output', 'models'), '/', toupper(model), '.RData'))){
      print('model exists and skip == T, so skipping this file!')
      #return(NA)
    } else {
      print(paste0(spp, '- making model - ', Sys.time()))
      mod <- make_sdm(se = dfC, pa_col = 'value', xy_col = c('x', 'y'), month_col = 'month', year_col = 'year', model = model)
      save(mod, file = paste0(file.path(getwd(),spp, 'model_output', 'models'), '/', toupper(model), '.RData'))
    }
  } else {
    print(paste0(spp, '- making model - ', Sys.time()))
    mod <- make_sdm(se = dfC, pa_col = 'value', xy_col = c('x', 'y'), month_col = 'month', year_col = 'year', model = model)
    save(mod, file = paste0(file.path(getwd(),spp, 'model_output', 'models'), '/', toupper(model), '.RData'))
  }
    
  if(skip){
    if(file.exists(paste0(file.path(getwd(),spp, 'model_output', 'cvs'), '/', toupper(model), '.RData'))){
      print('cv exist and skip == T, so skipping this file!')
      #return(NA)
    } else {
      print(paste0(spp, '- performing CV - ', Sys.time()))
      load(paste0(file.path(getwd(),spp, 'model_output', 'models'), '/', toupper(model), '.RData'))
      cv <- sdm_cv(mod = mod, se = dfC, pa_col = 'value', xy_col = c('x', 'y'), month_col = 'month', year_col = 'year', model = model)
      save(cv, file = paste0(file.path(getwd(),spp, 'model_output', 'cvs'), '/', toupper(model), '.RData'))
    }
  } else {
    print(paste0(spp, '- performing CV - ', Sys.time()))
    load(paste0(file.path(getwd(),spp, 'model_output', 'models'), '/', toupper(model), '.RData'))
    cv <- sdm_cv(mod = mod, se = dfC, pa_col = 'value', xy_col = c('x', 'y'), month_col = 'month', year_col = 'year', model = model)
    save(cv, file = paste0(file.path(getwd(),spp, 'model_output', 'cvs'), '/', toupper(model), '.RData'))
  }

  if(skip){
    if(file.exists(paste0(file.path(getwd(),spp, 'model_output', 'preds'), '/',toupper(model), '.RData'))){
      print('preds exist and skip == T, so skipping this file!')
      #return(NA)
    } else {
      print(paste0(spp, '- Getting Preds - ', Sys.time()))
      load(paste0(file.path(getwd(),spp, 'model_output', 'cvs'), '/', toupper(model), '.RData'))
      preds <- sdm_preds(cv = cv, model = model)
      save(preds, file = paste0(file.path(getwd(),spp, 'model_output', 'preds'), '/',toupper(model), '.RData'))
    }
  } else {
    print(paste0(spp, '- Getting Preds - ', Sys.time()))
    load(paste0(file.path(getwd(),spp, 'model_output', 'cvs'), '/', toupper(model), '.RData'))
    preds <- sdm_preds(cv = cv, model = model)
    save(preds, file = paste0(file.path(getwd(),spp, 'model_output', 'preds'), '/',toupper(model), '.RData'))
  }
    
  if(skip){
    if(file.exists(paste0(file.path(getwd(),spp, 'model_output', 'eval_metrics'), '/',toupper(model), '.RData'))){
      print('eval metrics exist and skip == T, so skipping this file!')
      #return(NA)
    } else {
      print(paste0(spp, '- Evaluating Model - ', Sys.time()))
      load(paste0(file.path(getwd(),spp, 'model_output', 'preds'), '/',toupper(model), '.RData'))
      ev <- sdm_eval(preds = preds, metric = 'auc', model = model)
      save(ev, file = paste0(file.path(getwd(),spp, 'model_output', 'eval_metrics'), '/',toupper(model), '.RData'))
    }
  } else {
    print(paste0(spp, '- Evaluating Model - ', Sys.time()))
    load(paste0(file.path(getwd(),spp, 'model_output', 'preds'), '/',toupper(model), '.RData'))
    ev <- sdm_eval(preds = preds, metric = 'auc', model = model)
    save(ev, file = paste0(file.path(getwd(),spp, 'model_output', 'eval_metrics'), '/',toupper(model), '.RData'))
  }
    
  if(skip){
    if(file.exists(paste0(file.path(getwd(),spp, 'model_output', 'importance'), '/',toupper(model), '.RData'))){
      print('importance exists and skip == T, so skipping this file!')
      return(NA)
    } else {
      print(paste0(spp, '- Getting Variable Importance - ', Sys.time()))
      load(paste0(file.path(getwd(),spp, 'model_output', 'models'), '/', toupper(model), '.RData'))
      imp <- sdm_importance(mod = mod, se = dfC, pa_col = 'value', xy_col = c('x', 'y'), month_col = 'month', year_col = 'year', model = model)
      save(imp, file = paste0(file.path(getwd(),spp, 'model_output', 'importance'), '/',toupper(model), '.RData'))
    }
  } else {
    print(paste0(spp, '- Getting Variable Importance - ', Sys.time()))
    load(paste0(file.path(getwd(),spp, 'model_output', 'models'), '/', toupper(model), '.RData'))
    imp <- sdm_importance(mod = mod, se = dfC, pa_col = 'value', xy_col = c('x', 'y'), month_col = 'month', year_col = 'year', model = model)
    save(imp, file = paste0(file.path(getwd(),spp, 'model_output', 'importance'), '/',toupper(model), '.RData'))
  }
    

  #} #end y 
    return(summary(mod))
}

plan(multisession, workers = 2)
#sink(file = 'rasters.log', append = T)
checks <- future_pmap(list(..1 = args$spp, ..2 = args$model, ..3 = args$skip), ~ makeMods(spp = ..1, model = ..2, skip = ..3), .progress = T,  .options = furrr_options(seed = 2025))
#sink()
plan(sequential)

#sdmtmb/brt can't be run together in parallel on container
for(m in models){
  print(m)
  for(x in 1:nrow(spp.list)){
    makeMods(spp = spp.list$Name[x], skip = T, model = m)
  }
}


load('norm_MOM6_082025.RData') #norm
predictMods <- function(spp, model){
  #open log file
  sink(file = file.path(getwd(), 'logs', paste0(model, '_prediction.log')), append = T)
  #sink(file = file.path(getwd(), 'logs', paste0(csvName, '.log')), append = T, type = 'message')
  
  # Ensure the sinks are closed when the function exits, regardless of how it exits.
  on.exit({
    #sink(type = "message")
    sink()
  })
  
  print(Sys.time())
  print(spp)
  
  load(paste0(file.path(getwd(),spp, 'model_output', 'models'), '/', toupper(model), '.RData')) #load model - mod

  #predict model
  print('predicting model...')
  abund <- make_predictions(mod = mod, model = model, rasts = norm, mask = T, bathy_nm = 'bathy', bathy_max = 1000, se = dfC, staticData = 'staticVariables_cropped.RData', xy_col = c('x', 'y'), month_col = 'month', year_col = 'year')
  save(abund, file = paste(file.path(getwd(),spp, 'output_rasters'), toupper(model), '.RData', sep = '/'))
  
  return(paste(file.path(getwd(),spp, 'output_rasters'), toupper(model), '.RData', sep = '/'))
}

##############################

##############################
##### MAKE ENSEMBLE  #########
##############################


for(x in 1:nrow(spp.list)){
  print(Sys.time())
  print(spp.list$Name[x])
  
  print('generating weights...')
  #load in evaluation metrics
  evalFlist <- dir(file.path(getwd(),spp.list$Name[x], 'model_output', 'eval_metrics'), pattern = '.RData')
  eval <- vector(length = length(evalFlist))
  for(y in 1:length(evalFlist)){
    load(evalFlist[y])
    eval[y] <- ev
  }
  
  #generate weights
  weights <- EFHSDM::MakeEnsemble(rmse = eval)
  save(weights, file = paste(file.path(getwd(),spp.list$Name[x], 'model_output'), 'ensemble_weights.RData', sep = '/'))
  
  print('making ensemble...')
  #pull preds 
  predFlist <- dir(file.path(getwd(),spp.list$Name[x], 'model_output', 'preds'), pattern = '.RData')
  pds <- vector(mode = 'list', length = length(predFlist))
  for(y in 1:length(predFlist)){
    load(predFlist[y])
    pds[[y]] <- preds
  }
  
  #make ensemble 
  ens <- make_sdm(model = 'ens', ensembleWeights = weights, ensemblePreds = pds)
  save(ens, file = paste(file.path(getwd(),spp.list$Name[x], 'model_output', 'models'), 'ENSEMBLE.RData', sep = '/'))
  
  print('predicting ensemble...')
  
  abundFlist <- dir(file.path(getwd(),spp.list$Name[x], 'output_rasters'), pattern = '.RData')
  abds <- vector(mode = 'list', length = length(abundFlist))
  for(y in 1:length(abundFlist)){
    load(abundFlist[y])
    abds[[y]] <- abund
  }
  
  
  abund <- make_predictions(model = 'ens', rasts = abds, weights = weights)
  save(abund, file = paste(file.path(getwd(),spp.list$Name[x], 'output_rasters'), 'ENSEMBLE.RData', sep = '/'))
}
##############################