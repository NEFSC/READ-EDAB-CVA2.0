##summer flounder test for methods workshop

###load libraries
library(ncdf4)
library(caret)
library(DescTools)
library(fields)
library(abind)
library(sf)
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


####### BUILD DIRECTORY ######
setwd("~/ClimateVulnerabilityAssessment2.0)/SpeciesModels")

#name directory - use common name for species in camel case
dirName <- 'SummerFlounder'

#make directory 
dir.create(file.path(getwd(),dirName), showWarnings = T) #warning will show up if directory already exists but nothing bad will happen

#set working directory to new folder 
setwd(file.path(getwd(),dirName))

#######

#### GET FISHERIES DATA ####

#target species
trg <- 	"PARALICHTHYS DENTATUS" #summer flounder example for methods workshop

dirName 

#open connection - needs vpn
channel <- dbutils::connect_to_database(server="NEFSC_USERS",uid="KGALLAGHER")

#load create_pa_rast function
source('~/ClimateVulnerabilityAssessment2.0)/Code/create_pa_rast.R')

summerFlounder <- create_pa_rast(target = trg)

#save rasters in directory
sppRast <- summerFlounder[[1]]
save(sppRast, file = 'rasters.RData')

#######

####### GET ENVIRONMENTAL DATA ####
##load normalized MOM6 output 
load("~/ClimateVulnerabilityAssessment2.0)/Data/variable_normalized_list_1993_2019.RData")
#'surfaceT', 'bottomT', 'surfaceS', 'bottomS', 'bottomO2', 'surfacepH', 'NPP', 'MLD', 'zoops',  'argSol', 'daizPP', 'smPP', 'mdPP', 'lgPP', 'smZoo', 'mdZoo', 'lgZoo', 'poc'

#get distance to shore 
load('~/ClimateVulnerabilityAssessment2.0)/Data/dist_to_coast.RData') #d2c - already a raster 

##get bathy & rugosity
#bathymetry
bathy <- raster('~/ClimateVulnerabilityAssessment2.0)/Data/NES_bathy.tiff') #load bathy raster

#rugosity 
bathySpat <- rast(bathy) #convert to terra to use terra function
rug <- terrain(bathy, v = 'TRI', unit = 'degrees') #calculate rugosity (techically Terrain Ruggedness Index)

#put bathy/rugosity on MOM6 grid for predictions 
bathy.grid <- rasterToPoints(bathy)
rug.grid <- rasterToPoints(rug)
#coordinates(bathy.grid) <- ~x + y
#proj4string(bathy.grid) <- CRS("+proj=longlat +datum=WGS84 +no_defs ")

r <- subset(normVars[[1]], 1)

bathyC <- rasterize(x = bathy.grid[,1:2], y = r, field = bathy.grid[,3], fun = mean)
rugC <- rasterize(x = rug.grid[,1:2], y = r, field = rug.grid[,3], fun = mean)

#######

#### MATCH FISHERIES & ENVIRONMENTAL DATA ####
##get data frame 
spDF <- summerFlounder[[2]]
save(spDF, file = 'fisheries_dataframe.RData')

#match data to rasters - do for all covariates and then just use the ones that are important 
for(x in 1:length(normVars)){
  v <- normVars[[x]]
  #nv <- names(normVars)[x]
  
  vc <- vector(length = nrow(spDF))
  #matching by layer
  for(y in 1:dim(v)[3]){
    i <- which(spDF$variable == names(v)[y])#match name of raster to variable (raster name from spp stack) 
    vc[i] <- extract(subset(v,y), spDF[i,c('x','y')], fun = mean) #extract value 
    #print(x)
  }
  spDF <- cbind(spDF, vc)
  print(x)
}
colnames(spDF)[8:ncol(spDF)] <- names(normVars) #name columns 

#match static variables 
#bathymetry
spDF$bathy <- raster::extract(bathy, spDF[,1:2]) #extract values 
#rugosity
spDF$rugosity <- terra::extract(rug, spDF[,1:2]) #extract values  
#distance to coast 
spDF$coast_dist <- raster::extract(d2c, spDF[,1:2]) #extract values 

#subset spDF by removing 0s (0 = not surveyed, 1 = surveyed but target spp not found, 2 = surveyed and species found)
sppDF2 <- spDF[spDF$value != 0,]
#switch it back to 0/1 now that we only have grid cells that were sampled 
sppDF2$value <- replace(sppDF2$value, sppDF2$value == 1, 0)
sppDF2$value <- replace(sppDF2$value, sppDF2$value == 2, 1)

save(sppDF2, file = 'fisheries_environment_dataframe.RData')
#######

#### DETERMINE GUILDS ####
#load in CSV with feeding and habitat guilds defined 
guilds <- read.csv('~/ClimateVulnerabilityAssessment2.0)/CVA2.0 Species List.csv')
g <- which(guilds$SCI_NAME == trg) #find row associated with species

#isolate feeding guild
fGuild <- guilds$Feeding.Guild[g]

#isolate habitat guild
hGuild <- guilds$Habitat.Guild[g]

##isolate covariates for each guild type 
#Feeding guilds: 
  #Planktivores - diazotroph, small, medium, and large phytoplankton primary productivity integrated in top 100 m (new NO3-based) 
  #Piscivores - small, medium, large zooplankton Nitrogen biomass in upper 100 m
  #Bethos/benthivores - Downware particulate organic carbon flux, net primary production? 
  #Apex predators - net primary production

if(fGuild == 'Planktivore'){
  fInd <- which(colnames(sppDF2) == 'daizPP' | 
                  colnames(sppDF2) == 'smPP' | 
                  colnames(sppDF2) == 'mdPP' |
                  colnames(sppDF2) == 'lgPP')
}
if(fGuild == 'Piscivore'){
  fInd <- which(colnames(sppDF2) == 'smZoo' | 
                  colnames(sppDF2) == 'mdZoo' | 
                  colnames(sppDF2) == 'lgZoo')
}
if(fGuild == 'Benthos' | fGuild == 'Benthivore'){
  fInd <- which(colnames(sppDF2) == 'poc' | 
                  colnames(sppDF2) == 'NPP')
}
if(fGuild == 'Apex Predator'){
  fInd <- which(colnames(sppDF2) == 'NPP')
}

#Habitat guilds: 
#  Groundfish/benthos - bottom temperature, salinity, oxygen, aragonite solubility (proxy for pH)
#   Pelagic/migratory - surface temperature, salinity, surface pH, MLD

if(hGuild == 'Groundfish' | hGuild == 'Benthic'){
  hInd <- which(colnames(sppDF2) == 'bottomT' | 
                  colnames(sppDF2) == 'bottomS' | 
                  colnames(sppDF2) == 'bottomO2' |
                  colnames(sppDF2) == 'argSol' )
}
if(hGuild == 'Pelagic' | hGuild == 'Pelagic Migratory'){
  hInd <- which(colnames(sppDF2) == 'surfaceT' | 
                  colnames(sppDF2) == 'surfaceS' | 
                  colnames(sppDF2) == 'surfacepH' |
                  colnames(sppDF2) == 'MLD' )
}

#######

#### REMOVE CORRELATED VARIABLES ####
spDF2 <- sppDF2[,c(1:2,4,6:7, 26:28, hInd, fInd)] #subset dataframe to x,y,value, month, year, static variables, and relevant environmental covariates

corInd <- findCorrelation(cor(spDF2[,-c(1:5)]), names = T) #find correlated variables 

spDF2 <- spDF2[,-which(colnames(spDF2) == corInd)]

#######

####### RUN MODELS #######
#### EFHSDM - GAM ####

#formula - write out static variables that will be consistent across all feeding and habitat guilds
form <- "value ~ s(x,y, bs = 'ts', k=10) + s(bathy, bs = 'ts', k=6) + s(rugosity, bs = 'ts', k=6) + s(coast_dist, bs = 'ts', k=6) + s(month_num, bs = 'cc', k=6) + s(coast_dist, month_num, bs = 'ts', k = 3) + s(year, bs = 'ts', k=6)"
#add covariates
for(x in 6:ncol(spDF2)){
    form <- paste(form, ' + s(', colnames(spDF2)[x], ", bs = 'ts', k=6)", sep = '')
}

gam.form <- formula(form)

#run model
bi.model <- FitGAM(gam.formula = gam.form, data = spDF2, family.gam = "binomial", select = T, reduce = T)

#cross validation to get RMSE
###get RMSE
bigam.cv <- CrossValidateModel(model = bi.model, data = spDF2, folds = 10, model.type = "gam")
bigam.preds <- bigam.cv[[1]]
gam.rmse <- RMSE(obs = bigam.preds$abund, pred = bigam.preds$cvpred)

#make abundance rasters for entire timeseries

#make static vars (month/year) into rasters
r <- subset(bathyC, 1)
rlon<-rlat<-r #copy r to rlon and rlat rasters [1]][1]which will contain the longitude and latitude
xy<-xyFromCell(r,1:length(r)) #matrix of longitudes (x) and latitudes(y)
rlon[]<-xy[,1] #raster of longitudes
rlat[]<-xy[,2] #raster of latitides
rMonth <- rYear <- r

abundGAM <- vector(mode = 'list', length = nlayers(normVars[[1]]))
for(x in 1:nlayers(normVars[[1]])){ #all the rasters in normVars have the same number of layers so it doesn't matter which one we call
  
  mm.year <- strsplit(names(normVars[[1]])[x], split = '[.]') #all the rasterbricks in normVars also have the same names so again, doesn't matter which one we call
  
  mm <- match(mm.year[[1]][1],month.abb) 
  rMonth[] <- mm
  yr <- as.numeric(mm.year[[1]][2])
  rYear[] <- yr
  
  nStack <- vector(mode = 'list', length = length(normVars))
  for(n in 1:length(normVars)){
    nStack[[n]] <- subset(normVars[[n]], x)
  }
  nStack <- stack(nStack)
  crs(nStack) <- crs(bathyC)
  extent(nStack) <- extent(bathyC)
  
  sr <- stack(rlon, rlat, rMonth, rYear, bathyC, rugC, d2c, nStack)
  names(sr) <- c("x", "y", "month_num", "year", 'bathy', 'rugosity', 'coast_dist', names(normVars))
  
  #mask off waters deeper than 1000 m
  sr <- replace(sr, bathyC < -1000 | bathyC > 0, NA)
  
  abundGAM[[x]] <- MakeGAMAbundance(model = bi.model, r.stack = sr)
  print(x)
}
names(abundGAM) <- names(normVars[[1]])

###save everything
save(bi.model, bigam.cv, bigam.preds, gam.rmse, abundGAM, file = "GAMoutputs.RData")


#### EFHSDM - MAXENT ####
#run model - no formula needed
mxnt <- FitMaxnet(data = spDF2, species = 'value', vars = c('x', 'y', 'month_num', 'year', colnames(spDF2)[6:ncol(spDF2)]), reduce = T)

#cross validate and get RMSE
mxnt.cv <- CrossValidateModel(model = mxnt, data = spDF2, folds = 10, model.type = "maxnet", species = 'value', scale.preds = F)
mxnt.preds <- mxnt.cv[[1]]
mxnt.rmse <- RMSE(obs = mxnt.preds$abund, pred = mxnt.preds$cvpred)

#make abundance rasters for entire timeseries 
abundMXNT <- vector(mode = 'list', length = nlayers(normVars[[1]]))
for(x in 1:nlayers(normVars[[1]])){
  
  mm.year <- strsplit(names(normVars[[1]])[x], split = '[.]')
  
  mm <- match(mm.year[[1]][1],month.abb) 
  rMonth[] <- mm
  yr <- as.numeric(mm.year[[1]][2])
  rYear[] <- yr
  
  nStack <- vector(mode = 'list', length = length(normVars))
  for(n in 1:length(normVars)){
    nStack[[n]] <- subset(normVars[[n]], x)
  }
  nStack <- stack(nStack)
  crs(nStack) <- crs(bathyC)
  extent(nStack) <- extent(bathyC)
  
  sr <- stack(rlon, rlat, rMonth, rYear, bathyC, rugC, d2c, nStack)
  names(sr) <- c("x", "y", "month_num", "year", 'bathy', 'rugosity', 'coast_dist', names(normVars))
  
  #mask off waters deeper than 1000 m
  sr <- replace(sr, bathyC < -1000 | bathyC > 0, NA)
  
  abundMXNT[[x]] <- raster(MakeMaxEntAbundance(model = mxnt, maxent.stack = sr, type = 'maxnet'))
  print(x)
}
#abundMXNT <- stack(abundMXNT)
names(abundMXNT) <- names(normVars[[1]])

###save everything
save(mxnt, mxnt.cv, mxnt.preds, mxnt.rmse, abundMXNT, file = "MAXENToutputs.RData")


#### METEO - RANDOM FOREST SPATIAL INTERPOLATION (RFSI) ####

spDF2$staid <- 1:nrow(spDF2) #standin station ids 
#spDF2$MMYY <- my(paste(spDF2$month_num, spDF2$year, sep = '-'))
spDF2$year <- year(spDF2$year)

#convert dataframe to spatial object 
stDF = st_as_sf(spDF2, coords = c("x", "y"), crs = 4326, agr = "constant")
stDF = st_sftime(stDF, time_column_name = 'time')

#create formula 
form <- "value ~ month_num"
#add covariates
for(x in 6:ncol(spDF2)){
  form <- paste(form, ' + ', colnames(spDF2)[x], sep = '')
}

#tuning parameters for model
n.obs <- 5:10
min.node.size <- 2:10
sample.fraction <- seq(1, 0.632, -0.05) # 0.632 without / 1 with replacement
splitrule <- "variance"
ntree <- 200 # 500
mtry <- 3:(2+2*max(n.obs))
tgrid = expand.grid(min.node.size=min.node.size, num.trees=ntree,
                    mtry=mtry, n.obs=n.obs, sample.fraction=sample.fraction)


# cross-validate and tune the RFSI model (all included in single function)
rfsiCV <- cv.rfsi(formula = formula(form),
                   data = stDF,
                   data.staid.x.y.z = c('staid', 'X', 'Y', 'year'),
                   cpus = 10,
                   progress = TRUE,
                   importance = "impurity",
                   seed = 42,
                   acc.metric = 'RMSE',
                   tgrid = tgrid,
                   tgrid.n= length(tgrid),
                   tune.type = "LLO", # Leave-Location-Out CV
                   write.forest = T, 
                   s.crs = st_crs(stDF), 
                   k = 10)

#save best model 
rfsi_model <- rfsiCV$final.model

#get RMSE
rf.rmse <- EFHSDM::RMSE(obs = rfsiCV$obs, pred = rfsiCV$pred)

#build predictions 
r <- subset(normVars[[1]], 1)
rMonth <- rYear <- r

abundRF <- vector(mode = 'list', length = nlayers(normVars[[1]]))
for(x in 1:nlayers(normVars[[1]])){
  
  mm.year <- strsplit(names(normVars[[1]])[x], split = '[.]')
  
  mm <- match(mm.year[[1]][1],month.abb) 
  rMonth[] <- mm
  yr <- as.numeric(mm.year[[1]][2])
  rYear[] <- yr
  
  nStack <- vector(mode = 'list', length = length(normVars))
  for(n in 1:length(normVars)){
    nStack[[n]] <- subset(normVars[[n]], x)
  }
  nStack <- stack(nStack)
  crs(nStack) <- crs(bathyC)
  extent(nStack) <- extent(bathyC)
  
  sr <- stack(rlon, rlat, rMonth, rYear, bathyC, rugC, d2c, nStack)
  names(sr) <- c("x", "y", "month_num", "year", 'bathy', 'rugosity', 'coast_dist', names(normVars))
  
  #mask off waters deeper than 1000 m
  sr <- replace(sr, bathyC < -1000 | bathyC > 0, NA)
  
  abundRF[[x]] <- raster(pred.rfsi(model = rfsi_tune$final.model,
                            data = stDF, 
                            obs.col = "value",
                            newdata = terra::rast(sr), 
                            output.format = "SpatRaster", # "sf", # "SpatVector", 
                            cpus = 10, # detectCores()-1,
                            progress = TRUE))
  print(x)
}
names(abundRF) <- names(normVars[[1]])

###save everything
save(rfsi_model, rfsiCV, rf.rmse, abundRF, file = "RFoutputs.RData")


#### Boosted Regression Tree ####
mod <- gbm.step(data = spDF2, gbm.x = c('x', 'y', 'month_num', 'year', names(spDF2)[6:ncol(spDF2)]), 
                gbm.y = 'value', 
                family = 'bernoulli', tree.complexity = 5,
                learning.rate = 0.005, bag.fraction = 0.75, n.folds = 10)

#k-fold validation code from Camrin Brawn (WHOI): https://zenodo.org/records/7971532
#code for paper Camrin et al 2023: https://esajournals.onlinelibrary.wiley.com/doi/10.1002/eap.2893
source('~/ClimateVulnerabilityAssessment2.0)/Code/eval_kfold_brt_wpredictions.r')
source('~/ClimateVulnerabilityAssessment2.0)/Code/eval_brt.r')
source('~/ClimateVulnerabilityAssessment2.0)/Code/pseudoR2.brt.r')
source('~/ClimateVulnerabilityAssessment2.0)/Code/saveTSS.r')
source('~/ClimateVulnerabilityAssessment2.0)/Code/saveAUC.r')
source('~/ClimateVulnerabilityAssessment2.0)/Code/bhattacharyya.stat.r')
source('~/ClimateVulnerabilityAssessment2.0)/Code/bhatt.coef.r')

brt_kfold <- eval_kfold_brt(dataInput = spDF2, 
                            gbm.x = c('x', 'y', 'month_num', 'year', names(spDF2)[6:ncol(spDF2)]), 
                            gbm.y = 'value',  
                            learning.rate = 0.005, 
                            bag.fraction = 0.75, 
                            tree.complexity = 5, 
                            k_folds = 10,
                            is_fixed = F)

#this is all informative but doesn't provide the information to calculate RMSE
summary_stats <- as.data.frame(eval_brt(mod, spDF2, 'value', plot = FALSE))

brt_kfold2 <- lapply(brt_kfold[[1]], FUN=function(x) x) %>% do.call(rbind, .) %>% as.data.frame()
brt_kfold2 <- data.frame(as.list(colMeans(brt_kfold2)))
for (i in 1:ncol(brt_kfold2)) names(brt_kfold2)[i] <- paste0('kfold_', names(brt_kfold2)[i])
summary_stats <- cbind(summary_stats, brt_kfold2)
summary_stats$kfold_nfolds <- 10
summary_stats$n_pres <- length(which(spDF2[,response] == 1))
summary_stats$n_abs <- length(which(spDF2[,response] == 0))

##edited eval_kfold_brt to provide all training data subsets with predictions/obs to calculate RMSE 
# eval_kfold_brt will now return the summary_stats table as item 1 in the list, 
# and a datatable with the observations and predictions in item 2 of the list
#calculate RMSE
(brt.rmse <- EFHSDM::RMSE(obs = brt_kfold[[2]]$value, pred = brt_kfold[[2]]$preds))

##put obs and predicted together for ensemble 
brt.preds <- data.frame(obs = brt_kfold[[2]]$value, preds = brt_kfold[[2]]$preds)


##abundance 
#get lat/lon from raster
r <- subset(normVars[[1]], 1)
rlon<-rlat<-r #copy r to rlon and rlat rasters [1]][1]which will contain the longitude and latitude
xy<-xyFromCell(r,1:length(r)) #matrix of logitudes (x) and latitudes(y)
rlon[]<-xy[,1] #raster of longitudes
rlat[]<-xy[,2] #raster of latitides
rMonth <- rYear <- r

abundBRT <- vector(mode = 'list', length = nlayers(normVars[[1]]))
for(x in 1:nlayers(normVars[[1]])){
  
  mm.year <- strsplit(names(normVars[[1]])[x], split = '[.]')
  
  mm <- match(mm.year[[1]][1],month.abb) 
  rMonth[] <- mm
  yr <- as.numeric(mm.year[[1]][2])
  rYear[] <- yr
  
  nStack <- vector(mode = 'list', length = length(normVars))
  for(n in 1:length(normVars)){
    nStack[[n]] <- subset(normVars[[n]], x)
  }
  nStack <- stack(nStack)
  crs(nStack) <- crs(bathyC)
  extent(nStack) <- extent(bathyC)
  
  sr <- stack(rlon, rlat, rMonth, rYear, bathyC, rugC, d2c, nStack)
  names(sr) <- c("x", "y", "month_num", "year", 'bathy', 'rugosity', "coast_dist", names(normVars))
  
  #mask off waters deeper than 1000 m
  sr <- replace(sr, bathyC < -1000 | bathyC > 0, NA)
  
  p <- dismo::predict(object = mod, x = sr,type="response")
  abundBRT[[x]] <- raster::rasterize(x = spDF2[,c(1,2)], y = rYear, field = p)
  print(x)
}
names(abundBRT) <- names(normVars[[1]])

###save everything
save(mod, brt.rmse, brt.preds, abundBRT, file = "BRToutputs.RData")


#### SDMTMB - STATE-SPACE GLMM ####
## WARNING - THIS IS THE LONG ONE - 2/12 deg mesh takes about an hour to run the Cross Validation and ~2 hours to predict the entire timeseries

#make mesh
mesh <- make_mesh(spDF2, xy_cols = c('x', 'y'), cutoff = 6/12) #using lon/lat since this is on the reprojected regular lat/lon grid, and the domain crosses multiple UTM zones  
#MOM6 resolution is 1/12 = ~8 km 

form <- "value ~ s(bathy, k=6) + s(rugosity, k=6) + s(month_num, k = 6) + s(year, k = 6) + s(coast_dist, k = 6) + s(month_num, coast_dist, k = 3)"
#add feeding guild covariates
for(x in 1:length(fInd)){
  form <- paste(form, ' + s(', colnames(spDF2)[fInd[x]], ', k = 6) ', sep = '')
}
#add habitat guild covariates
for(x in 1:length(hInd)){
  form <- paste(form, ' + s(', colnames(spDF2)[hInd[x]], ', k = 6) ', sep = '')
}

fit <- sdmTMB(
  formula = formula(form),
  data = spDF2,
  mesh = mesh,
  family = binomial(link = 'logit'),
  spatial = "off",
  time = 'year',
  spatiotemporal = 'ar1',
  reml = T
)
save(mesh, fit, file = "sdmTMBoutputs.RData")


#cross-validation
plan(multisession, workers = 5)
fitCV <- sdmTMB_cv(
  formula = formula(form),
  data = spDF2,
  mesh = mesh,
  family = binomial(link = 'logit'),
  spatial = "off",
  time = 'year',
  spatiotemporal = 'ar1',
  k_folds = 10,
  reml =T
)
plan(sequential)

#calculate RMSE
(sdm.rmse <- EFHSDM::RMSE(obs = fitCV$data$value, pred = fitCV$data$cv_predicted))

save(mesh, fit, fitCV, sdm.rmse, file = "sdmTMBoutputs.RData")

#change names from sdmTMB preds to make it work with EFHSDM
sdm.preds <- fitCV$data
colnames(sdm.preds)[9:11] <- c('fold', 'pred', 'loglik')

##get abundances
abundSDM <- vector(mode = 'list', length = nlayers(normVars[[1]]))
for(x in 1:nlayers(normVars[[1]])){ #stopped randomly in the middle so make sure the number is reset
  
  rS <- vector(mode = 'list', length = length(normVars))
  for(y in 1:length(normVars)){
    rS[[y]] <- subset(normVars[[y]], x)
  }
  rS <- stack(rS)
  names(rS) <- names(normVars)
  
  #turn raster into data frame since that is the predict works with sdmTMB
  rastDF <- as.data.frame(rasterToPoints(rS))
  
  #get month/year
  rastDF$Mon.Yr <- names(normVars[[1]])[x]
  my <- strsplit(x = rastDF$Mon.Yr, split = '[.]')
  rastDF$month <- unlist(lapply(my, FUN = function(x){x[1]}))
  rastDF$month_num <- match(rastDF$month, month.abb)
  rastDF$year <- as.numeric(unlist(lapply(my, FUN = function(x){x[2]})))
  
  ##extract bathy/rug points
  rastDF$bathy <- raster::extract(bathy, rastDF[,c('x', 'y')])
  rastDF$rugosity <- raster::extract(rug, rastDF[,c('x', 'y')])
  rastDF$coast_dist <- raster::extract(d2c, rastDF[,c('x', 'y')])
  
  abundSDM[[x]] <- predict(fit, newdata = rastDF)
  print(x)
}
names(abundSDM) <- names(normVars[[1]])

###save everything
save(mesh, fit, fitCV, sdm.preds, sdm.rmse, abundSDM, file = "sdmTMBoutputs.RData")



#######

####### CREATE ENSEMBLE #######
enWeights <- MakeEnsemble(rmse = c(gam.rmse, mxnt.rmse, sdm.rmse, rf.rmse, brt.rmse)) #make weights 

ens <- ValidateEnsemble(pred.list = list(bigam.preds, mxnt.preds, sdm.preds, rfsiCV, brt.preds), model.weights = enWeights, latlon = F)

#make predictions 
abundENS <- vector(mode = 'list', length = dim(bottomT)[3])
for(x in 1:dim(bottomT)[3]){
  #convert gam/MAXENT/RF/BRT to raster
  gamRast <- terra::rast(abund[[x]])
  mxntRast <- terra::rast(abundMXNT[[x]])
  rfRast <- terra::rast(abundRF[[x]])
  brtRast <- terra::rast(abundBRT[[x]])
  
  #manipulate sdmTMB DF back to a raster 
  sdm <- abundSDM[[x]]
  sdm$prob <- exp(sdm$est)/(1+exp(sdm$est))
  coordinates(sdm) <- ~x + y
  proj4string(sdm) <- CRS("+proj=longlat +datum=WGS84 +no_defs ")
  sdmRast <- raster::rasterize(sdm, abund[[x]], 'prob')
  sdmRast <- terra::rast(sdmRast)
  
  abundENS[[x]] <- raster(MakeEnsembleAbundance(model.weights = enWeights, abund.list = list(gamRast, mxntRast, sdmRast, rfRast, brtRast)))
  print(x)
}
names(abundENS) <- names(bottomT)

save(enWeights, ens, abundENS, file = 'ensemble.RData')


#################
### test plots

#monthly averages
par(mfrow=c(3,4))
for(x in 1:12){
  mn <- seq(x, length(abundENS), by = 12)
  MNS <- abundENS[mn]
  mm <- mean(stack(MNS))
  plot(mm, zlim = c(0,1), main = month.abb[x])
}

#changes across time series 
r93 <- stack(abundENS[1:12])
r19 <- stack(abundENS[313:324])

rDiff <- r19-r93
plot(rDiff, zlim = c(-0.5,1), main = month.abb)

#mask off offshore 
r19.2 <- replace(r19, bathyC < -1000, 0)
