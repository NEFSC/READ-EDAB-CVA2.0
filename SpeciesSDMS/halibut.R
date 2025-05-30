##summer flounder test for methods workshop

###load libraries
library(ncdf4)
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


####### BUILD DIRECTORY ######
setwd("C:/Users/katherine.gallagher/Documents/SpeciesModels")

#name directory - use common name for species in camel case
dirName <- 'Halibut'

#make directory 
dir.create(file.path(getwd(),dirName), showWarnings = F)

#set working directory to new folder 
setwd(file.path(getwd(),dirName))

#######

#### GET FISHERIES DATA ####

#target species
trg <- 	"HIPPOGLOSSUS HIPPOGLOSSUS" 

dirName 

#open connection - needs vpn
channel <- dbutils::connect_to_database(server="NEFSC_pw_oraprod",uid="KGALLAGHER")

#load create_pa_rast function
source('C:/Users/katherine.gallagher/Documents/Code/create_pa_rast.R')

halibut <- create_pa_rast(target = trg)

#save rasters in directory
sppRast <- halibut[[1]]
save(sppRast, file = 'rasters.RData')

#######

####### GET ENVIRONMENTAL DATA ####
##load MOM6 output 
load('C:/Users/katherine.gallagher/Documents/MOM6/Data/regrid_NEsubset_MOM6output_May2025.RData')

##get bathy & rugosity
#bathymetry
bathy <- raster('C:/Users/katherine.gallagher/Documents/NES_bathy.tiff') #load bathy raster

#rugosity 
bathySpat <- rast(bathy) #convert to terra to use terra function
rug <- terrain(bathy, v = 'TRI', unit = 'degrees') #calculate rugosity (techically Terrain Ruggedness Index)

#put bathy/rugosity on MOM6 grid for predictions 
bathy.grid <- rasterToPoints(bathy)
rug.grid <- rasterToPoints(rug)
#coordinates(bathy.grid) <- ~x + y
#proj4string(bathy.grid) <- CRS("+proj=longlat +datum=WGS84 +no_defs ")

r <- subset(bottomT, 1)

bathyC <- rasterize(x = bathy.grid[,1:2], y = r, field = bathy.grid[,3], fun = mean)
rugC <- rasterize(x = rug.grid[,1:2], y = r, field = rug.grid[,3], fun = mean)

#######

#### MATCH FISHERIES & ENVIRONMENTAL DATA ####
##get data frame 
spDF <- halibut[[2]]

#match data to rasters
#bottom temperature - needs to be done by layer
spDF$botTemp <- NA #create blank column
for(x in 1:dim(bottomT)[3]){
  i <- which(spDF$variable == names(bottomT)[x])#match name of raster to variable (raster name from spp stack) 
  spDF$botTemp[i] <- extract(subset(bottomT,x), spDF[i,c('x','y')], fun = mean)
  print(x)
}

#bathymetry
spDF$bathy <- raster::extract(bathy, spDF[,1:2]) #extract values 
#rugosity
spDF$rugosity <- terra::extract(rug, spDF[,1:2]) #extract values  

#subset spDF by removing 0s (1 = surveyed but target spp not found, 2 = surveyed and species found)
sppDF2 <- spDF[spDF$value != 0,]
#switch it back to 0/1 now that we only have grid cells that were sampled 
sppDF2$value <- replace(sppDF2$value, sppDF2$value == 1, 0)
sppDF2$value <- replace(sppDF2$value, sppDF2$value == 2, 1)

save(sppDF2, file = 'fisheries_environment_dataframe.RData')
#######

####### RUN MODELS #######
#### EFHSDM - GAM ####
#formula
gam.form <- formula("value ~ s(x,y, bs = 'ds', k=10) + s(botTemp, k=4) + s(bathy, k=4) + s(rugosity, k=4) + s(month.num, k = 4) + s(year, k = 4) ")
#run model
bi.model <- FitGAM(gam.formula = gam.form, data = sppDF2, family.gam = "binomial", select = T, reduce = T)

#cross validation to get RMSE
###get RMSE
bigam.cv <- CrossValidateModel(model = bi.model, data = sppDF2, folds = 10, model.type = "gam")
bigam.preds <- bigam.cv[[1]]
gam.rmse <- RMSE(obs = bigam.preds$abund, pred = bigam.preds$cvpred)

#make abundance rasters for entire timeseries

#get static vars (month/year) to rasters
r <- subset(bottomT, 1)
rlon<-rlat<-r #copy r to rlon and rlat rasters [1]][1]which will contain the longitude and latitude
xy<-xyFromCell(r,1:length(r)) #matrix of logitudes (x) and latitudes(y)
rlon[]<-xy[,1] #raster of longitudes
rlat[]<-xy[,2] #raster of latitides
rMonth <- rYear <- r

abund <- vector(mode = 'list', length = dim(bottomT)[3])
for(x in 1:dim(bottomT)[3]){
  
  mm.year <- strsplit(names(bottomT)[x], split = '[.]')
  
  mm <- match(mm.year[[1]][1],month.abb) 
  rMonth[] <- mm
  yr <- as.numeric(mm.year[[1]][2])
  rYear[] <- yr
  
  sr <- stack(rlon, rlat, rMonth, rYear, subset(bottomT, x), bathyC, rugC)
  names(sr) <- c("x", "y", "month.num", "year", "botTemp", 'bathy', 'rugosity')
  abund[[x]] <- MakeGAMAbundance(model = bi.model, r.stack = sr)
  print(x)
}
names(abund) <- names(bottomT)

###save everything
save(bi.model, bigam.cv, bigam.preds, gam.rmse, abund, file = "GAMoutputs.RData")


#### EFHSDM - MAXENT ####
#run model - no formula needed
mxnt <- FitMaxnet(data = sppDF2, species = 'value', vars = c('x', 'y', 'month.num', 'year', 'botTemp', 'bathy', 'rugosity'), reduce = T)

#cross validate and get RMSE
mxnt.cv <- CrossValidateModel(model = mxnt, data = sppDF2, folds = 10, model.type = "maxnet", species = 'value', scale.preds = F)
mxnt.preds <- mxnt.cv[[1]]
mxnt.rmse <- RMSE(obs = mxnt.preds$abund, pred = mxnt.preds$cvpred)

#make abundance rasters for entire timeseries 
abundMXNT <- vector(mode = 'list', length = dim(bottomT)[3])
for(x in 1:dim(bottomT)[3]){
  
  mm.year <- strsplit(names(bottomT)[x], split = '[.]')
  
  mm <- match(mm.year[[1]][1],month.abb) 
  rMonth[] <- mm
  yr <- as.numeric(mm.year[[1]][2])
  rYear[] <- yr
  
  sr <- stack(rlon, rlat, rMonth, rYear, subset(bottomT, x), bathyC, rugC)
  names(sr) <- c("x", "y", "month.num", "year", "botTemp", 'bathy', 'rugosity')
  abundMXNT[[x]] <- raster(MakeMaxEntAbundance(model = mxnt, maxent.stack = sr, type = 'maxnet'))
  print(x)
}
#abundMXNT <- stack(abundMXNT)
names(abundMXNT) <- names(bottomT)

###save everything
save(mxnt, mxnt.cv, mxnt.preds, mxnt.rmse, abundMXNT, file = "MAXENToutputs.RData")


#### MATEO - RANDOM FOREST SPATIAL INTERPOLATION (RFSI) ####

#convert dataframe to spatial object 
stDF = st_as_sf(sppDF2, coords = c("x", "y"), crs = 4326, agr = "constant")

# fit the RFSI model
rfsi_model <- rfsi(formula = value ~ botTemp + bathy + rugosity + year + month.num,
                   data = stDF, 
                   zero.tol = 0,
                   n.obs = 10, # number of nearest observations
                   cpus = 4,
                   progress = TRUE,
                   importance = "impurity",
                   seed = 42,
                   num.trees = 250,
                   mtry = 5,
                   splitrule = "variance",
                   min.node.size = 5,
                   sample.fraction = 0.95,
                   quantreg = FALSE)

#cross validation 
## now do cross validation 
# making tgrid - defaults from examples 
n.obs <- 1:6
min.node.size <- 2:10
sample.fraction <- seq(1, 0.632, -0.05) # 0.632 without / 1 with replacement
splitrule <- "variance"
ntree <- 250 # 500
mtry <- 3:(2+2*max(n.obs))
tgrid = expand.grid(min.node.size=min.node.size, num.trees=ntree,
                    mtry=mtry, n.obs=n.obs, sample.fraction=sample.fraction)

rfsiCV <- cv.rfsi(formula = value ~ botTemp + bathy + rugosity + year + month.num,
                  data = stDF, 
                  tgrid = tgrid, # combinations for tuning
                  tgrid.n = 2, # number of randomly selected combinations from tgrid for tuning
                  tune.type = "LLO", # Leave-Location-Out CV
                  k = 10, # number of folds
                  seed = 42,
                  acc.metric = "RMSE", # R2, CCC, MAE
                  output.format = "data.frame", # "data.frame", # "SpatVector",
                  cpus=4, # detectCores()-1,
                  progress=1,
                  importance = "impurity") # ranger parameter

#get RMSE
rf.rmse <- EFHSDM::RMSE(obs = rfsiCV$obs, pred = rfsiCV$pred)

#build predictions 
r <- subset(bottomT, 1)
rMonth <- rYear <- r

abundRF <- vector(mode = 'list', length = dim(bottomT)[3])
for(x in 1:dim(bottomT)[3]){
  
  mm.year <- strsplit(names(bottomT)[x], split = '[.]')
  
  mm <- match(mm.year[[1]][1],month.abb) 
  rMonth[] <- mm
  yr <- as.numeric(mm.year[[1]][2])
  rYear[] <- yr
  
  sr <- stack(rMonth, rYear, subset(bottomT, x), bathyC, rugC)
  names(sr) <- c("month.num", "year", "botTemp", 'bathy','rugosity')
  abundRF[[x]] <- raster(pred.rfsi(model = rfsi_model,
                            data = stDF, 
                            obs.col = "value",
                            newdata = terra::rast(sr), 
                            output.format = "SpatRaster", # "sf", # "SpatVector", 
                            zero.tol = 0,
                            cpus = 1, # detectCores()-1,
                            progress = TRUE))
  print(x)
}
names(abundRF) <- names(bottomT)

###save everything
save(rfsi_model, rfsiCV, rf.rmse, abundRF, file = "RFoutputs.RData")


#### SDMTMB - STATE-SPACE GLMM ####
## WARNING - THIS IS THE LONG ONE - 2/12 deg mesh takes about an hour to run the Cross Validation and ~2 hours to predict the entire timeseries

#make mesh
mesh <- make_mesh(sppDF2, xy_cols = c('x', 'y'), cutoff = 3/12) #using lon/lat since this is on the reprojected regular lat/lon grid, and the domain crosses multiple UTM zones  
#MOM6 resolution is 1/12 = ~8 km 

fit <- sdmTMB(
  value ~ s(botTemp) + s(bathy) + s(rugosity) + s(month.num, k = 4),
  data = sppDF2,
  mesh = mesh,
  family = binomial(link = 'logit'),
  spatial = "off",
  time = 'year',
  spatiotemporal = 'ar1',
  reml = T
)

#cross-validation
plan(multisession)
fitCV <- sdmTMB_cv(
  value ~ s(botTemp) + s(bathy) + s(rugosity) + s(month.num, k = 4),
  data = sppDF2,
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

#change names from sdmTMB preds to make it work with EFHSDM
sdm.preds <- fitCV$data
colnames(sdm.preds)[9:11] <- c('fold', 'pred', 'loglik')

##get abundances
abundSDM <- vector(mode = 'list', length = dim(bottomT)[3])
for(x in 1:dim(bottomT)[3]){ #stopped randomly in the middle so make sure the number is reset
  
  #turn raster into data frame since that is the predict works with sdmTMB
  rastDF <- as.data.frame(rasterToPoints(subset(bottomT,x)))
  rastDF$Mon.Yr <- colnames(rastDF)[3]
  colnames(rastDF)[3] <- 'botTemp'
  #get month/year
  my <- strsplit(x = rastDF$Mon.Yr, split = '[.]')
  rastDF$month <- unlist(lapply(my, FUN = function(x){x[1]}))
  rastDF$month.num <- match(rastDF$month, month.abb)
  rastDF$year <- as.numeric(unlist(lapply(my, FUN = function(x){x[2]})))
  
  ##extract bathy/rug points
  rastDF$bathy <- raster::extract(bathy, rastDF[,c('x', 'y')])
  rastDF$rugosity <- raster::extract(rug, rastDF[,c('x', 'y')])
  
  abundSDM[[x]] <- predict(fit, newdata = rastDF)
  print(x)
}
names(abundSDM) <- names(rast)

###save everything
save(mesh, fit, fitCV, sdm.preds, sdm.rmse, abundSDM, file = "sdmTMBoutputs.RData")



#######

####### CREATE ENSEMBLE #######
enWeights <- MakeEnsemble(rmse = c(gam.rmse, mxnt.rmse, sdm.rmse, rf.rmse)) #make weights 

ens <- ValidateEnsemble(pred.list = list(bigam.preds, mxnt.preds, sdm.preds, rfsiCV), model.weights = enWeights, latlon = F)

#make predictions 
abundENS <- vector(mode = 'list', length = dim(bottomT)[3])
for(x in 1:dim(bottomT)[3]){
  #convert gam/MAXENT/RF to raster
  gamRast <- terra::rast(abund[[x]])
  mxntRast <- terra::rast(abundMXNT[[x]])
  rfRast <- terra::rast(abundRF[[x]])
  
  #manipulate sdmTMB DF back to a raster 
  sdm <- abundSDM[[x]]
  sdm$prob <- exp(sdm$est)/(1+exp(sdm$est))
  coordinates(sdm) <- ~x + y
  proj4string(sdm) <- CRS("+proj=longlat +datum=WGS84 +no_defs ")
  sdmRast <- raster::rasterize(sdm, abund[[x]], 'prob')
  sdmRast <- terra::rast(sdmRast)
  
  abundENS[[x]] <- raster(MakeEnsembleAbundance(model.weights = enWeights, abund.list = list(gamRast, mxntRast, sdmRast, rfRast)))
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
