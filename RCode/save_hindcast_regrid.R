### save regridded hindcast MOM6 output
#done on 5/23 - will need to redo once new hindcast is out later this summer?

##load packages
library(ncdf4)
library(raster)

##### get static vars

#get static MOM6 vars to define extent
url_static <-"http://psl.noaa.gov/thredds/dodsC/Projects/CEFI/regional_mom6/cefi_portal/northwest_atlantic/full_domain/hindcast/monthly/raw/r20230520/ocean_static.nc"
ncopendap_static <- nc_open(url_static)
lon <- ncvar_get(ncopendap_static, "geolon")
lat <- ncvar_get(ncopendap_static, "geolat")

e <- extent(min(lon), max(lon), min(lat), max(lat)) #extent
se <- extent(-78,-65, 35,45) #extent to subset to

#generate layer names
md <- expand.grid(1:12, 1993:2019)
MD <- paste(month.abb[md[,1]], md[,2], sep = '.')

############ load data & crop 
#### surface Temperature ####
# Specify the OPeNDAP server URL (using regular grid output)
url <- "http://psl.noaa.gov/thredds/dodsC/Projects/CEFI/regional_mom6/cefi_portal/northwest_atlantic/full_domain/hindcast/monthly/regrid/r20230520/tos.nwa.full.hcast.monthly.regrid.r20230520.199301-201912.nc"
surfaceT <- stack(url)

extent(surfaceT) <- e #set extent

names(surfaceT) <- MD #set names to each layer as Mon.Year

surfaceT <- crop(surfaceT, se) #takes ~30 minutes

#### bottom Temperature ####
# Specify the OPeNDAP server URL (using regular grid output)
url <- "http://psl.noaa.gov/thredds/dodsC/Projects/CEFI/regional_mom6/cefi_portal/northwest_atlantic/full_domain/hindcast/monthly/regrid/r20230520/tob.nwa.full.hcast.monthly.regrid.r20230520.199301-201912.nc"
bottomT <- stack(url)

extent(bottomT) <- e

names(bottomT) <- MD #set names to each layer as Mon.Year

bottomT <- crop(bottomT, se) #takes ~30 minutes

#### surface salinity ####
# Specify the OPeNDAP server URL (using regular grid output)
url <- "http://psl.noaa.gov/thredds/dodsC/Projects/CEFI/regional_mom6/cefi_portal/northwest_atlantic/full_domain/hindcast/monthly/regrid/r20230520/sos.nwa.full.hcast.monthly.regrid.r20230520.199301-201912.nc"
surfaceS <- stack(url)

extent(surfaceS) <- e

names(surfaceS) <- MD #set names to each layer as Mon.Year

surfaceS <- crop(surfaceS, se) #takes ~30 minutes

#### bottom salinity ####
# Specify the OPeNDAP server URL (using regular grid output)
url <- "http://psl.noaa.gov/thredds/dodsC/Projects/CEFI/regional_mom6/cefi_portal/northwest_atlantic/full_domain/hindcast/monthly/regrid/r20230520/sob.nwa.full.hcast.monthly.regrid.r20230520.199301-201912.nc"
bottomS <- stack(url)

extent(bottomS) <- e

names(bottomS) <- MD #set names to each layer as Mon.Year

bottomS <- crop(bottomS, se) #takes ~30 minutes

#### bottom oxygen ####
# Specify the OPeNDAP server URL (using regular grid output)
url <- "http://psl.noaa.gov/thredds/dodsC/Projects/CEFI/regional_mom6/cefi_portal/northwest_atlantic/full_domain/hindcast/monthly/regrid/r20230520/btm_o2.nwa.full.hcast.monthly.regrid.r20230520.199301-201912.nc"
bottomO2 <- stack(url)

extent(bottomO2) <- e

names(bottomO2) <- MD #set names to each layer as Mon.Year

bottomO2 <- crop(bottomO2, se) #takes ~30 minutes

#### surface pH ####
# Specify the OPeNDAP server URL (using regular grid output)
url <- "http://psl.noaa.gov/thredds/dodsC/Projects/CEFI/regional_mom6/cefi_portal/northwest_atlantic/full_domain/hindcast/monthly/regrid/r20230520/phos.nwa.full.hcast.monthly.regrid.r20230520.199301-201912.nc"
surfacepH <- stack(url)

extent(surfacepH) <- e

names(surfacepH) <- MD #set names to each layer as Mon.Year

surfacepH <- crop(surfacepH, se) #takes ~30 minutes

#### net primary productivity - vertically integrated ####
# Specify the OPeNDAP server URL (using regular grid output)
url <- "http://psl.noaa.gov/thredds/dodsC/Projects/CEFI/regional_mom6/cefi_portal/northwest_atlantic/full_domain/hindcast/monthly/regrid/r20230520/wc_vert_int_npp.nwa.full.hcast.monthly.regrid.r20230520.199301-201912.nc"
NPP <- stack(url)

extent(NPP) <- e

names(NPP) <- MD #set names to each layer as Mon.Year

NPP <- crop(NPP, se) #takes ~30 minutes

#### MLD ####
# Specify the OPeNDAP server URL (using regular grid output)
url <- "http://psl.noaa.gov/thredds/dodsC/Projects/CEFI/regional_mom6/cefi_portal/northwest_atlantic/full_domain/hindcast/monthly/regrid/r20230520/MLD_003.nwa.full.hcast.monthly.regrid.r20230520.199301-201912.nc"
MLD <- stack(url)

extent(MLD) <- e

names(MLD) <- MD #set names to each layer as Mon.Year

MLD <- crop(MLD, se) #takes ~30 minutes

#### Mesozooplankton biomass ####
# Specify the OPeNDAP server URL (using regular grid output)
url <- "http://psl.noaa.gov/thredds/dodsC/Projects/CEFI/regional_mom6/cefi_portal/northwest_atlantic/full_domain/hindcast/monthly/regrid/r20230520/mesozoo_200.nwa.full.hcast.monthly.regrid.r20230520.199301-201912.nc"
zoops <- stack(url)

extent(zoops) <- e

names(zoops) <- MD #set names to each layer as Mon.Year

zoops <- crop(zoops, se) #takes ~30 minutes

#### Bottom Aragonite Solubility ####
# Specify the OPeNDAP server URL (using regular grid output)
url <- "http://psl.noaa.gov/thredds/dodsC/Projects/CEFI/regional_mom6/cefi_portal/northwest_atlantic/full_domain/hindcast/monthly/regrid/r20230520/btm_co3_sol_arag.nwa.full.hcast.monthly.regrid.r20230520.199301-201912.nc"
argSol <- stack(url)

extent(argSol) <- e

names(argSol) <- MD #set names to each layer as Mon.Year

argSol <- crop(argSol, se) #takes ~30 minutes

#### Diaz. Phyto PP ####
# Specify the OPeNDAP server URL (using regular grid output)
url <- "http://psl.noaa.gov/thredds/dodsC/Projects/CEFI/regional_mom6/cefi_portal/northwest_atlantic/full_domain/hindcast/monthly/regrid/r20230520/jprod_ndi_new_100.nwa.full.hcast.monthly.regrid.r20230520.199301-201912.nc"

daizPP <- stack(url)

extent(daizPP) <- e

names(daizPP) <- MD #set names to each layer as Mon.Year

daizPP <- crop(daizPP, se) #takes ~30 minutes

#### Small Phyto PP ####
# Specify the OPeNDAP server URL (using regular grid output)
url <- "http://psl.noaa.gov/thredds/dodsC/Projects/CEFI/regional_mom6/cefi_portal/northwest_atlantic/full_domain/hindcast/monthly/regrid/r20230520/jprod_nsmp_new_100.nwa.full.hcast.monthly.regrid.r20230520.199301-201912.nc"
smPP <- stack(url)

extent(smPP) <- e

names(smPP) <- MD #set names to each layer as Mon.Year

smPP <- crop(smPP, se) #takes ~30 minutes

#### Medium Phyto PP ####
# Specify the OPeNDAP server URL (using regular grid output)
url <- "http://psl.noaa.gov/thredds/dodsC/Projects/CEFI/regional_mom6/cefi_portal/northwest_atlantic/full_domain/hindcast/monthly/regrid/r20230520/jprod_nmdp_new_100.nwa.full.hcast.monthly.regrid.r20230520.199301-201912.nc"
mdPP <- stack(url)

extent(mdPP) <- e

names(mdPP) <- MD #set names to each layer as Mon.Year

mdPP <- crop(mdPP, se) #takes ~30 minutes

#### Large Phyto PP ####
# Specify the OPeNDAP server URL (using regular grid output)
url <- "http://psl.noaa.gov/thredds/dodsC/Projects/CEFI/regional_mom6/cefi_portal/northwest_atlantic/full_domain/hindcast/monthly/regrid/r20230520/jprod_nlgp_new_100.nwa.full.hcast.monthly.regrid.r20230520.199301-201912.nc"
lgPP <- stack(url)

extent(lgPP) <- e

names(lgPP) <- MD #set names to each layer as Mon.Year

lgPP <- crop(lgPP, se) #takes ~30 minutes

#### Small Zooplankton Nitrogen ####
# Specify the OPeNDAP server URL (using regular grid output)
url <- "http://psl.noaa.gov/thredds/dodsC/Projects/CEFI/regional_mom6/cefi_portal/northwest_atlantic/full_domain/hindcast/monthly/regrid/r20230520/nsmz_100.nwa.full.hcast.monthly.regrid.r20230520.199301-201912.nc"
smZoo <- stack(url)

extent(smZoo) <- e

names(smZoo) <- MD #set names to each layer as Mon.Year

smZoo <- crop(smZoo, se) #takes ~30 minutes


#### Medium Zooplankton Nitrogen ####
# Specify the OPeNDAP server URL (using regular grid output)
url <- "http://psl.noaa.gov/thredds/dodsC/Projects/CEFI/regional_mom6/cefi_portal/northwest_atlantic/full_domain/hindcast/monthly/regrid/r20230520/nmdz_100.nwa.full.hcast.monthly.regrid.r20230520.199301-201912.nc"
mdZoo <- stack(url)

extent(mdZoo) <- e

names(mdZoo) <- MD #set names to each layer as Mon.Year

mdZoo <- crop(mdZoo, se) #takes ~30 minutes


#### Large Zooplankton Nitrogen ####
# Specify the OPeNDAP server URL (using regular grid output)
url <- "http://psl.noaa.gov/thredds/dodsC/Projects/CEFI/regional_mom6/cefi_portal/northwest_atlantic/full_domain/hindcast/monthly/regrid/r20230520/nlgz_100.nwa.full.hcast.monthly.regrid.r20230520.199301-201912.nc"
lgZoo <- stack(url)

extent(lgZoo) <- e

names(lgZoo) <- MD #set names to each layer as Mon.Year

lgZoo <- crop(lgZoo, se) #takes ~30 minutes

#### Downward POC flux ####
# Specify the OPeNDAP server URL (using regular grid output)
url <- "http://psl.noaa.gov/thredds/dodsC/Projects/CEFI/regional_mom6/cefi_portal/northwest_atlantic/full_domain/hindcast/monthly/regrid/r20230520/epc100.nwa.full.hcast.monthly.regrid.r20230520.199301-201912.nc"
poc <- stack(url)

extent(poc) <- e

names(poc) <- MD #set names to each layer as Mon.Year

poc <- crop(poc, se) #takes ~30 minutes


############

############ SAVE
save(surfaceT, bottomT, surfaceS, bottomS, bottomO2, surfacepH, NPP, MLD, zoops, 
     argSol, daizPP, smPP, mdPP, lgPP, 
     smZoo, mdZoo, lgZoo, poc, file = 'C:/Users/katherine.gallagher/Documents/MOM6/Data/regrid_NEsubset_MOM6output_June2025.RData')

##############
##LOAD and CREATE MONTHLY AVERAGES 
load('C:/Users/katherine.gallagher/Documents/MOM6/Data/regrid_NEsubset_MOM6output_June2025.RData')

varList <- list(surfaceT, bottomT, surfaceS, bottomS, bottomO2, 
                surfacepH, NPP, MLD, zoops, 
                argSol, daizPP, smPP, mdPP, lgPP, 
                smZoo, mdZoo, lgZoo, poc)

##1993 - 2019 
varAvg93_19 <- vector(mode = 'list', length = length(varList))
for(v in 1:length(varList)){
  #complete time series average for each variable
  avg <- NULL
  for(x in 1:12){
    mn <- seq(x, nlayers(varList[[v]]), by = 12) #grab all month xs from timeseries by creating a sequence
    MNS <- subset(varList[[v]], mn) #subset raster brick
    mm <- mean(MNS) #average 
    avg <- abind(as.array(mm), avg, along = 3) #make array and bind together
    #print(x)
  }
  avg <- brick(avg) #turn it back into a rasterbrick
  names(avg) <- month.abb #give names
  extent(avg) <- extent(varList[[v]])
  
  varAvg93_19[[v]] <- avg
  print(v)
}

names(varAvg93_19) <- c('surfaceT', 'bottomT', 'surfaceS', 'bottomS', 
                        'bottomO2', 'surfacepH', 'NPP', 'MLD', 'zoops',  
                        'argSol', 'daizPP', 'smPP', 'mdPP', 'lgPP', 
                        'smZoo', 'mdZoo', 'lgZoo', 'poc')

save(varAvg93_19, file = 'variable_average_list_1993_2019.RData')

##normalize data 
mth <- rep(1:12, times = 27) #27 is for the 27 years of the time series 
normVars <- vector(mode = 'list', length = length(varList))
for(v in 1:length(varList)){
  #subset both lists to appropriate rasterbrick - they are in the same order so this will be the same variable 
  vr <- varList[[v]]
  av <- varAvg93_19[[v]]
  
  #normalize data in loop 
  n <- NULL 
  for(x in 1:nlayers(vr)){
    #subset both rasterbricks 
    subX <- subset(vr,x)
    subA <- subset(av,mth[x])
    
    nm <- (subX - subA) / (subX@data@max - subX@data@min)
    n <- abind(as.array(nm), n, along = 3)
  }
  n <- brick(n) #turn it back into a rasterbrick
  names(n) <- names(varList[[v]]) #give names
  extent(n) <- extent(varList[[v]])
  
  normVars[[v]] <- n
  print(v)
}

names(normVars) <- c('surfaceT', 'bottomT', 'surfaceS', 'bottomS', 
                     'bottomO2', 'surfacepH', 'NPP', 'MLD', 'zoops',  
                     'argSol', 'daizPP', 'smPP', 'mdPP', 'lgPP', 
                     'smZoo', 'mdZoo', 'lgZoo', 'poc')

save(normVars, file = 'variable_normalized_list_1993_2019.RData')

### make distance to shore 
grid.pts <- as.data.frame(rasterToPoints(surfaceT)[,1:2])

#load coastal shape file
library(sf)
coast <- sf::read_sf('C:/Users/katherine.gallagher/Documents/gshhg-shp-2.3.7/GSHHS_shp/f/GSHHS_f_L1.shp', crs = 4326)
geo.box <- c(xmin = -78, xmax = -65, ymin = 35, ymax = 45)
sf_use_s2(FALSE)
coast2 <- st_crop(coast, geo.box)
grid.pts$coast_dist <- dist2land(grid.pts, lon = 'x', lat = 'y', binary = F, bind = F, shapefile = coast)

d2c <- rasterize(x = grid.pts[,1:2], y = surfaceT, field = grid.pts[,3], fun = mean)
save(d2c, file = 'C:/Users/katherine.gallagher/Documents/MOM6/Data/dist_to_coast.RData')
