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

############

############ SAVE
save(surfaceT, bottomT, surfaceS, bottomS, bottomO2, surfacepH, NPP, MLD, zoops, file = 'C:/Users/katherine.gallagher/Documents/MOM6/Data/regrid_NEsubset_MOM6output_May2025.RData')


