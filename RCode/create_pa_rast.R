##build presence/absence rasters from fisheries independent and dependent surveys 

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

#open connection - needs vpn
#channel <- dbutils::connect_to_database(server="NEFSC_pw_oraprod",uid="KGALLAGHER")

#scientific name of species of interest
#trg <- "MUSTELUS CANIS"  #smooth dogfish

create_pa_rast <- function(target){
  #function to create presence/absence (PA) raster brick for a target species
  #target is the scientific name of the species of interest
  #connection to NEFSC database should be established prior to running the function
#### step 1 ####
#get grid 

#load in regrid lat lons to map to 
url <- "http://psl.noaa.gov/thredds/dodsC/Projects/CEFI/regional_mom6/cefi_portal/northwest_atlantic/full_domain/hindcast/monthly/regrid/r20230520/tob.nwa.full.hcast.monthly.regrid.r20230520.199301-201912.nc"
ncopendap_regrid <- nc_open(url)
#these are vectors of the unique lon/lats on the regridded MOM6 grid
lonR <- ncvar_get(ncopendap_regrid, "lon") 
latR <- ncvar_get(ncopendap_regrid, "lat")
exm <- ncvar_get(ncopendap_regrid, 'tob', start = c(1,1,1), count = c(-1,-1,1)) #pull example slice of regridded data to use to map NAs for land 

###FISHERIES INDEPENDENT DATA 
####pull data
## pull fisheries-independent data
data <- get_survdat_data(channel, getWeightLength = F, getLengths = F, getBio = F, conversion.factor = T)
surv <- data$survdat
rm(data) #clear out object to save room

#pull species list to help with matching 
spp.qry <-  paste0("select SCINAME, COMNAME, SVSPP
      from svdbs.SVSPECIES_LIST
      order by SVSPP")
survSPP <- data.table::as.data.table(DBI::dbGetQuery(channel, spp.qry))

surv <- merge(surv, survSPP)

####FISHERIES DEPENDENT DATA
## pull fisheries-dependent data - a bit more intense since there isn't a nice function to do it, and needs to be subset, but follows the same basic steps as survdat
#observer data 
obs.qry <-  paste0("select YEAR, MONTH, TRIPID, HAULNUM, LONHBEG, LATHBEG, NESPP4, HAILWT
      from obdbs.OBSPP
      where YEAR between 1993 and 2019
      order by YEAR, MONTH, TRIPID")
obs <- data.table::as.data.table(DBI::dbGetQuery(channel, obs.qry))

#at sea monitor data
asm.qry <- paste0("select YEAR, MONTH, TRIPID, HAULNUM, LONHBEG, LATHBEG, NESPP4, HAILWT
      from obdbs.ASMSPP
      where YEAR between 1993 and 2019
      order by YEAR, MONTH, TRIPID")
asm <- data.table::as.data.table(DBI::dbGetQuery(channel, asm.qry))

#combine them 
ob.asm <- rbind(obs, asm)

#species query
spp.qry <- paste0("select NESPP4, COMNAME, SCINAME
      from obdbs.OBSPEC")
obsSPP <- data.table::as.data.table(DBI::dbGetQuery(channel, spp.qry))

#merge datasets 
obsdat <- merge(ob.asm, obsSPP, by = 'NESPP4')
#remove NAs
obsdat <- na.omit(obsdat)

#fix formatting for YEAR & MONTH (other variables could/should be fixed too but these are the ones we're using at the moment)
obsdat$YEAR <- as.numeric(obsdat$YEAR)
obsdat$MONTH <- as.numeric(obsdat$MONTH)

###converting position to decimal degrees 
#split out components
#longitude
obsdat$LON_DEG <- substr(obsdat$LONHBEG, 1, 2)
obsdat$LON_MIN <- substr(obsdat$LONHBEG, 3, 4)
obsdat$LON_SEC <- substr(obsdat$LONHBEG, 5, 6)
#convert to decimal degrees 
obsdat$LONDD <- as.numeric(conv_unit(paste(obsdat$LON_DEG, obsdat$LON_MIN, obsdat$LON_SEC, sep = ' '), from = 'deg_min_sec', to = 'dec_deg'))  
obsdat$LONDD <- obsdat$LONDD * -1 #multiply by negative 1 to get W

#latitude
obsdat$LAT_DEG <- substr(obsdat$LATHBEG, 1, 2)
obsdat$LAT_MIN <- substr(obsdat$LATHBEG, 3, 4)
obsdat$LAT_SEC <- substr(obsdat$LATHBEG, 5, 6)
#convert to decimal degrees
obsdat$LATDD <- as.numeric(conv_unit(paste(obsdat$LAT_DEG, obsdat$LAT_MIN, obsdat$LAT_SEC, sep = ' '), from = 'deg_min_sec', to = 'dec_deg'))

##subset to NE 
i <- which(obsdat$LONDD >= -78 & obsdat$LONDD <= -65 & obsdat$LATDD >= 35 & obsdat$LATDD <= 45)
obsdatNE <- obsdat[i,]
###from data pull to here takes ~12 minutes; NE subset is much easier to work with 
rm(obsdat) #remove larger dataframe to save space/memory

#### step 2 ####
#loop through years and months to build rasters for target species

yrs <- 1993:2019
mths <- 1:12

#### FISHERIES DEPENDENT DATA
#check if we can use the fisheries dependent data, using the same thresholds as McHenry et al 2019
oInd <- obsSPP$NESPP4[which(obsSPP$SCINAME == trg)]
iSPP <- obsdatNE$NESPP4 %in% oInd #was the species caught?  

m <- unique(obsdatNE$MONTH[iSPP]) #in which months has the species has been caught? 

#per McHenry et al 2019, fisheries dependent data were only considered if the species was caught at least 30 times across at least 6 different months 
sppDep <- NULL
if(length(which(iSPP == T)) >= 30 & length(m) >= 6){
  
  for(y in yrs){
    for(mn in mths){
      #create matrix for each month
      mat <- matrix(0, nrow = length(lonR), ncol = length(latR))
      
      #subset observer data to month and year 
      sub <- obsdatNE[obsdatNE$YEAR == y & obsdatNE$MONTH == mn, ]
      sub$ID <- paste0(sub$TRIPID, sub$HAULNUM) #create unique ids to loop through 
      ids <- unique(sub$ID)
      #if there are data
      if(nrow(sub) != 0){
        for(i in ids){ #for each unique trip and haul
          tow <- sub[sub$ID == i, ]
          
          #find closest lat/lon - using just first entry because they should all be the same  
          iLon <- Closest(x = lonR-360, a = tow$LONDD[1], which = T)
          iLat <- Closest(x = latR, a = tow$LATDD[1], which = T)
          
          if(trg %in% tow$SCINAME){ #if species is in tow
            mat[iLat, iLon] <- 2 #replace grid cell with a 2
          } else {
            mat[iLat, iLon] <- 1 #if species is not found but grid cell is sampled, 1
          }
        } #end for i
      } #end if nrow(sub)
      
      #put in NAs for land following example file 
      #mat <- replace(mat, is.na(t(exm)), NA)
      
      #add to array
      sppDep <- abind(sppDep, t(mat), along = 3)
    } #end mn
    print(y)
  } #end y
 
} #end if species was detected enough


#### FISHERIES INDEPENDENT DATA
surv$MONTH <- month(surv$EST_TOWDATE) #create month column since surv data doesn't come with one 

sppInd <- NULL #initialize array 
  
for(y in yrs){
    for(mn in mths){
      #create matrix for each month
      mat <- matrix(0, nrow = length(lonR), ncol = length(latR))
      
      #subset observer data to month and year 
      sub <- surv[surv$YEAR == y & surv$MONTH == mn, ]
      sub$ID <- paste0(sub$CRUISE6, sub$STRATUM, sub$TOW, sub$STATION) #create unique ids to loop through 
      ids <- unique(sub$ID)
      #if there are data
      if(nrow(sub) != 0){
        for(i in ids){ #for each unique trip and haul
          tow <- sub[sub$ID == i, ]
          
          #find closest lat/lon - using just first entry because they should all be the same  
          iLon <- Closest(x = lonR-360, a = tow$LON[1], which = T)
          iLat <- Closest(x = latR, a = tow$LAT[1], which = T)
          
          if(trg %in% tow$SCINAME){ #if species is in tow
            mat[iLon, iLat] <- 2 #replace grid cell with a 2
          } else {
            mat[iLon, iLat] <- 1 #1 if area was surveyed, but species was absent
          }
        } #end for i
      } #end if nrow(sub)
      
      #put in NAs for land following example file 
      mat <- replace(mat, is.na(exm), NA)
      
      #add to array
      sppInd <- abind(sppInd, t(mat), along = 3)
    } #end mn
    print(y)
  } #end y

#### step 3 #####
# combine independent and dependent data
if(is.null(sppDep) == F){ #if fisheries dependent data exists
  sppAll <- sppDep + sppInd
} else {
  sppAll <- sppInd #sppAll is just equal to the fisheries independent data
}
sppAll <- replace(sppAll, sppAll > 2, 2) #set any values greater 2 (there aren't many) - indicating that the target species was found in both independent and dependent surveys within the same month - to 1

#make into rasterBrick
sppBrick <- raster::brick(sppAll)
e <- extent(min(lonR), max(lonR), min(latR), max(latR)) #lon/lat from MOM6 regrided above
extent(sppBrick) <- e
proj4string(sppBrick) <- CRS('+proj=longlat +datum=WGS84 +no_defs')
sppBrick <- rotate(flip(sppBrick)) #rotate to get orientation right 

#subset to NE
se <- extent(-78,-65, 35,45) 
sppBrickNE <- crop(sppBrick, se) 

#make names for raster layers
nm <- expand.grid(yrs, month.abb[mths])
nm <- nm[order(nm$Var1),] #order by year
nms <- paste(nm$Var2, nm$Var1, sep = '-')
names(sppBrickNE) <- nms

#### convert to data frame for modeling 
sppDF <- as.data.frame(rasterToPoints(sppBrickNE))
sppDF <- melt(sppDF, id.vars = 1:2)

#split out month and year
sppDF$variable <- as.character(sppDF$variable)
my <- strsplit(x = sppDF$variable, split = '[.]')
sppDF$month <- unlist(lapply(my, FUN = function(x){x[1]}))
sppDF$month.num <- match(sppDF$month, month.abb)
sppDF$year <- as.numeric(unlist(lapply(my, FUN = function(x){x[2]})))




return(list(sppBrickNE, sppDF)) #return both raster brick and data.frame

} #end function 