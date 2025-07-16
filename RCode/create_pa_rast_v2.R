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

#surveyCSV <- c('~/TrawlData/MaineDMR_Trawl_Survey_Tow_Catch_2025-06-30.csv')
#surveyColumns <- matrix(c('Start_Longitude', 'Start_Latitude', 'Start_Date', 'Number_Caught', 'Common_Name'), nrow = 1, ncol = 5)

standardize_data <- function(dataType, channel = channel, csv, csvCols){
  #this function handles all the data pulls and returns a data frame with standard column names to be made into a raster stack
  #dataType is one of the following - 'Surveys', "Observer", or 'CSV'
  #if NOAA Surveys/Obs - will use either survdat or ROracle commands to get data, connection defined with channel 
  #if csv, a file path must be supplied with csv, and relevant columns must be listed in csvCols in the following order: 'longitude', 'latitude', 'date (must be some sort of time vector that can be converted to Posix with as.POSIXct', 'count (can be count/abundance/density, etc)', 'name'
  
  
  ##now move onto different processing pipelines
  if(dataType == 'Surveys'){
    print('Making rasters from Survey Data...')
    ## pull fisheries-independent data
    data <- get_survdat_data(channel, getWeightLength = F, getLengths = F, getBio = F, conversion.factor = T)
    surv <- data$survdat
    rm(data) #we only need surv, so clearing out everything else to save room
    
    #pull species list to help with matching 
    spp.qry <-  paste0("select SCINAME, COMNAME, SVSPP
        from svdbs.SVSPECIES_LIST
        order by SVSPP")
    survSPP <- data.table::as.data.table(DBI::dbGetQuery(channel, spp.qry))
    
    surv <- merge(surv, survSPP)
    surv$MONTH <- month(surv$EST_TOWDATE) #create month column since surv data doesn't come with one 
    dat <- surv[,c('YEAR', "MONTH", "LON", "LAT", "ABUNDANCE", "SCINAME")] #subset to necessary columns 
    names(dat) <- c('year', 'month', 'lon', 'lat', 'count', 'name') #standardize names
  } #end if survey 
  
  if(dataType == 'Observer'){
    print('Making rasters from Observer Data...')
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
    dat <- obsdat[,c('YEAR', "MONTH", 'LONDD', 'LATDD', 'HAILWT', 'SCINAME')] #subset to important columns
    names(dat) <- c('year', 'month', 'lon', 'lat', 'count', 'name') #standardize names
  } #end if observer
  
  if(dataType == 'CSV'){ 
    ### build raster from CSV data from state/other surveys
    print('Making rasters from CSV Data...')
    
    s <- read.csv(csv) #read in each csv
    s$time <- as.POSIXct(s[,csvCols[3]]) #pull time
    s$month <- month(s$time) #make month
    s$year <- year(s$time) #make year
    
    s <- s[,c('time', 'month', 'year', csvCols)] #subset to important columns
    colnames(s) <- c('time', 'month', 'year', 'lon', 'lat', 'date', 'count', 'name')
    dat <- s
    names(dat) <- c('year', 'month', 'lon', 'lat', 'count', 'name') #standardize names
  }
  
  return(dat)
}

create_rast <- function(data, dataType, grid, tmMult = 24 * 60 * 60, origin = '1993-01-01', targetVec){
  #function to create presence/absence (PA) raster brick for a target species
  #data is the dataframe output from standardize_data
  #grid is a ncdcf object with the variables lon, lat, time, and gridVar - can be link to remote data - just must be able to be read with nc_open
      #grid should be **square** and lon/lat should be vectors when pulled 
      #this netcdf should cover the extent of the desired study area
  #gridVar is a variable name in grid - values won't be used, but this is helpful to define any land masks 
  #tmMult - a value to multiply by to get grid timestamp into seconds (for days, this would be 24 * 60 * 60)
  #origin - start of timeseries as YYYY-MM-DD
  #targetVec is a character vector containing the name of the target species and any possible variations 
  
  #subRast - option to subset final raster brick, if TRUE, sub must also be supplied
  #sub is an extent object used to subset raster brick
  
#### step 1 ####
#get grid to map to
gridNC <- nc_open(grid)
#these are vectors of the unique lon/lats on the regridded MOM6 grid
lonR <- ncvar_get(gridNC, "lon") 
latR <- ncvar_get(gridNC, "lat")
tm <- as.POSIXct(ncvar_get(gridNC, 'time') * tmMult, origin = origin)
nc_close(gridNC)

#### step 2 ####
#loop through years & months to build rasters

if(dataType == 'Surveys' | dataType == 'Observer'){
  data <- read.csv(data) 
}

#first, build vector of months/years throughout grid timeseries 
my <- expand.grid(1:12, unique(year(tm)))
my$month.year <- paste(my$Var1, my$Var2, sep = '.')

#create month.year in data
data$month.year <- paste(data$month, data$year, sep = '.')

spRast <- NULL
dataOK <- TRUE

#add check for Fisheries dependent data (dataType = Observer); fisheries independent surveys do not need to do this
dataOK <- TRUE #assume data is good to go
if(dataType == 'Observer'){
  iSPP <- data$name %in% targetVec #was the species caught?  
  
  m <- unique(data$month[iSPP]) #in which months has the species has been caught? 
  
  #per McHenry et al 2019, fisheries dependent data were only considered if the species was caught at least 30 times across at least 6 different months 
  if(length(which(iSPP == T)) >= 30 & length(m) >= 6){ #so if these requirements are met
    print('Observer data meet minimum thresholds...continuing with raster')
  } else {
    print('Observer data do NOT meet minimum thresholds to build raster...make sure that all possible variations of the species name (scientific, common, alternative common names) are included in targetVec')
    dataOK <- FALSE 
  }
} #end if observer

if(dataOK){
  for(x in 1:nrow(my)){
  #create matrix for each month
  mat <- matrix(0, nrow = length(lonR), ncol = length(latR))
  
  #subset observer data to month and year 
  sub <- data[data$month.year == my$month.year[x], ]

  #if there are data
  if(nrow(sub) != 0){
    for(i in 1:nrow(sub)){ 
      
      #find closest lat/lon   
      iLon <- Closest(x = lonR-360, a = sub$lon[i], which = T)
      iLat <- Closest(x = latR, a = sub$lat[i], which = T)
      
      if(sub$name[i] %in% targetVec & sub$count[i] != 0){ #if species is in tow
        mat[iLat, iLon] <- 2 #replace grid cell with a 2
      } else {
        mat[iLat, iLon] <- 1 #if species is not found but grid cell is sampled, 1
      }
    } #end for i
  } #end if nrow(sub)
  
  #add to array
  spRast <- abind(spRast, t(mat), along = 3)
  }
}

#### step 3 #####
#turn into raster brick

#make into rasterBrick
sppBrick <- raster::brick(spRast)
e <- extent(min(lonR), max(lonR), min(latR), max(latR)) #lon/lat from grid
extent(sppBrick) <- e
proj4string(sppBrick) <- CRS('+proj=longlat +datum=WGS84 +no_defs')
sppBrick <- rotate(flip(sppBrick)) #rotate to get orientation right 


#make names for raster layers
names(sppBrick) <- paste(sprintf("%02d", my$Var1), my$Var2, sep = '.') #use unique month.year combinations  

return(sppBrick) #return rasterbrick

} #end function 

merge_rasts <- function(rastList){
  #merge rasters generated from multiple data types 
  #rastList is a list of rasterStacks generated by create_rast
  
  #remove null objects in list 
  rastList <- rastList[!sapply(rastList, is.null)]
  
  rastAll <- rastList[[1]] #make raster to fill (all should have the same extent so it really doesn't matter )
  rastAll[] <- 0 #make it empty 
  for(x in 1:length(rastList)){ #for each member of rastList
    rastAll <- replace(rastAll, rastList[[x]] == 1, 1) #if that member of Rastlist has a 1, replace the value in rastAll with a 1 (indicating that the area was surveyed, but the species was not found)
    rastAll <- replace(rastAll, rastList[[x]] == 2, 2)#if that member of Rastlist has a 2, replace the value in rastAll with a 2 (indicating that the area was surveyed AND the species was found)
  } #this should iterate over each member of the list and update for each layer
  
  return(rastAll)

}