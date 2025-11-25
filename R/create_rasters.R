#' @title Fisheries Data Functions
#' @description Functions to standardize fisheries independent and dependent datasets and build and combine rasters of presence/absence for target species
#' \itemize{
#' \item \code{standardize_data} pulls NEFSC trawl and observer data from the NEFSC databases, or opens a local CSV file, and standardizes the output for \code{create_rast}
#' \item \code{create_rast} builds a raster on the grid provided of effort and presence/absences for the target species from the provided source data
#' \item \code{merge_rasts} combines rasters from \code{create_rast} across multiple sources
#' }
#' @param dataType type of data to be standardize. Must be one of the following: 'NESurveys', 'NEObserver', or 'CSV'
#' @param channel connection to remote databases. Only required for 'NESurveys' or 'NEObserver'
#' @param csv path to local CSV file
#' @param csvCols Column names in csv file in the following order: 'towid', 'longitude', 'latitude', 'date', 'count (can be count/abundance/density, etc)', 'name'. Date must be in a format that can be converted to POSIX with as.POSIXct
#' @param yrRange A vector with length of 2 indicating the start and end, inclusive, year of the desired time series
#' @param data data frame to convert to presence/absence raster
#' @param isObs TRUE/FALSE indicating whether or not the data is observer or similar fisheries-dependent data. Will force function to check if species is observed at least 30 times throughout timeseries before creating raster
#' @param grid static link to a ncdcf object with the variables lon, lat, time - can be link to remote data - must be able to be read with nc_open
#' @param tmMult multiplier to help convert timestep to POSIX (seconds since origin), defaults to 86400 (number of seconds in a day)
#' @param origin Origin of time series
#' @param targetVec a vector containing all possible names for the target species. Must have a length >= 1
#' @param rastList list of outputs from \code{create_rast} to combine
#' @return \code{standardize_data} returns a data frame. It is recommended to save this as a csv file as pulling survey and observer datasets does take time
#' @return \code{create_rast} returns a rasterBrick with the same extent as the provided grid, and a number of layers equal to the timeseries associated with the provided model data. Values in the rasterBrick will be 0 for no effort in the cell, 1 for fishing effort but no catch, or 2 for fishing effort with catch of target species 
#' @return \code{merge_rasts} returns a single rasterBrick with the same resolution and range as the \code{create_rast} outputs

standardize_data <- function(dataType, channel = channel, csv, csvCols, yrRange){

if(dataType %in% c('NESurveys', 'NEObserver', 'CSV')){
  
  ##now move onto different processing pipelines
  if(dataType == 'Surveys'){
    print('Standardizing Survey Data...')
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
    surv$towID <- paste(surv$CRUISE6, surv$STRATUM, surv$TOW, surv$STATION)
    dat <- surv[,c('towID', 'YEAR', "MONTH", "LON", "LAT", "ABUNDANCE", "SCINAME")] #subset to necessary columns 
    names(dat) <- c('towID', 'year', 'month', 'lon', 'lat', 'count', 'name') #standardize names
    dat <- dat[dat$year >= yrRange[1] & dat$year <= yrRange[2],]
  } #end if survey 
  
  if(dataType == 'Observer'){
    print('Standardizing Observer Data...')
    ## pull fisheries-dependent data - a bit more intense since there isn't a nice function to do it, and needs to be subset, but follows the same basic steps as survdat
    #observer data 
    obs.qry <-  paste0("select YEAR, MONTH, TRIPID, HAULNUM, LONHBEG, LATHBEG, NESPP4, HAILWT
        from obdbs.OBSPP
        where YEAR between ", yrRange[1], " and ", yrRange[2], "
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
    
    obsdat$ID <- paste0(obsdat$TRIPID, obsdat$HAULNUM)
    
    dat <- obsdat[,c('ID', 'YEAR', "MONTH", 'LONDD', 'LATDD', 'HAILWT', 'SCINAME')] #subset to important columns
    names(dat) <- c('towID', 'year', 'month', 'lon', 'lat', 'count', 'name') #standardize names
  } #end if observer
  
  if(dataType == 'CSV'){ 
    ### build raster from CSV data from state/other surveys
    print('Standardizing CSV Data...')
    
    s <- read.csv(csv) #read in each csv
    s$time <- as.POSIXct(s[,csvCols[4]]) #pull time
    s$month <- month(s$time) #make month
    s$year <- year(s$time) #make year
    
    dat <- s[,c('time', 'month', 'year', csvCols)] #subset to important columns
    colnames(dat) <- c('time', 'month', 'year', 'towID', 'lon', 'lat', 'date', 'count', 'name')
	dat <- dat[dat$year >= yrRange[1] & dat$year <= yrRange[2],]
  }

  return(dat)

  }  else {
    stop('dataType is not NESurveys, NEObserver, or CSV')
}

}

create_rast <- function(data, isObs = FALSE, grid, tmMult = 24 * 60 * 60, origin = '1993-01-01', targetVec){
  
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

#first to make this a little bit easier, subset the dataset to have the same year extent as the grid
tmInd <- which(data$year >= min(year(tm)) & data$year <= max(year(tm)))
data <- data[tmInd,]

#first, build vector of months/years throughout grid timeseries 
my <- expand.grid(1:12, unique(year(tm)))
my$month.year <- paste(my$Var1, my$Var2, sep = '.')

#create month.year in data
data$month.year <- paste(data$month, data$year, sep = '.')

spRast <- NULL
dataOK <- TRUE

#add check for Fisheries dependent data (dataType = Observer); fisheries independent surveys do not need to do this
dataOK <- TRUE #assume data is good to go
if(isObs == TRUE){
  iSPP <- data$name %in% targetVec #was the species caught?  
  
 # m <- unique(data$month[iSPP]) #in which months has the species has been caught? 
  
  #per McHenry et al 2019, fisheries dependent data were only considered if the species was caught at least 30 times across at least 6 different months
	#we lightened that threshold to just 30 observations to account for the fact that some observer programs only go for 5 months of the year 
  #if(length(which(iSPP == T)) >= 30 & length(m) >= 6){ #so if these requirements are met
  if(length(which(iSPP == T)) >= 30){
    print('Observer data meet minimum thresholds...continuing with raster')
  } else {
    print('Observer data do NOT meet minimum thresholds to build raster...make sure that all possible variations of the species name (scientific, common, alternative common names) are included in targetVec')
    dataOK <- FALSE 
  }
} #end if observer

if(dataOK){
  for(x in 1:nrow(my)){
  #create matrix for each month
  mat <- matrix(0, nrow = length(latR), ncol = length(lonR))
  
  #subset observer data to month and year 
  sub <- data[data$month.year == my$month.year[x], ]

  #if there are data
  if(nrow(sub) != 0){
    ids <- unique(sub$towID)
    for(i in ids){ 
      tow <- sub[sub$towID == i, ]
      
      #find closest lat/lon   
      iLon <- Closest(x = lonR-360, a = tow$lon[1], which = T)
      iLat <- Closest(x = latR, a = tow$lat[1], which = T)
      
      if(any(targetVec %in% tow$name)){ #if species is in tow
        mat[iLat, iLon] <- 2 #replace grid cell with a 2
      } else {
        mat[iLat, iLon] <- 1 #if species is not found but grid cell is sampled, 1
      }
    } #end for i
  } #end if nrow(sub)
  
  #add to array
  spRast <- abind(spRast, mat, along = 3)
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

} else {
  sppBrick <- NULL
}

return(sppBrick) #return rasterbrick

} #end function 

merge_rasts <- function(rastList){
  require(raster)
  #merge rasters generated from multiple data types 
  #rastList is a list of rasterStacks generated by create_rast
  
  #remove null objects in list 
  i <- unlist(lapply(rastList,is.null))
  if(length(which(i)) != 0){
    rastList <- rastList[!i]
  }
  
  rastAll <- vector(mode = 'list', length = raster::nlayers(rastList[[1]]))
  for(x in 1:raster::nlayers(rastList[[1]])){
    rs <- stack(lapply(rastList, FUN = function(y){subset(y, x)}))
    rastAll[[x]] <- calc(rs, max) 
  }
  rastAll <- raster::brick(rastAll)
  extent(rastAll) <- extent(rastList[[1]])
  proj4string(rastAll) <- CRS('+proj=longlat +datum=WGS84 +no_defs')
  names(rastAll) <- names(rastList[[1]])
  return(rastAll)
  
}
