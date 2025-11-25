### function to pull MOM6 environmental data and create a list of normalized data rasters 

########## INPUTS 
#varURL <- "https://psl.noaa.gov/cefi_portal/data_index/cefi_data_indexing.Projects.CEFI.regional_mom6.cefi_portal.northwest_atlantic.full_domain.hindcast.json" #url pointing to variable lists for desired run type (hindcast, decadal forecast) and domain (nwa) in question
#reqVars <- c('Bottom Temperature') #list of vars to pull 
#gt <- 'regrid' #desired grid type (options = regrid or raw)
#of <- 'monthly' #desired output frequency (options = daily or monthly)
#urlHead <- "https://psl.noaa.gov/thredds/ncss/grid/Projects/CEFI/regional_mom6/" #first part of the URL to access ncss server to allow for subsetting 
##############

pull_hind <- function(varURL, reqVars, shortNames, gt = 'regrid', of = 'monthly', bounds = c(-78,-65, 35,45), static = "http://psl.noaa.gov/thredds/dodsC/Projects/CEFI/regional_mom6/cefi_portal/northwest_atlantic/full_domain/hindcast/monthly/raw/r20230520/ocean_static.nc", release){
  require(jsonlite)
  require(ncdf4)
  
  vars <- fromJSON(varURL) #turn json file into a list 
  
  long.name <- url <- grid.type <- out.freq <- rl <- NULL  #pull the long names, full opendap urls, grid types, and output frequency for indexing which files to pull 
  for(x in 1:length(vars)){
    long.name <- c(long.name, vars[[x]]$cefi_long_name)
    grid.type <- c(grid.type, vars[[x]]$cefi_grid_type)
    out.freq <- c(out.freq, vars[[x]]$cefi_output_frequency)
    url <- c(url, vars[[x]]$cefi_opendap)
    rl <- c(rl, vars[[x]]$cefi_release)
  }
  
  varList <- vector(mode = 'list', length = length(reqVars)) #initalize empty lists to store all the data 
  
  #get info for subsetting
  #putting subsetting back because everything else takes too long otherwise 
  stat <- nc_open(static)
  lon <- ncvar_get(stat, "geolon")
  lat <- ncvar_get(stat, "geolat")
  
  e <- extent(min(lon), max(lon), min(lat), max(lat)) #extent
  se <- extent(bounds) #extent to subset to
  
  for(y in 1:length(reqVars)){
    ind <- which(long.name == reqVars[y] & grid.type == gt & out.freq == of  & rl == release) #find appropriate url for the variable 
    
    #load url 
    v <- stack(url[ind])

    #create and set names
    n <- matrix(unlist(strsplit(names(v), split = '[.]')), ncol =3, nrow = nlayers(v), byrow = T)
    n[,1] <- gsub('X', replacement = '', n[,1])
    
    names(v) <- paste(n[,2], n[,1], sep = '.') #set names 
    extent(v) <- e #set extent
    #subset
    v <- crop(v, se) #this is the rate limiting step
    
    varList[[y]] <- v #save raw data in list 
  }
  names(varList) <- shortNames
  return(varList)
}

avg_env <- function(rastList){

  avgList <- vector(mode = 'list', length = length(rastList))
  
  for(x in 1:length(rastList)){
    v <- rastList[[x]]
    
    ## create monthly averages 
    avgs <- NULL 
    for(m in 1:12){
      mn <- seq(m, nlayers(v), by = 12) #grab all month xs from timeseries by creating a sequence
      MNS <- subset(v, mn) #subset raster brick
      mm <- mean(MNS) #average 
      avgs <- abind(as.array(mm), avgs, along = 3) #make array and bind together
      #print(m)
    } #end m
    #convert to rasterBrick
    
    avgs <- brick(avgs)
    extent(avgs) <- extent(v) #make the extent the same as v
    crs(avgs) <- crs(v)
    names(avgs) <- 1:12
    avgList[[x]] <- avgs

  } #end x
  
  names(avgList) <- names(rastList)
  return(avgList)
} #end function

sd_env <- function(rastList){
  
  sdList <- vector(mode = 'list', length = length(rastList))
  
  for(x in 1:length(rastList)){
    v <- rastList[[x]]
    
    ## create monthly averages 
    sds <- NULL 
    for(m in 1:12){
      mn <- seq(m, nlayers(v), by = 12) #grab all month xs from timeseries by creating a sequence
      MNS <- raster::subset(v, mn) #subset raster brick
      mm <- raster::calc(MNS, fun = sd, na.rm = T) #average 
      sds <- abind::abind(as.array(mm), sds, along = 3) #make array and bind together
      #print(m)
    } #end m
    #convert to rasterBrick
    
    sds <- brick(sds)
    extent(sds) <- extent(v) #make the extent the same as v
    crs(sds) <- crs(v)
    names(sds) <- 1:12
    sdList[[x]] <- sds
    
  } #end x
  
  #names(sdList) <- names(rastList)
  return(sdList)
} #end function
    
norm_env <- function(rawList, avgList, sdList, shortNames){
 normList <- vector(mode = 'list', length = length(rawList))
  
  for(x in 1:length(rawList)){
    v <- rawList[[x]]
    avgs <- avgList[[x]]
    sds <- sdList[[x]]
    
    ##normalize data to monthly averages 
    mth <- rep(1:12, times = nlayers(v)/12) #creates repeating list of 1:12 for each year
    #normalize data
    norm <- NULL 
    for(m in 1:nlayers(v)){
      #subset both rasterbricks 
      subX <- subset(v,m)
      subA <- subset(avgs,mth[m])
      subS <- subset(sds,mth[m])
      
      nm <- (subX - subA) / subS
      norm <- abind(as.array(nm), norm, along = 3)
    } #end m
    norm <- brick(norm)
    extent(norm) <- extent(v) #make the extent the same as v
    crs(norm) <- crs(v)
    names(norm) <- names(v) #set names 
    normList[[x]] <- norm #save
  } #end x 
  
  names(normList) <- shortNames #set names of lists 
  
  return(normList)

} #end function 


###pull forecast data

pull_forecast <- function(varURL, reqVars, shortNames, gt = 'regrid', of = 'monthly', bounds = c(-78,-65, 35,45), static = "http://psl.noaa.gov/thredds/dodsC/Projects/CEFI/regional_mom6/cefi_portal/northwest_atlantic/full_domain/hindcast/monthly/raw/r20230520/ocean_static.nc", release, init){
  require(jsonlite)
  require(ncdf4)
  
  vars <- fromJSON(varURL) #turn json file into a list 
  
  long.name <- url <- grid.type <- out.freq <- rl <- init.date <- NULL  #pull the long names, full opendap urls, grid types, and output frequency for indexing which files to pull 
  for(x in 1:length(vars)){
    long.name <- c(long.name, vars[[x]]$cefi_long_name)
    grid.type <- c(grid.type, vars[[x]]$cefi_grid_type)
    out.freq <- c(out.freq, vars[[x]]$cefi_output_frequency)
    url <- c(url, vars[[x]]$cefi_opendap)
    rl <- c(rl, vars[[x]]$cefi_release)
      init.date <- c(init.date, vars[[x]]$cefi_init_date)
  }
  
  varList <- vector(mode = 'list', length = length(reqVars)) #initalize empty lists to store all the data 
  
  #get info for subsetting
  #putting subsetting back because everything else takes too long otherwise 
  stat <- nc_open(static)
  lon <- ncvar_get(stat, "geolon")
  lat <- ncvar_get(stat, "geolat")
    nc_close(stat)
  
  e <- extent(min(lon), max(lon), min(lat), max(lat)) #extent
  se <- extent(bounds) #extent to subset to
  
  for(y in 1:length(reqVars)){
    ind <- which(long.name == reqVars[y] & grid.type == gt & out.freq == of  & rl == release & init.date == init) #find appropriate url for the variable 
    
    #load url with netcdf to account for ensemble members
      var <- NULL
    r <- nc_open(url[ind])
      tm <- ncvar_get(r, 'lead')
      for(m in 1:10){
          vm <- NULL
          for(z in 1:length(tm)){
              v <- ncvar_get(r, names(r$var), start = c(1, 1, z, m), count = c(-1, -1, 1, 1))
              vm <- abind(vm, v, along = 3)
            } #end z
          var <- abind(var, vm, along = 4)
        } #end m
      nc_close(r)

      ##take average of ensemble members
      varAvg <- apply(var, MARGIN = c(1:3), FUN = mean, na.rm = T)

        #flip to get orientation right
      varFlip <- aperm(varAvg, c(2, 1, 3))
        varFlip <- varFlip[nrow(varFlip):1, , ]
      
      #convert to raster
      v <- brick(varFlip)
      extent(v) <- e
      crs(v) <- "+proj=longlat +datum=WGS84 +no_defs"

    #create and set names
    yr <- as.numeric(substr(init, 2, 5))
      yr10 <- yr+9
      nms <- expand.grid(1:12, yr:yr10)
    
    names(v) <- paste(nms[,1], nms[,2], sep = '.') #set names 
    extent(v) <- e #set extent
    #subset
    v <- crop(v, se) #this is the rate limiting step
    
    varList[[y]] <- v #save raw data in list 
  }
  names(varList) <- shortNames
  return(varList)
}

