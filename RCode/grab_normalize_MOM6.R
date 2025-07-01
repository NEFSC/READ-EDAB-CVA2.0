### function to pull MOM6 environmental data and create a list of normalized data rasters 

########## INPUTS 
varURL <- "https://psl.noaa.gov/cefi_portal/data_index/cefi_data_indexing.Projects.CEFI.regional_mom6.cefi_portal.northwest_atlantic.full_domain.hindcast.json" #url pointing to variable lists for desired run type (hindcast, decadal forecast) and domain (nwa) in question
reqVars <- c('Bottom Temperature') #list of vars to pull 
gt <- 'regrid' #desired grid type (options = regrid or raw)
of <- 'monthly' #desired output frequency (options = daily or monthly)
nms <- paste(month.abb[expand.grid(1:12, 1993:2019)[,1]], expand.grid(1:12, 1993:2019)[,2], sep = '.') #names of layers - will depend on length of time series and frequency of output
north = 45
south = 35
east = -65
west = -78
#urlHead <- "https://psl.noaa.gov/thredds/ncss/grid/Projects/CEFI/regional_mom6/" #first part of the URL to access ncss server to allow for subsetting 
##############

pull_normalize_env <- function(varURL, reqVars, gt, of, nms, north, south, east, west, url_static = "http://psl.noaa.gov/thredds/dodsC/Projects/CEFI/regional_mom6/cefi_portal/northwest_atlantic/full_domain/hindcast/monthly/raw/r20230520/ocean_static.nc"){
  require(jsonlite)
  
  vars <- fromJSON(varURL) #turn json file into a list 
  
  long.name <- url <- nc <- rp <- grid.type <- out.freq <- NULL  #pull the long names, full opendap urls, grid types, and output frequency for indexing which files to pull 
  for(x in 1:length(vars)){
    long.name <- c(long.name, vars[[x]]$cefi_long_name)
    grid.type <- c(grid.type, vars[[x]]$cefi_grid_type)
    out.freq <- c(out.freq, vars[[x]]$cefi_output_frequency)
    #rp <- c(rp, vars[[x]]$cefi_rel_path)
   # nc <- c(nc, vars[[x]]$cefi_filename)
    url <- c(url, vars[[x]]$cefi_opendap)
  }
  
  varList <- avgList <- normList <- vector(mode = 'list', length = length(reqVars)) #initalize empty lists to store all the data 
  
  ##get static vars to define subset
  #get static MOM6 vars to define extent
  ncopendap_static <- nc_open(url_static)
  lon <- ncvar_get(ncopendap_static, "geolon")
  lat <- ncvar_get(ncopendap_static, "geolat")
  e <- extent(min(lon), max(lon), min(lat), max(lat)) #create extent object to define rasterStack extents
  
  se <- extent(west,east, south,north) #define extent to subset to
  
  for(y in 1:length(reqVars)){
    ind <- which(long.name == reqVars[y] & grid.type == gt & out.freq == of) #find appropriate url for the variable 
    
    if(length(ind) == 1){
      u <- nc[ind]
    } else {
      #put error message here about checking inputs 
    }
    
    #load url 
    v <- stack(url[ind]) 
    names(v) <- nms #set names 
    varList[[y]] <- v #save raw data in list 
    
    #subset
    v <- crop(v, se) #this is the rate limiting step
    
    ## create monthly averages 
    avgs <- NULL 
    for(m in 1:12){
      mn <- seq(m, nlayers(v), by = 12) #grab all month xs from timeseries by creating a sequence
      MNS <- subset(v, mn) #subset raster brick
      mm <- mean(MNS) #average 
      avgs <- abind(as.array(mm), avgs, along = 3) #make array and bind together
      #print(x)
    } #end m
    names(avgs) <- month.abb
    avgList[[y]] <- avgs
    
    ##normalize data to monthly averages 
    mth <- rep(1:12, times = length(v)/12) #creates repeating list of 1:12 for each year
    #normalize data
    norm <- NULL 
    for(m in 1:nlayers(vr)){
      #subset both rasterbricks 
      subX <- subset(v,m)
      subA <- subset(avgs,mth[m])
      
      nm <- (subX - subA) / (subX@data@max - subX@data@min)
      norm <- abind(as.array(nm), norm, along = 3)
    } #end m
    names(norm) <- names(v) #set names 
    normList[[y]] <- norm #save
    
    print(reqVars[y])
  } #end y 
  
  names(varList) <- names(avgList) <- names(normList) <- reqVars #set names of lists 
  
  return(list(varList, avgList, normList))

} #end function 

