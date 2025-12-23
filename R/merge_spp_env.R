#' @title Merge Fisheries and Environmental Datasets
#' @description
#' Combine fisheries and environmental rasters into a data frame, with optional cleaning functions
#' \itemize{
#' \item \code{merge_spp_env} combines the fisheries and environmental rasters into a data frame
#' \item \code{match_guilds} reduces the number of environmental covariates based on defined feeding and habitat guilds 
#' \item \code{remove_corr} removed correlated covariates
#' \item \code{clean_data} converts the presence/absence/effort (0-2) data to presence/absence data (0-1)
#' }

#' @param rastStack rasterStack representing fisheries presence, absence, and effort with a layer for each timestep. Layers should be named by timestamp, for example, month.year
#' @param envData output from one of the MOM6 data functions, preferably \code{norm_env}. A list containing rasterStacks of each environmental variable to be matched. Each member of the list should have a name which will correspond to the environmental variable's column name in the resulting dataframe. Each rasterStack should have the same number of layers as \code{rastStack}, and each layer should have a name that can be matched to the names in \code{rastStack}
#' @param addStatic TRUE/FALSE to add static variables such as bathymetry. If TRUE, staticData needs to be supplied.
#' @param staticData a file path to a list of rasters containing static variables such as bathymetry. Names of the listed rasters will become the column names. 
#' @param spp_env,se data.frame of fisheries and environmental data
#' @param spp Species name to add to log files and save data to correct directory (see vignette for recommended directory set up)
#' @param spp_col name of column with species name in spp_guild key
#' @param spp_guild file path to csv file containing a list of species and corresponding feeding and habitat guilds. See vignette for recommended set up of this file
#' @param feeding_key,habitat_key file path to csv for feeding and habitat guild keys respectively. See vignette for recommended set up of these files. 
#' @param feeding_col,habitat_col columns names indicating the guild column in the feeding and habitat guild csvs
#' @param static_vars a vector with the column names to be retained
#' @param pa_col column name for presence/absence column
#' @param xy_col a vector with a length of 2 indicating the longitude and latitude column names
#' @param month_col,year_col column names for month and year columns respectively
#' @param name Species name to add to log files and save data to correct directory (see vignette for recommended directory set up)
#' @param skip TRUE/FALSE indicating whether to skip creating the output file if file already exists
#' @param mMin,mMax,yMin,yMax (m)onth and (y)ear ranges to help create names for combined rasterBrick since they aren't always retained and these are necessary for matching fisheries and environmental data

#' @return \code{merge_spp_env} returns a data frame with the following columns: 
#' #' \itemize{
#' \item \code{x, y, month, year} indicating longitude, latitude, month, and year
#' \item \code{value} - the column indicating presence, absence, and effort from the fisheries raster
#' \item A number of columns equal to the length of \code{envData} with the same names containing the corresponding environmental data at each location
#' \item A number of columns equal to the length of \code{staticData}, if \code{addStatic = T}, with column names matching the names of the rasters in the list
#' }
#' @return \code{match_guilds} returns a data frame that contains the static variables listed and only environmental covariates associated with the species' feeding and habitat guilds
#' @return \code{remove_corr} returns a data frame with correlated covariates removed
#' @return \code{clean_data} returns a data frame where presence/absence has been set to 0 for absent and 1 for present


merge_spp_env <- function(rastStack, envData, addStatic = TRUE, staticData){

  ###step 1 - convert rasterstack to DF 
  ##rasterStacks are bigger than environmental covariates, so subset to match 
  rastSub <- raster::crop(rastStack, extent(envData[[1]][[1]]))
  #### convert to data frame for modeling 
  sppDF <- as.data.frame(rasterToPoints(rastSub))
  sppDF <- melt(sppDF, id.vars = 1:2)
  
  #split out month and year
  sppDF$variable <- as.character(sppDF$variable)
  my <- strsplit(x = sppDF$variable, split = '[.]')
  sppDF$month <- unlist(lapply(my, FUN = function(x){x[1]}))
  sppDF$month <- as.numeric(gsub('X', '', sppDF$month))
  sppDF$year <- as.numeric(unlist(lapply(my, FUN = function(x){x[2]})))
  
  #step 2 - match to environmental data
  nms <- vector(length = length(envData))
  for(x in 1:length(envData)){
    v <- envData[[x]][[1]] #isolate one environmental variable 
    
    vc <- vector(length = nrow(sppDF)) #create empty vector
    #matching by layer
    for(y in 1:nlayers(v)){
      i <- which(sppDF[,'variable'] == names(v)[y])#match name of raster layer [which correlates to time] to variable (raster name from spp stack) 
      vc[i] <- extract(subset(v,y), sppDF[i,c('x','y')], fun = mean) #extract value 
      #print(x)
    } #end y 
    nms[x] <- names(envData[[x]]) #name vector with name of rasterStack 
    sppDF <- cbind(sppDF, vc)
    print(names(envData)[x])
  } #end x
  colnames(sppDF)[which(colnames(sppDF) == 'vc')] <- nms
  
  if(addStatic){ #if addStatic == True 
    load(staticData) #staticVars
    nmSD <- vector(length = length(staticVars))
    for(x in 1:length(staticVars)){
      s <- extract(x = staticVars[[x]], y = sppDF[,c('x','y')], fun = 'mean', method = 'simple')
      nmSD[x] <- names(staticVars)[x]
      sppDF <- cbind(sppDF, s)
      print(names(staticVars)[x])
    } #end x
    colnames(sppDF)[which(colnames(sppDF) == 's')] <- nmSD
  } #end if(addStatic)
  
  return(sppDF)
  
} #end function

match_guilds <- function(spp_env, spp, spp_col = 'name', spp_guild, feeding_key, feeding_col = 'Feeding.Guild', habitat_key,  habitat_col = 'Habitat.Guild', static_vars = c('x', 'y', 'month', 'year', 'bathy', 'rugosity', 'dist2coast'), pa_col = 'value'){

  guilds <- read.csv(spp_guild) #load in species list 
  
  g <- which(guilds[,spp_col] %in% spp) #find row associated with target species
  
  #isolate feeding guild
  fGuild <- guilds[g, feeding_col]
  
  #isolate habitat guild
  hGuild <- guilds[g, habitat_col]
  
  ##isolate covariates for each guild type based on keys 
  #Feeding guilds: 
  feeding <- read.csv(feeding_key)
  i <- which(colnames(feeding) == fGuild)
  fInd <- feeding[nzchar(feeding[,i]), i] #get covariates 
  
  
  #Habitat guilds: 
  habitat <- read.csv(habitat_key)
  i <- which(colnames(habitat) == hGuild)
  hInd <- habitat[nzchar(habitat[,i]), i] #get covariates 
  
  ### subset spp_env
  ind <- colnames(spp_env) %in% c(pa_col, static_vars, fInd, hInd)
  se <- spp_env[,ind]
  
  return(se)
  
}

remove_corr <- function(se, pa_col, xy_col, month_col, year_col){
  ind <- which(colnames(se) == pa_col | colnames(se) == xy_col[1] | colnames(se) == xy_col[2] | colnames(se) == month_col | colnames(se) == year_col)
  corInd <- findCorrelation(cor(se[,-ind]), names = T) #find correlated variables 
  if(length(corInd) != 0){
    se <- se[,-which(colnames(se) == corInd)]
  }
  return(se)
}

clean_data <- function(se, pa_col){
  paDF <- se[se[,pa_col] != 0,] #remove unsampled cells (pa_col == 0)
  #switch it back to 0/1 now that we only have grid cells that were sampled 
  paDF[,pa_col] <- replace(paDF[,pa_col], paDF[,pa_col] == 1, 0)
  paDF[,pa_col] <- replace(paDF[,pa_col], paDF[,pa_col] == 2, 1)
  
  paDF2 <- paDF[complete.cases(paDF),]
  
  return(paDF)
}
