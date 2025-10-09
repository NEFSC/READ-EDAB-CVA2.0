#' @title Merge Fisheries and Environmental Datasets
#' @description
#' Combine fisheries and environmental rasters into a data frame
#' \itemize{
#' \item \code{merge_spp_env} combines the fisheries and environmental rasters into a data frame
#' \item \code{makeDF} is the wrapper function that produces the logs and has skip functionality
#' }
#' @param rastStack rasterStack representing fisheries presence, absence, and effort with a layer for each timestep. Layers should be named by timestamp, for example, month.year
#' @param envData output from one of the MOM6 data functions, preferably \code{norm_env}. A list containing rasterStacks of each environmental variable to be matched. Each member of the list should have a name which will correspond to the environmental variable's column name in the resulting dataframe. Each rasterStack should have the same number of layers as \code{rastStack}, and each layer should have a name that can be matched to the names in \code{rastStack}
#' @param addStatic TRUE/FALSE to add static variables such as bathymetry. If TRUE, staticData needs to be supplied.
#' @param staticData a file path to a list of rasters containing static variables such as bathymetry. Names of the listed rasters will become the column names. 
#' @param name Species name to add to log files and save data to correct directory (see vignette for recommended directory set up)
#' @param skip TRUE/FALSE indicating whether to skip creating the raster file if file already exists
#' @param mMin,mMax,yMin,yMax (m)onth and (y)ear ranges to help create names for combined rasterBrick since they aren't always retained and these are necessary for matching fisheries and environmental data
#' @return \code{merge_spp_env} returns a data frame with the following columns: 
#' #' \itemize{
#' \item \code{x, y, month, year} indicating longitude, latitude, month, and year
#' \item \code{value} - the column indicating presence, absence, and effort from the fisheries raster
#' \item A number of columns equal to the length of \code{envData} with the same names containing the corresponding environmental data at each location
#' \item A number of columns equal to the length of \code{staticData}, if \code{addStatic = T}, with column names matching the names of the rasters in the list
#' }
#' @return \code{makeDF} does not return anything to the environment; it only saves the final dataset and writes to the log file

merge_spp_env <- function(rastStack, envData, addStatic = TRUE, staticData){
  #rastStack is the output of create_rast - a raster stack with a layer for each time step
  
  #envData is a list of environmental rasterStacks - output from pull_normalize_env
  
  #addStatic - binary option to add any additional static environmental variables such as bathymetry, rugosity, etc
  
  #staticData is a list of additional rasters to be added 
  
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

makeDF <- function(name, skip, mMin, mMax, yMin, yMax){
  sink(file.path(getwd(), 'logs', 'build_data_frames.log'), append = T)
  # Ensure the sinks are closed when the function exits, regardless of how it exits.
  on.exit({
    #sink(type = "message")
    sink()
  })
  
  print(Sys.time())
  print(name)
  
  if(skip){
    if(file.exists(paste(file.path(getwd(),name), 'pa_clean.RData', sep = '/'))){
      print('file exists and skip == T, so skipping this species!')
      return(NA)
    } else {
      combinedRasts <- brick(paste(file.path(getwd(),name, 'input_rasters'), 'combined_rasters.nc', sep = '/')) #combinedRasts
      
      my <- expand.grid(mMin:mMax, yMin:yMax)
      names(combinedRasts) <- paste(sprintf("%02d", my$Var1), my$Var2, sep = '.')
      
      print(paste0(name, '- merging spp and env data -', Sys.time()))
      df <- merge_spp_env(rastStack = combinedRasts, envData = norm, addStatic = TRUE, staticData = './Data/staticVariables_cropped_normZ.RData')
      save(df, file = paste(file.path(getwd(),name), 'presence_absence_environment.RData', sep = '/'))
      
      print(paste0(name, '- matching guild -', Sys.time()))
      
      dfG <- match_guilds(spp_env = df, spp = c('Summer Flounder', 'Fluke', "PARALICHTHYS DENTATUS"), spp_col = 'SCI_NAME', spp_guild = 'spp_list.csv', feeding_key = 'feeding_guilds.csv', feeding_col = 'Feeding.Guild', habitat_key = 'habitat_guilds.csv',  habitat_col = 'Habitat.Guild', static_vars = c('x', 'y', 'month', 'year', 'bathy', 'rugosity', 'dist2coast'), pa_col = 'value')
      save(dfG, file = paste(file.path(getwd(),name), 'pa_guild.RData', sep = '/'))
      
      print(paste0(name, '- converting to binary -', Sys.time()))
      dfB <- clean_data(dfG, pa_col = 'value')
      save(dfB, file = paste(file.path(getwd(),name), 'pa_binary.RData', sep = '/'))
      
      print(paste0(name, ' - removing correlated variables -', Sys.time()))
      dfC <- remove_corr(dfB, pa_col = 'value', xy_col = c('x', 'y'), month_col = 'month', year_col = 'year')
      save(dfC, file = paste(file.path(getwd(),name), 'pa_clean.RData', sep = '/'))
    }
  } else {
    combinedRasts <- brick(paste(file.path(getwd(),name, 'input_rasters'), 'combined_rasters.nc', sep = '/')) #combinedRasts
    
    my <- expand.grid(mMin:mMax, yMin:yMax)
    names(combinedRasts) <- paste(sprintf("%02d", my$Var1), my$Var2, sep = '.')
    
    print(paste0(name, '- merging spp and env data -', Sys.time()))
    df <- merge_spp_env(rastStack = combinedRasts, envData = norm, addStatic = TRUE, staticData = './Data/staticVariables_cropped_normZ.RData')
    save(df, file = paste(file.path(getwd(),name), 'presence_absence_environment.RData', sep = '/'))
    
    print(paste0(name, '- matching guild -', Sys.time()))
    
    dfG <- match_guilds(spp_env = df, spp = c('Summer Flounder', 'Fluke', "PARALICHTHYS DENTATUS"), spp_col = 'SCI_NAME', spp_guild = 'spp_list.csv', feeding_key = 'feeding_guilds.csv', feeding_col = 'Feeding.Guild', habitat_key = 'habitat_guilds.csv',  habitat_col = 'Habitat.Guild', static_vars = c('x', 'y', 'month', 'year', 'bathy', 'rugosity', 'dist2coast'), pa_col = 'value')
    save(dfG, file = paste(file.path(getwd(),name), 'pa_guild.RData', sep = '/'))
    
    print(paste0(name, '- converting to binary -', Sys.time()))
    dfB <- clean_data(dfG, pa_col = 'value')
    save(dfB, file = paste(file.path(getwd(),name), 'pa_binary.RData', sep = '/'))
    
    print(paste0(name, ' - removing correlated variables -', Sys.time()))
    dfC <- remove_corr(dfB, pa_col = 'value', xy_col = c('x', 'y'), month_col = 'month', year_col = 'year')
    save(dfC, file = paste(file.path(getwd(),name), 'pa_clean.RData', sep = '/'))
  }
}