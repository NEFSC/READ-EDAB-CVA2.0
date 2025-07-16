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
      s <- extract(x = staticVars[[x]], y = sppDF[,c('x','y')], method = 'simple')
      nmSD[x] <- names(staticVars)[x]
      sppDF <- cbind(sppDF, s)
      print(names(staticVars)[x])
    } #end x
    colnames(sppDF)[which(colnames(sppDF) == 's')] <- nmSD
  } #end if(addStatic)
  
  return(sppDF)
  
} #end function