merge_spp_env <- function(rastStack, envData, addStatic = TRUE, staticData){
  #rastStack is the output of create_rast - a raster stack with a layer for each time step
  
  #envData is a list of environmental rasterStacks - output from pull_normalize_env
  
  #addStatic - binary option to add any additional static environmental variables such as bathymetry, rugosity, etc
  
  #staticData is a list of additional rasters to be added 
  
  ###step 1 - convert rasterstack to DF 
  #### convert to data frame for modeling 
  sppDF <- as.data.frame(rasterToPoints(rastStack))
  sppDF <- melt(sppDF, id.vars = 1:2)
  
  #split out month and year
  sppDF$variable <- as.character(sppDF$variable)
  my <- strsplit(x = sppDF$variable, split = '[.]')
  sppDF$month <- unlist(lapply(my, FUN = function(x){x[1]}))
  sppDF$year <- as.numeric(unlist(lapply(my, FUN = function(x){x[2]})))
  
  #step 2 - match to environmental data
  for(x in 1:length(envData)){
    v <- envData[[x]] #isolate one environmental variable 
    
    vc <- vector(length = nrow(sppDF)) #create empty vector
    #matching by layer
    for(y in 1:dim(v)[3]){
      i <- which(sppDF[,'variable'] == names(v)[y])#match name of raster layer [which correlates to time] to variable (raster name from spp stack) 
      vc[i] <- extract(subset(v,y), sppDF[i,c('x','y')], fun = mean) #extract value 
      #print(x)
    } #end y 
    names(vc) <- names(envData)[x] #name vector with name of rasterStack 
    sppDF <- cbind(sppDF, vc)
    print(names(envData)[x])
  } #end x
  
  if(addStatic){ #if addStatic == True 
    load(staticData)
    for(x in 1:length(staticData)){
      s <- extract(x = staticData[[x]], y = sppDF[,c('x','y')], method = 'simple')
      names(s) <- names(staticData)[x]
      sppDF <- cbind(sppDF, s)
      print(names(staticData)[x])
    } #end x
    
  } #end if(addStatic)
  
  return(sppDF)
  
} #end function