merge_spp_env <- function(sppData, name2match, envData, addStatic = TRUE, staticData){
  #sppData is a data.frame of species presence/absence 
  #name2match - name of a column in sppData that matches the names of the layers in envData (all alyers should have the same names)
  #envData is a list of environmental rasterStacks - output from pull_normalize_env
  #addStatic - binary option to add any additional static environmental variables such as bathymetry, rugosity, etc
  #staticData is a list of additional rasters to be added 
  
  for(x in 1:length(envData)){
    v <- envData[[x]] #isolate one environmental variable 
    
    vc <- vector(length = nrow(sppData)) #create empty vector
    #matching by layer
    for(y in 1:dim(v)[3]){
      i <- which(sppData[,'name2match'] == names(v)[y])#match name of raster layer [which correlates to time] to variable (raster name from spp stack) 
      vc[i] <- extract(subset(v,y), sppData[i,c('x','y')], fun = mean) #extract value 
      #print(x)
    } #end y 
    names(vc) <- names(envData)[x] #name vector with name of rasterStack 
    sppData <- cbind(sppData, vc)
    print(names(envData)[x])
  } #end x
  
  if(addStatic){ #if addStatic == True 
    for(x in 1:length(staticData)){
      s <- extract(x = staticData[[x]], y = sppData[,c('x','y')], method = 'simple')
      names(s) <- names(staticData)[x]
      sppData <- cbind(sppData, s)
      print(names(staticData)[x])
    } #end x
    
  } #end if(addStatic)
  
  return(sppData)
  
} #end function