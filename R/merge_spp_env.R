#' @title Merge Species and Environmental Rasters into a Data Frame
#' @description
#' Match species and environmental rasters to build data frame for modeling
#'
#' @param rastStack rasterStack representing fisheries presence, absence, and effort with a layer for each timestep. Layers should be named by timestamp, for example, month.year
#' @param envData output from one of the MOM6 data functions, preferably \code{norm_env}. A list containing rasterStacks of each environmental variable to be matched. Each member of the list should have a name which will correspond to the environmental variable's column name in the resulting dataframe. Each rasterStack should have the same number of layers as \code{rastStack}, and each layer should have a name that can be matched to the names in \code{rastStack}
#' @param addStatic TRUE/FALSE to add static variables such as bathymetry. If TRUE, staticData needs to be supplied.
#' @param staticVars list of rasters containing static variables such as bathymetry. Names of the listed rasters will become the column names.
#'
#' @return \code{merge_spp_env} returns a data frame with the following columns:
#' \itemize{
#' \item \code{x, y, month, year} indicating longitude, latitude, month, and year
#' \item \code{value} - the column indicating presence, absence, and effort from the fisheries raster
#' \item A number of columns equal to the length of \code{envData} with the same names containing the corresponding environmental data at each location
#' \item A number of columns equal to the length of \code{staticData}, if \code{addStatic = T}, with column names matching the names of the rasters in the list
#' }
#'

merge_spp_env <- function(rastStack, envData, addStatic = TRUE, staticVars){

  ###step 1 - convert rasterstack to DF
  ##rasterStacks are bigger than environmental covariates, so subset to match
  rastSub <- raster::crop(rastStack, raster::extent(envData[[1]][[1]]))
  #### convert to data frame for modeling
  sppDF <- as.data.frame(raster::rasterToPoints(rastSub))
  sppDF <- reshape2::melt(sppDF, id.vars = 1:2)

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
    for(y in 1:raster::nlayers(v)){
      i <- which(sppDF[,'variable'] == names(v)[y])#match name of raster layer [which correlates to time] to variable (raster name from spp stack)
      vc[i] <- raster::extract(subset(v,y), sppDF[i,c('x','y')], fun = mean) #extract value
      #print(x)
    } #end y
    nms[x] <- names(envData[[x]]) #name vector with name of rasterStack
    sppDF <- cbind(sppDF, vc)
    print(names(envData)[x])
  } #end x
  colnames(sppDF)[which(colnames(sppDF) == 'vc')] <- nms

  if(addStatic){ #if addStatic == True
    nmSD <- vector(length = length(staticVars))
    for(x in 1:length(staticVars)){
      s <- raster::extract(x = staticVars[[x]], y = sppDF[,c('x','y')], fun = 'mean', method = 'simple')
      nmSD[x] <- names(staticVars)[x]
      sppDF <- cbind(sppDF, s)
      print(names(staticVars)[x])
    } #end x
    colnames(sppDF)[which(colnames(sppDF) == 's')] <- nmSD
  } #end if(addStatic)

  return(sppDF)

} #end function
