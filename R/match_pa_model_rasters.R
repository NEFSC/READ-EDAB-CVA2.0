#' @title Match Species and Model Rasters and Create a Data Frame
#' @description
#' Match species and environmental rasters to build data frame for modeling
#'
#' @param pa_rasters rasterStack representing fisheries presence, absence, and effort with a layer for each timestep. Layers should be named by timestamp, for example, month.year
#' @param model_data  A list containing rasterStacks of each environmental variable to be matched, or the output of \code{norm_env}. Each member of the list should have a name which will correspond to the environmental variable's column name in the resulting dataframe. Each rasterStack should have the same number of layers as \code{pa_rasters}, and each layer should have a name that can be matched to the names in \code{pa_rasters}
#' @param add_static TRUE/FALSE to add static variables such as bathymetry. If TRUE, staticData needs to be supplied.
#' @param static_variables list of rasters containing static variables such as bathymetry. Names of the listed rasters will become the column names.
#'
#' @return \code{merge_pa_model_rasters} returns a data frame with the following columns:
#' \itemize{
#' \item \code{x, y, month, year} indicating longitude, latitude, month, and year
#' \item \code{value} - the column indicating presence, absence, and effort from the fisheries raster
#' \item A number of columns equal to the length of \code{model_data} with the same names containing the corresponding environmental data at each location
#' \item A number of columns equal to the length of \code{static_variables}, if \code{add_static = T}, with column names matching the names of the rasters in the list
#' }
#'

match_pa_model_rasters <- function(pa_rasters, model_data, add_static = TRUE, static_variables){

  ###step 1 - convert rasterstack to DF
  ##rasterStacks are bigger than environmental covariates, so subset to match
  rastSub <- raster::crop(pa_rasters, raster::extent(model_data[[1]][[1]]))
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
  nms <- vector(length = length(model_data))
  for(x in 1:length(model_data)){
    v <- model_data[[x]][[1]] #isolate one environmental variable

    vc <- vector(length = nrow(sppDF)) #create empty vector
    #matching by layer
    for(y in 1:raster::nlayers(v)){
      i <- which(sppDF[,'variable'] == names(v)[y])#match name of raster layer [which correlates to time] to variable (raster name from spp stack)
      vc[i] <- raster::extract(subset(v,y), sppDF[i,c('x','y')], fun = mean) #extract value
      #print(x)
    } #end y
    nms[x] <- names(model_data[[x]]) #name vector with name of rasterStack
    sppDF <- cbind(sppDF, vc)
    print(names(model_data)[x])
  } #end x
  colnames(sppDF)[which(colnames(sppDF) == 'vc')] <- nms

  if(add_static){ #if add_static == True
    nmSD <- vector(length = length(static_variables))
    for(x in 1:length(static_variables)){
      s <- raster::extract(x = static_variables[[x]], y = sppDF[,c('x','y')], fun = 'mean', method = 'simple')
      nmSD[x] <- names(static_variables)[x]
      sppDF <- cbind(sppDF, s)
      print(names(static_variables)[x])
    } #end x
    colnames(sppDF)[which(colnames(sppDF) == 's')] <- nmSD
  } #end if(add_static)

  return(sppDF)

} #end function
