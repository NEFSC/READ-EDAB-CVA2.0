#' @title Pull MOM6 Hindcast Data
#' @description
#' Pull hindcast data from the MOM6 model output based on provided URL from the CEFI portal
#'
#' @param varURL URL pointing to JSON table variable lists for desired MOM6 hindcast and domain
#' @param reqVars vector of variable names to pull. Must match names in the 'cefi_long_name' column provided JSON table
#' @param shortNames vector of simplified variable names to help name resulting raster files. Must be the same length as reqVars.
#' @param gt desired grid type. Must match one of the options in the 'cefi_grid_type' column in provided JSON table
#' @param of desired output frequency. Must match one of the options in the 'cefi_output_frequency' column in provided JSON table
#' @param bounds xmin, xmax, ymin, ymax of desired output raster
#' @param static URL to static grid for MOM6
#' @param release release code. Must match one of the options in the 'cefi_release' column in provided JSON table
#'
#' @return a list whose length is equal to the number of variables supplied, where each item in the list is a rasterStack of data associated with that variable

pull_hind <- function(varURL, reqVars, shortNames, gt = 'regrid', of = 'monthly', bounds = c(-78,-65, 35,45), static, release){

  vars <- jsonlite::fromJSON(varURL) #turn json file into a list

  long.name <- url <- grid.type <- out.freq <- rl <- NULL  #pull the long names, full opendap urls, grid types, and output frequency for indexing which files to pull
  for(x in 1:length(vars)){
    long.name <- c(long.name, vars[[x]]$cefi_long_name)
    grid.type <- c(grid.type, vars[[x]]$cefi_grid_type)
    out.freq <- c(out.freq, vars[[x]]$cefi_output_frequency)
    url <- c(url, vars[[x]]$cefi_opendap)
    rl <- c(rl, vars[[x]]$cefi_release)
  }

  rawList <- vector(mode = 'list', length = length(reqVars)) #initalize empty lists to store all the data

  #get info for subsetting
  #putting subsetting back because everything else takes too long otherwise
  stat <- ncdf4::nc_open(static)
  lon <-ncdf4::ncvar_get(stat, "geolon")
  lat <- ncdf4::ncvar_get(stat, "geolat")
  ncdf4::nc_close(stat)

  e <- raster::extent(min(lon), max(lon), min(lat), max(lat)) #extent
  se <- raster::extent(bounds) #extent to subset to

  for(y in 1:length(reqVars)){
    ind <- which(long.name == reqVars[y] & grid.type == gt & out.freq == of  & rl == release) #find appropriate url for the variable

    #load url
    v <- raster::stack(url[ind])

    #create and set names
    n <- matrix(unlist(strsplit(names(v), split = '[.]')), ncol =3, nrow = raster::nlayers(v), byrow = T)
    n[,1] <- gsub('X', replacement = '', n[,1])

    names(v) <- paste(n[,2], n[,1], sep = '.') #set names
    raster::extent(v) <- e #set extent
    #subset
    v <- raster::crop(v, se) #this is the rate limiting step

    rawList[[y]] <- v #save raw data in list
  }
  names(rawList) <- shortNames
  return(rawList)
}

