#' @title Pull MOM6 Hindcast Data
#' @description
#' Pull hindcast data from the MOM6 model output based on provided URL from the CEFI portal
#'
#' @param var_url URL pointing to JSON table variable lists for desired MOM6 hindcast and domain
#' @param req_var name of variable to pull. Must match names in the 'cefi_long_name' column provided JSON table
#' @param gt desired grid type. Must match one of the options in the 'cefi_grid_type' column in provided JSON table
#' @param of desired output frequency. Must match one of the options in the 'cefi_output_frequency' column in provided JSON table
#' @param bounds xmin, xmax, ymin, ymax of desired output raster
#' @param static_grid URL to static grid for MOM6
#' @param release release code. Must match one of the options in the 'cefi_release' column in provided JSON table
#'
#' @return  a spatRaster of data associated with the requested variable
#'
#'@export

pull_mom6_hindcast <- function(
  var_url,
  req_var,
  gt = 'regrid',
  of = 'monthly',
  bounds = c(-78, -65, 35, 45),
  static_grid,
  release
) {

  #get grid info for subsetting to help reduce computation time
  stat <- ncdf4::nc_open(static)
  lon <- ncdf4::ncvar_get(stat, "geolon")
  lat <- ncdf4::ncvar_get(stat, "geolat")
  ncdf4::nc_close(stat)

  e <- terra::ext(min(lon), max(lon), min(lat), max(lat)) #define grid extent
  se <- terra::ext(bounds) #define extent to subset to
  
  vars <- jsonlite::fromJSON(var_url) #turn json file into a list
  
  #pull the long names, full opendap urls, grid types, and output frequency for indexing which files to pull
  long.name <- url <- grid.type <- out.freq <- rl <- NULL 
  for (x in 1:length(vars)) {
    long.name <- c(long.name, vars[[x]]$cefi_long_name)
    grid.type <- c(grid.type, vars[[x]]$cefi_grid_type)
    out.freq <- c(out.freq, vars[[x]]$cefi_output_frequency)
    url <- c(url, vars[[x]]$cefi_opendap)
    rl <- c(rl, vars[[x]]$cefi_release)
  }
  
  #find appropriate url for requested variable
    ind <- which(
      long.name == req_vars &
        grid.type == gt &
        out.freq == of &
        rl == release
    )

    #load url
    #v <- raster::stack(url[ind])
    v <- terra::rast(url[ind])

    #create and set names using month and year 
    n <- matrix(
      unlist(strsplit(names(v), split = '[.]')),
      ncol = 3,
      nrow = raster::nlayers(v),
      byrow = T
    )
    n[, 1] <- gsub('X', replacement = '', n[, 1])

    names(v) <- paste(n[, 2], n[, 1], sep = '.') #set names
    terra::ext(v) <- e #set extent
    #subset
    v <- terra::crop(v, se) #this is the rate limiting step

  return(v)
}
