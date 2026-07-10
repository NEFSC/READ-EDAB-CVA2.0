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
  stat <- ncdf4::nc_open(static_grid)
  
  ncdf4::nc_close(stat)

  #e <- terra::ext(min(lon), max(lon), min(lat), max(lat)) #define grid extent
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
      long.name == req_var &
        grid.type == gt &
        out.freq == of &
        rl == release
    )

    #load url
    #v <- raster::stack(url[ind])
    v <- ncdf4::nc_open(url[ind])
    
    #get dimensions
    lon <- ncdf4::ncvar_get(v, "lon")
    lat <- ncdf4::ncvar_get(v, "lat")
    tm <- as.POSIXct(ncdf4::ncvar_get(v, 'time')*60*60*24, origin = '1993-01-01')
    
    #find indexes for lon/lat to crop to bounding box
    lonInd <- which(lon >= bounds[1] & lon <= bounds[2])
    latInd <- which(lat >= bounds[3] & lat <= bounds[4])
    
    #pull variable
    var <- ncdf4::ncvar_get(v, 
                            names(v$var),
                            start = c(lonInd[1], latInd[1], 1),
                            count = c(length(lonInd), length(latInd), -1))
    ncdf4::nc_close(v)
    
    # Convert the array to a SpatRaster
    # Because ncdf4 loads arrays as [Lon, Lat, Time], we transpose it to [Lat, Lon, Time] 
    # so terra reads the rows and columns correctly.
    r_list <- lapply(1:dim(var)[3], function(i) {
      terra::rast(t(var[,,i]))
    })
    cropped_rast <- terra::rast(r_list)
    
    # Apply the correct spatial metadata
    terra::ext(cropped_rast) <- c(min(lon[lonInd]), max(lon[lonInd]), min(lat[latInd]), max(lat[latInd]))
    terra::crs(cropped_rast) <- "EPSG:4326" # Or whatever coordinate system the data uses

    #create and set names using month and year 
    m <- lubridate::month(tm)
    yr <- lubridate::year(tm)

    names(cropped_rast) <- paste(m, yr, sep = '.') #set names
    #terra::ext(v) <- e #set extent
    
    #flip it
    cropped_rast <- terra::flip(cropped_rast, direction="vertical")
   
  return(cropped_rast)
}
