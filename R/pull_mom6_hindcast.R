#' @title Pull MOM6 Hindcast Data
#' @description
#' Pull hindcast data from the MOM6 model output based on provided URL from the CEFI portal
#'
#' @param var_url URL pointing to JSON table variable lists for desired MOM6 hindcast and domain
#' @param req_var name of variable to pull. Must match names in the 'cefi_long_name' column provided JSON table
#' @param gt desired grid type. Must match one of the options in the 'cefi_grid_type' column in provided JSON table
#' @param of desired output frequency. Must match one of the options in the 'cefi_output_frequency' column in provided JSON table
#' @param bounds xmin, xmax, ymin, ymax of desired output raster
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
    release
) {
  
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
  
  if(length(ind) > 1){ #if ind matches multiple files (which is the case for MLD because the names aren't unique)
    #max/min MLD are provided on regridded products, find where those are and remove them. 
    iMin <- grep('min', url[ind])
    iMax <- grep('max', url[ind])
    ind <- ind[-c(iMin, iMax)]
  }
  
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
  
  # Define a chunk size (number of time steps to pull per request).
  # If you still get the DATADDS error, lower this number (e.g., 12 or 24).
  chunk_size <- 50 
  var <- NULL
  # Loop through time using chunks
  for (start_t in seq(1, length(tm), by = chunk_size)) {
    
    # Calculate how many time steps to pull in this specific chunk
    # (Prevents overshooting the end of the time series)
    count_t <- min(chunk_size, length(tm) - start_t + 1)
    
    # Pull the chunk
    v_chunk <- ncdf4::ncvar_get(
      v,
      names(v$var),
      start = c(lonInd[1], latInd[1], start_t),
      count = c(length(lonInd), length(latInd), count_t)
    )
    
    # Bind the chunk along the 3rd dimension (Time)
    var <- abind::abind(var, v_chunk, along = 3)
  }
  
  # Close the NetCDF connection once finished
  ncdf4::nc_close(v)
  
  # Convert the array to a SpatRaster
  # Because ncdf4 loads arrays as [Lon, Lat, Time], we transpose it to [Lat, Lon, Time] 
  # so terra reads the rows and columns correctly.
  r_list <- lapply(1:dim(varArr)[3], function(i) {
    terra::rast(t(varArr[,,i]))
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
