#' @title Pull MOM6 Decadal Forecast Data
#' @description
#' Pull decadal forecast data from the MOM6 model output based on provided URL from the CEFI portal. Decadal forecasts currently have 10 ensemble members, which have slight variations in the forcing mechanisms and are not forced by climate projections. As such ensemble members are currently averaged. 
#'
#' @param var_url URL pointing to JSON table variable lists for desired MOM6 forecast and domain
#' @param req_var variable name to pull. Must match names in the 'cefi_long_name' column provided JSON table
#' @param gt desired grid type. Must match one of the options in the 'cefi_grid_type' column in provided JSON table
#' @param of desired output frequency. Must match one of the options in the 'cefi_output_frequency' column in provided JSON table
#' @param bounds xmin, xmax, ymin, ymax of desired output raster
#' @param release release code. Must match one of the options in the 'cefi_release' column in provided JSON table
#' @param init initialization code. Must match one of the options in the 'cefi_init_date' column in provided JSON table. For forecast only
#'
#' @return a spatRaster of data associated with the requested variable, averaged across the 10 ensemble members
#'
#'@export

pull_mom6_forecast <- function(
    var_url,
    req_var,
    gt = 'regrid',
    of = 'monthly',
    bounds = c(-78, -65, 35, 45),
    release,
    init
) {
  vars <- jsonlite::fromJSON(var_url) #turn json file into a list
  
  long.name <- url <- grid.type <- out.freq <- rl <- init.date <- NULL #pull the long names, full opendap urls, grid types, and output frequency for indexing which files to pull
  for (x in 1:length(vars)) {
    long.name <- c(long.name, vars[[x]]$cefi_long_name)
    grid.type <- c(grid.type, vars[[x]]$cefi_grid_type)
    out.freq <- c(out.freq, vars[[x]]$cefi_output_frequency)
    url <- c(url, vars[[x]]$cefi_opendap)
    rl <- c(rl, vars[[x]]$cefi_release)
    init.date <- c(init.date, vars[[x]]$cefi_init_date)
  }
  
  ind <- which(
    long.name == req_var &
      grid.type == gt &
      out.freq == of &
      rl == release &
      init.date == init
  ) #find appropriate url for the variable
  
  if(length(ind) > 1){ #if ind matches multiple files (which is the case for MLD because the names aren't unique)
    #max/min MLD are provided on regridded products, find where those are and remove them. 
    iMin <- grep('min', url[ind])
    iMax <- grep('max', url[ind])
    ind <- ind[-c(iMin, iMax)]
  }
  
  #load url with netcdf to account for ensemble members
  var <- NULL
  r <- ncdf4::nc_open(url[ind])
  #get lon/lat first for subsetting
  lon <- ncdf4::ncvar_get(r, "lon")
  lat <- ncdf4::ncvar_get(r, "lat")
  
  #find indexes for lon/lat to crop to bounding box
  lonInd <- which(lon >= bounds[1] & lon <= bounds[2])
  latInd <- which(lat >= bounds[3] & lat <= bounds[4])
  
  tm <- ncdf4::ncvar_get(r, 'lead')
  for (m in 1:10) {
    vm <- NULL
    for (z in 1:length(tm)) {
      v <- ncdf4::ncvar_get(
        r,
        names(r$var),
        start = c(lonInd[1], latInd[1], z, m),
        count = c(length(lonInd), length(latInd), 1, 1)
      )
      vm <- abind::abind(vm, v, along = 3)
    } #end z
    var <- abind::abind(var, vm, along = 4)
  } #end m
  ncdf4::nc_close(r)
  
  ##take average of ensemble members
  varAvg <- apply(var, MARGIN = c(1:3), FUN = mean, na.rm = T)
  
  # Convert the array to a SpatRaster
  # Because ncdf4 loads arrays as [Lon, Lat, Time], we transpose it to [Lat, Lon, Time] 
  # so terra reads the rows and columns correctly.
  r_list <- lapply(1:dim(varAvg)[3], function(i) {
    terra::rast(t(varAvg[,,i]))
  })
  cropped_rast <- terra::rast(r_list)
  
  # Apply the correct spatial metadata
  terra::ext(cropped_rast) <- c(min(lon[lonInd]), max(lon[lonInd]), min(lat[latInd]), max(lat[latInd]))
  terra::crs(cropped_rast) <- "EPSG:4326" # Or whatever coordinate system the data uses
  
  #create and set names
  yrInit <- as.numeric(substr(init, 2, 5))
  d <- as.POSIXct(tm * 60 * 60 * 24, origin = paste(yrInit, '01', '01', sep = '-'))
  nms <- cbind(lubridate::month(d), lubridate::year(d))
  names(cropped_rast) <- paste(nms[, 1], nms[, 2], sep = '.') #set names
  
  return(cropped_rast)
}
