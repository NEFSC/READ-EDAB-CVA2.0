#' @title Pull MOM6 Decadal Forecast Data
#' @description
#' Pull decadal forecast data from the MOM6 model output based on provided URL from the CEFI portal. Decadal forecasts currently have 10 ensemble members, which have slight variations in the forcing mechanisms and are not forced by climate projections. As such ensemble members are currently averaged. 
#'
#' @param var_url URL pointing to JSON table variable lists for desired MOM6 forecast and domain
#' @param req_var variable name to pull. Must match names in the 'cefi_long_name' column provided JSON table
#' @param gt desired grid type. Must match one of the options in the 'cefi_grid_type' column in provided JSON table
#' @param of desired output frequency. Must match one of the options in the 'cefi_output_frequency' column in provided JSON table
#' @param bounds xmin, xmax, ymin, ymax of desired output raster
#' @param static URL to static grid for MOM6
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
  static,
  release,
  init,
  ens
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

  #get info for subsetting
  #putting subsetting back because everything else takes too long otherwise
  stat <- ncdf4::nc_open(static)
  lon <- ncdf4::ncvar_get(stat, "geolon")
  lat <- ncdf4::ncvar_get(stat, "geolat")
  ncdf4::nc_close(stat)

  e <- terra::::ext(min(lon), max(lon), min(lat), max(lat)) #extent
  se <- terra::ext(bounds) #extent to subset to

    ind <- which(
      long.name == req_var &
        grid.type == gt &
        out.freq == of &
        rl == release &
        init.date == init
    ) #find appropriate url for the variable

    #load url with netcdf to account for ensemble members
    var <- NULL
    r <- ncdf4::nc_open(url[ind])
    tm <- ncdf4::ncvar_get(r, 'lead')
    for (m in 1:10) {
      vm <- NULL
      for (z in 1:length(tm)) {
        v <- ncdf4::ncvar_get(
          r,
          names(r$var),
          start = c(1, 1, z, m),
          count = c(-1, -1, 1, 1)
        )
        vm <- abind::abind(vm, v, along = 3)
      } #end z
      var <- abind::abind(var, vm, along = 4)
    } #end m
    ncdf4::nc_close(r)

    ##take average of ensemble members
    varAvg <- apply(var, MARGIN = c(1:3), FUN = mean, na.rm = T)

    #flip to get orientation right
    varFlip <- aperm(varAvg, c(2, 1, 3))
    varFlip <- varFlip[nrow(varFlip):1, , ]

    #convert to raster
    v <- terra::rast(varFlip)
    terra::ext(v) <- e
    terra::crs(v) <- "+proj=longlat +datum=WGS84 +no_defs"

    #create and set names
    yr <- as.numeric(substr(init, 2, 5))
    yr10 <- yr + 9
    nms <- expand.grid(1:12, yr:yr10)

    names(v) <- paste(nms[, 1], nms[, 2], sep = '.') #set names
    terra::ext(v) <- e #set extent
    #subset
    v <- terra::crop(v, se) #this is the rate limiting step
  
  return(v)
}
