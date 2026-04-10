#' @title Pull MOM6 Forecast Data
#' @description
#' Pull forecast data from the MOM6 model output based on provided URL from the CEFI portal.
#'
#' @param var_url URL pointing to JSON table variable lists for desired MOM6 forecast and domain
#' @param req_vars vector of variable names to pull. Must match names in the 'cefi_long_name' column provided JSON table
#' @param short_names vector of simplified variable names to help name resulting raster files. Must be the same length as req_vars.
#' @param gt desired grid type. Must match one of the options in the 'cefi_grid_type' column in provided JSON table
#' @param of desired output frequency. Must match one of the options in the 'cefi_output_frequency' column in provided JSON table
#' @param bounds xmin, xmax, ymin, ymax of desired output raster
#' @param static URL to static grid for MOM6
#' @param release release code. Must match one of the options in the 'cefi_release' column in provided JSON table
#' @param init initialization code. Must match one of the options in the 'cefi_init_date' column in provided JSON table. For forecast only
#' @param ens ensemble member. Must be equal to 1-10. For decadal forecasts, different ensemble members represent slightly different forcing scenarios. For forecast only.
#'
#' @return a list whose length is equal to the number of variables supplied, where each item in the list is a rasterStack of data associated with that variable

pull_mom6_forecast <- function(
  var_url,
  req_vars,
  short_names,
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

  rawList <- vector(mode = 'list', length = length(req_vars)) #initalize empty lists to store all the data

  #get info for subsetting
  #putting subsetting back because everything else takes too long otherwise
  stat <- ncdf4::nc_open(static)
  lon <- ncdf4::ncvar_get(stat, "geolon")
  lat <- ncdf4::ncvar_get(stat, "geolat")
  ncdf4::nc_close(stat)

  e <- raster::extent(min(lon), max(lon), min(lat), max(lat)) #extent
  se <- raster::extent(bounds) #extent to subset to

  for (y in 1:length(req_vars)) {
    ind <- which(
      long.name == req_vars[y] &
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
    v <- raster::brick(varFlip)
    raster::extent(v) <- e
    raster::crs(v) <- "+proj=longlat +datum=WGS84 +no_defs"

    #create and set names
    yr <- as.numeric(substr(init, 2, 5))
    yr10 <- yr + 9
    nms <- expand.grid(1:12, yr:yr10)

    names(v) <- paste(nms[, 1], nms[, 2], sep = '.') #set names
    raster::extent(v) <- e #set extent
    #subset
    v <- raster::crop(v, se) #this is the rate limiting step

    rawList[[y]] <- v #save raw data in list
  }
  names(rawList) <- short_names
  return(rawList)
}
