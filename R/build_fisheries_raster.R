#' @title Convert Standardized Fisheries Data into Presence/Absence/Effort Raster
#' @description
#' Takes output from \code{standardize_data} and converts it into a raster with the spatial and temporal resolution of the provided grid ranging from 0-2 where:
#' \itemize{
#' \item 0 - No fishing effort in the cell
#' \item 1 - Fishing effort present for the target species in the cell
#' \item 2 - Target species caught in the cell
#' }
#'
#' @param data data frame to convert to presence/absence raster
#' @param is_obs TRUE/FALSE indicating whether or not the data is observer or similar fisheries-dependent data. Will force function to check if species is observed at least 30 times throughout timeseries before creating raster
#' @param grid static link to a ncdcf object with the variables lon, lat, time - can be link to remote data - must be able to be read with nc_open
#' @param tm_multiplier multiplier to help convert timestep to POSIX (seconds since origin), defaults to 86400 (number of seconds in a day)
#' @param origin Origin of time series to be used by POSIXct
#' @param all_names a vector containing all possible names for the target species. Must have a length >= 1
#'
#' @return a rasterBrick with the same extent as the provided grid, and a number of layers equal to the timeseries associated with the provided model data
#'

build_fisheries_raster <- function(
  data,
  is_obs = FALSE,
  grid,
  tm_multiplier = 24 * 60 * 60,
  origin = '1993-01-01',
  all_names
) {
  #### step 1 ####
  #get grid to map to
  gridNC <- ncdf4::nc_open(grid)
  #these are vectors of the unique lon/lats on the regridded MOM6 grid
  lonR <- ncdf4::ncvar_get(gridNC, "lon")
  latR <- ncdf4::ncvar_get(gridNC, "lat")
  tm <- as.POSIXct(
    ncdf4::ncvar_get(gridNC, 'time') * tm_multiplier,
    origin = origin
  )
  ncdf4::nc_close(gridNC)

  #### step 2 ####
  #loop through years & months to build rasters

  #first to make this a little bit easier, subset the dataset to have the same year extent as the grid
  tmInd <- which(
    data$year >= min(lubridate::year(tm)) &
      data$year <= max(lubridate::year(tm))
  )
  data <- data[tmInd, ]

  #first, build vector of months/years throughout grid timeseries
  my <- expand.grid(1:12, unique(lubridate::year(tm)))
  my$month.year <- paste(my$Var1, my$Var2, sep = '.')

  #create month.year in data
  data$month.year <- paste(data$month, data$year, sep = '.')

  spRast <- NULL
  dataOK <- TRUE

  #add check for Fisheries dependent data (dataType = Observer); fisheries independent surveys do not need to do this
  dataOK <- TRUE #assume data is good to go
  if (is_obs == TRUE) {
    iSPP <- data$name %in% all_names #was the species caught?

    # m <- unique(data$month[iSPP]) #in which months has the species has been caught?

    #per McHenry et al 2019, fisheries dependent data were only considered if the species was caught at least 30 times across at least 6 different months
    #we lightened that threshold to just 30 observations to account for the fact that some observer programs only go for 5 months of the year
    #if(length(which(iSPP == T)) >= 30 & length(m) >= 6){ #so if these requirements are met
    if (length(which(iSPP == T)) >= 30) {
      print('Observer data meet minimum thresholds...continuing with raster')
    } else {
      print(
        'Observer data do NOT meet minimum thresholds to build raster...make sure that all possible variations of the species name (scientific, common, alternative common names) are included in targetVec'
      )
      dataOK <- FALSE
    }
  } #end if observer

  if (dataOK) {
    for (x in 1:nrow(my)) {
      #create matrix for each month
      mat <- matrix(0, nrow = length(latR), ncol = length(lonR))

      #subset observer data to month and year
      sub <- data[data$month.year == my$month.year[x], ]

      #if there are data
      if (nrow(sub) != 0) {
        ids <- unique(sub$towID)
        for (i in ids) {
          tow <- sub[sub$towID == i, ]

          #find closest lat/lon
          iLon <- DescTools::Closest(x = lonR - 360, a = tow$lon[1], which = T)
          iLat <- DescTools::Closest(x = latR, a = tow$lat[1], which = T)

          if (any(all_names %in% tow$name)) {
            #if species is in tow
            mat[iLat, iLon] <- 2 #replace grid cell with a 2
          } else {
            mat[iLat, iLon] <- 1 #if species is not found but grid cell is sampled, 1
          }
        } #end for i
      } #end if nrow(sub)

      #add to array
      spRast <- abind::abind(spRast, mat, along = 3)
    }

    #### step 3 #####
    #turn into raster brick

    #make into rasterBrick
    sppBrick <- raster::brick(spRast)
    e <- raster::extent(min(lonR), max(lonR), min(latR), max(latR)) #lon/lat from grid
    raster::extent(sppBrick) <- e
    sp::proj4string(sppBrick) <- sp::CRS('+proj=longlat +datum=WGS84 +no_defs')
    sppBrick <- raster::rotate(raster::flip(sppBrick)) #rotate to get orientation right

    #make names for raster layers
    names(sppBrick) <- paste(sprintf("%02d", my$Var1), my$Var2, sep = '.') #use unique month.year combinations
  } else {
    sppBrick <- NULL
  }

  return(sppBrick) #return rasterbrick
} #end function
