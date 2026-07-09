#' @title Convert Standardized Fisheries Data into Presence/Absence/Effort Data Frame
#' @description
#' Takes output from \code{standardize_data} and matches it to the provided grid ranging from 0-2 where:
#' \itemize{
#' \item 0 - Fishing effort present for the target species in the cell
#' \item 1 - Target species caught in the cell
#' }
#'
#' @param data a data.frame containing the data to use from a specific source 
#' @param is_obs a TRUE/FALSE indicating whether or not the data source is observer or similar fisheries-dependent data. Will force function to check if species is observed at least 30 times throughout timeseries before creating using data
#' @param grid static link to a ncdcf object with the variables lon, lat, time - can be link to remote data - must be able to be read with nc_open
#' @param tm_multiplier multiplier to help convert timestep to POSIX (seconds since origin), defaults to 86400 (number of seconds in a day)
#' @param origin Origin of time series to be used by POSIXct
#' @param all_names a vector containing all possible names for the target species. Must have a length >= 1
#'
#' @return a data.frame with the same extent as the provided grid, and a number of layers equal to the timeseries associated with the provided model data
#'

build_fisheries_df <- function(data, is_obs, grid, tm_multiplier = 24 * 60 * 60, origin = '1993-01-01', all_names){
  
  #### step 1 ####
  #get grid to map to
  gridNC <- ncdf4::nc_open(grid)
  #these are vectors of the unique lon/lats on the regridded MOM6 grid
  lonR <- ncdf4::ncvar_get(gridNC, "lon")
  latR <- ncdf4::ncvar_get(gridNC, "lat")
  tm <- as.POSIXct(ncdf4::ncvar_get(gridNC, 'time') * tm_multiplier, origin = origin)
  ncdf4::nc_close(gridNC)
  
  #### step 2 ####
  depCheck <- TRUE
  #check for Fisheries dependent data (dataType = Observer); fisheries independent surveys do not need to do this
  if(is_obs == TRUE){
    iSPP <- data$name %in% all_names #was the species caught?

    #per McHenry et al 2019, fisheries dependent data were only considered if the species was caught at least 30 times across at least 6 different months
    #we lightened that threshold to just 30 observations to account for the fact that some observer programs only go for 5 months of the year
    if(length(which(iSPP == T)) <= 30){ #so if these requirements are not met, change flag & print warning
      print('Observer data do NOT meet minimum thresholds to include in data.frame...make sure that all possible variations of the species name (scientific, common, alternative common names) are included in targetVec')
      depCheck <- FALSE 
    }
  } #end if observer
   
  
  #### step 3 ####
  #loop through years & months to build presence-absence data frame
    if(depCheck){
    # remove any zero/NA counts - accounts for differences in sources providing presence/absence vs count data
    data <- data[!is.na(data$count) & data$count != 0,]
    
    #build vector of months/years throughout grid timeseries
    my <- expand.grid(1:12, unique(lubridate::year(tm)))
    my$month.year <- paste(my$Var1, my$Var2, sep = '.')
    
    #create month.year in data
    data$month.year <- paste(data$month, data$year, sep = '.')
  
    DF <- NULL
  
    for(x in 1:nrow(my)){ #for each month.year
      #subset data to month and year - we do this so we can merge data within the same grid cells by month.year
      sub <- data[data$month.year == my$month.year[x], ]
    
      #if there are data
      if(nrow(sub) != 0){
        ids <- unique(sub$towID)
        paM <- NULL #create null object for all data from one month
        for(i in ids){
          tow <- sub[sub$towID == i, ]
          
          #find closest grid cell lat/lon
          iLon <- DescTools::Closest(x = lonR, a = tow$lon[1], which = F)
          iLat <- DescTools::Closest(x = latR, a = tow$lat[1], which = F)
          
          #pull necessary information from tow, add corresponding grid lon/lat
          pa <- data.frame(year = tow$year[1], month = tow$month[1], tow.lon = tow$lon[1], tow.lat = tow$lat[1], grid.lon = iLon, grid.lat = iLat)
          
          #generate grid cell lat/lon ID to help combine data later 
          idLon <- DescTools::Closest(x = lonR-360, a = tow$lon[1], which = T)
          idLat <- DescTools::Closest(x = latR, a = tow$lat[1], which = T)
          pa$gridID <- paste(idLon, idLat, sep = '-') #add to pa data.frame
          
          #add presence/absence column - default to 0 (absent, but effort)
          pa$pa <- 0
          
          if(any(all_names %in% tow$name)){ #if species is in tow
            pa$pa <- 1 #replace with 1 
          }
          
          paM <- rbind(paM, pa) #combine
      } #end for i
        #remove bad gridIDs
        if('NA-NA' %in% paM$gridID){
          paM <- paM[-which(paM$gridID == 'NA-NA'),]
        }
        
        if(nrow(paM) != 0){
          paMax <- stats::aggregate(paM, by = list(paM$gridID), FUN = max) #combine by grid cell ID to get max presence/absence within each grid cell
          DF <- rbind(DF, paMax[,-1])
        }
       # print(my$month.year[x])
      } #end if nrow(sub)
      
    }  #end x
    
    DF <- DF[order(DF$year, DF$month),]
    
    return(DF) #return data.frame
  } #end if depCheck
  
} #end function

