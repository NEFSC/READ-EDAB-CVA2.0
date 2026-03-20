#' @title Standardize Fisheries Dependent and Independent Datasets
#' @description
#' Pulls and standardizes fisheries independent and dependent datasets from the NEFSC survey and observer programs using ROracle. The function also can open a local CSV file and standardize the output.
#'
#' @param dataType type of data to be standardize. Must be one of the following: 'NESurveys', 'NEObserver', or 'CSV'
#' @param channel connection to remote databases. Only required for 'NESurveys' or 'NEObserver'
#' @param csv path to local CSV file
#' @param csvCols Column names in csv file in the following order: 'towid', 'longitude', 'latitude', 'date', 'count (can be count/abundance/density, etc)', 'name'. Date must be in a format that can be converted to POSIX with as.POSIXct
#' @param yrRange A vector with length of 2 indicating the start and end, inclusive, year of the desired time series
#'
#' @return a data frame. It is recommended to save this as a csv file as pulling survey and observer datasets does take time
#'

standardize_data <- function(dataType, channel = channel, csv, csvCols, yrRange){

  if(dataType %in% c('NESurveys', 'NEObserver', 'CSV')){

    ##now move onto different processing pipelines
    if(dataType == 'Surveys'){
      print('Standardizing Survey Data...')
      ## pull fisheries-independent data
      data <- survdat::get_survdat_data(channel, getWeightLength = F, getLengths = F, getBio = F, conversion.factor = T)
      surv <- data$survdat
      rm(data) #we only need surv, so clearing out everything else to save room

      #pull species list to help with matching
      spp.qry <-  paste0("select SCINAME, COMNAME, SVSPP
        from svdbs.SVSPECIES_LIST
        order by SVSPP")
      survSPP <- data.table::as.data.table(DBI::dbGetQuery(channel, spp.qry))

      surv <- merge(surv, survSPP)
      surv$MONTH <- lubridate::month(surv$EST_TOWDATE) #create month column since surv data doesn't come with one
      surv$towID <- paste(surv$CRUISE6, surv$STRATUM, surv$TOW, surv$STATION)
      dat <- surv[,c('towID', 'YEAR', "MONTH", "LON", "LAT", "ABUNDANCE", "SCINAME")] #subset to necessary columns
      names(dat) <- c('towID', 'year', 'month', 'lon', 'lat', 'count', 'name') #standardize names
      dat <- dat[dat$year >= yrRange[1] & dat$year <= yrRange[2],]
    } #end if survey

    if(dataType == 'Observer'){
      print('Standardizing Observer Data...')
      ## pull fisheries-dependent data - a bit more intense since there isn't a nice function to do it, and needs to be subset, but follows the same basic steps as survdat
      #observer data
      obs.qry <-  paste0("select YEAR, MONTH, TRIPID, HAULNUM, LONHBEG, LATHBEG, NESPP4, HAILWT
        from obdbs.OBSPP
        where YEAR between ", yrRange[1], " and ", yrRange[2], "
        order by YEAR, MONTH, TRIPID")
      obs <- data.table::as.data.table(DBI::dbGetQuery(channel, obs.qry))

      #at sea monitor data
      asm.qry <- paste0("select YEAR, MONTH, TRIPID, HAULNUM, LONHBEG, LATHBEG, NESPP4, HAILWT
        from obdbs.ASMSPP
        where YEAR between 1993 and 2019
        order by YEAR, MONTH, TRIPID")
      asm <- data.table::as.data.table(DBI::dbGetQuery(channel, asm.qry))

      #combine them
      ob.asm <- rbind(obs, asm)

      #species query
      spp.qry <- paste0("select NESPP4, COMNAME, SCINAME
        from obdbs.OBSPEC")
      obsSPP <- data.table::as.data.table(DBI::dbGetQuery(channel, spp.qry))

      #merge datasets
      obsdat <- merge(ob.asm, obsSPP, by = 'NESPP4')
      #remove NAs
      obsdat <- na.omit(obsdat)

      #fix formatting for YEAR & MONTH (other variables could/should be fixed too but these are the ones we're using at the moment)
      obsdat$YEAR <- as.numeric(obsdat$YEAR)
      obsdat$MONTH <- as.numeric(obsdat$MONTH)

      ###converting position to decimal degrees
      #split out components
      #longitude
      obsdat$LON_DEG <- substr(obsdat$LONHBEG, 1, 2)
      obsdat$LON_MIN <- substr(obsdat$LONHBEG, 3, 4)
      obsdat$LON_SEC <- substr(obsdat$LONHBEG, 5, 6)
      #convert to decimal degrees
      obsdat$LONDD <- as.numeric(measurements::conv_unit(paste(obsdat$LON_DEG, obsdat$LON_MIN, obsdat$LON_SEC, sep = ' '), from = 'deg_min_sec', to = 'dec_deg'))
      obsdat$LONDD <- obsdat$LONDD * -1 #multiply by negative 1 to get W

      #latitude
      obsdat$LAT_DEG <- substr(obsdat$LATHBEG, 1, 2)
      obsdat$LAT_MIN <- substr(obsdat$LATHBEG, 3, 4)
      obsdat$LAT_SEC <- substr(obsdat$LATHBEG, 5, 6)
      #convert to decimal degrees
      obsdat$LATDD <- as.numeric(measurements::conv_unit(paste(obsdat$LAT_DEG, obsdat$LAT_MIN, obsdat$LAT_SEC, sep = ' '), from = 'deg_min_sec', to = 'dec_deg'))

      obsdat$ID <- paste0(obsdat$TRIPID, obsdat$HAULNUM)

      dat <- obsdat[,c('ID', 'YEAR', "MONTH", 'LONDD', 'LATDD', 'HAILWT', 'SCINAME')] #subset to important columns
      names(dat) <- c('towID', 'year', 'month', 'lon', 'lat', 'count', 'name') #standardize names
    } #end if observer

    if(dataType == 'CSV'){
      ### build raster from CSV data from state/other surveys
      print('Standardizing CSV Data...')

      s <- read.csv(csv) #read in each csv
      s$time <- as.POSIXct(s[,csvCols[4]]) #pull time
      s$month <- lubridate::month(s$time) #make month
      s$year <- lubridate::year(s$time) #make year

      dat <- s[,c('time', 'month', 'year', csvCols)] #subset to important columns
      colnames(dat) <- c('time', 'month', 'year', 'towID', 'lon', 'lat', 'date', 'count', 'name')
      dat <- dat[dat$year >= yrRange[1] & dat$year <= yrRange[2],]
    }

    return(dat)

  }  else {
    stop('dataType is not NESurveys, NEObserver, or CSV')
  }

}
