#' @title Make Timeseries and Maps of Species-Specific Exposure for each variable
#' @description A wrapper function for \code{makeExposureMaps} and \code{MakeExposureTimeseries} that handles object loading, and generating mean SDMs from timeseries
#'
#' @param spp species name. Used to pull correct data and save outputs in species-specific folders.
#' @param ensName name of ensemble species distribution model file to pull
#' @param sdmThreshold value between 0 and 1. Will remove values lower than this threshold from average ensemble model results to help reduce weird aliasing that can occur in workflow. Defaults to 0.1.
#' @param presentTime,futureTime character strings indicating the present and future time series to compare. Example: '1993-2019'. Used to pull correct ranked exposure values and save the data properly
#'
#' @return Nothing is returned. The outputs from \code{makeExposureMaps} and \code{MakeExposureTimeseries} are saved in the appropriate folders

makeVariableAverages <- function(spp, ensName, sdmThreshold = 0.1, presentTime, futureTime){
  #open log file
  sink(file = file.path(getwd(), 'logs', 'variable_averages.log'), append = T)

  # Ensure the sinks are closed when the function exits, regardless of how it exits.
  on.exit({
    #sink(type = "message")
    sink()
  })

  print(Sys.time())
  print(spp)

  #load SDM results and average monthly
  load(paste0('/home/kgallagher/ClimateVulnerabilityAssessment2.0/SDMs/', spp, '/output_rasters/', ensName, '.RData')) #abund

  #avg ensemble HSM
  avgHSM <- vector(mode = 'list', length = 12)
  for(y in 1:12){
    mn <- seq(y, length(abund), by = 12)
    MNS <- raster::stack(abund[mn])
    avgHSM[[y]] <- raster::calc(MNS, fun = mean, na.rm = T)
  } #end for
  avgHSM <- raster::stack(avgHSM)
  names(avgHSM) <- month.abb

  ###remove hsm with less than 0.1 to avoid weird aliasing
  avgHSM[avgHSM < sdmThreshold] <- NA

  ###to add: subset avgHSM here by stock area if needed to calculate exposure across time and space appropriately

  #load in ranked data
  load(paste0('./RawExposure/Data/', presentTime, ' vs ', futureTime, '_exposure_ranked.RData')) #expRanked

  #map
  mapExp <- makeExposureMaps(rankExp = expRanked, sdmRast = avgHSM)
  save(mapExp, file = paste0(file.path(getwd(),spp, 'Data'), '/', presentTime, ' vs ', futureTime, '/variable_exposure_maps.RData'))

  #timeseries
  vecExp <- makeExposureTimeseries(rankExp = expRanked, sdmRast = avgHSM)
  save(vecExp, file = paste0(file.path(getwd(),spp, 'Data'), '/', presentTime, ' vs ', futureTime, '/variable_exposure_timeseries.RData'))
}
