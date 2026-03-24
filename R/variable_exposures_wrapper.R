#' @title Make Timeseries and Maps of Species-Specific Exposure for each variable
#' @description A wrapper function for \code{make_variable_exposure} that handles object loading, and generating mean SDMs from timeseries, and produces both maps and timeseries
#'
#' @param spp species name. Used to pull correct data and save outputs in species-specific folders.
#' @param ens_name name of ensemble species distribution model file to pull
#' @param sdm_threshold value between 0 and 1. Will remove values lower than this threshold from average ensemble model results to help reduce weird aliasing that can occur in workflow. Defaults to 0.1.
#' @param present_time,future_time character strings indicating the present and future time series to compare. Example: '1993-2019'. Used to pull correct ranked exposure values and save the data properly
#'
#' @return Nothing is returned. The outputs from \code{make_variable_exposure(type = 'map')} and \code{make_variable_exposure(type = 'timeseries')} are saved in the appropriate folders

variable_exposures_wrapper <- function(spp, ens_name, sdm_threshold = 0.1, present_time, future_time){
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
  load(paste0('./', spp, '/output_rasters/', ens_name, '.RData')) #abund

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
  avgHSM[avgHSM < sdm_threshold] <- NA

  ###to add: subset avgHSM here by stock area if needed to calculate exposure across time and space appropriately

  #load in ranked data
  load(paste0('./RawExposure/Data/', present_time, ' vs ', future_time, '_exposure_ranked.RData')) #expRanked

  #map
  mapExp <- make_variable_exposure(type = 'map', ranked_exposure = expRanked, sdm_raster = avgHSM)
  save(mapExp, file = paste0(file.path(getwd(),spp, 'Data'), '/', present_time, ' vs ', future_time, '/variable_exposure_maps.RData'))

  #timeseries
  vecExp <- make_variable_exposure(type = 'timeseries', ranked_exposure = expRanked, sdm_raster = avgHSM)
  save(vecExp, file = paste0(file.path(getwd(),spp, 'Data'), '/', present_time, ' vs ', future_time, '/variable_exposure_timeseries.RData'))
}
