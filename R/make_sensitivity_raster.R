#' @title Make Sensitivity Rasters
#' @description
#' Map sensitivity scores to a raster
#'
#' @param species name of the species to plot. Must match folder name to pull correct data.
#' @param sensitivity data.frame with attribute scores and final sensitivity and certainty. Generated from \code{calculate.sensitivity} and \code{sensitivity.bootstrap}. Must include the column \code{Stock.Name} to subset data to correct species
#' @param stocks a shapefile containing the stocks to use; if NULL, the same value is applied to the entire raster
#' @param species_col,stock_col name of column with species and stock names to help subset data in \code{sensitivity}
#' @param stock_key a named vector containing the abbreviations and long names of stocks to use as labels
#' @param stock_order a vector containing the long stock names in the desired order
#' @param total_sens_col name of column containing total sensitivity scores in \code{sensitivity}
#' @param r_template a spatRaster of species distribution. Used as a template for the sensitivity raster and as a mask to set unoccupied spaces to NA
#' @param bathymetry a spatRaster of bathymetry data. Must be same extent as r_template
#' @param bathymetry_range a vector containing the minimum and maximum desired bathymetry values
#'
#' @return returns a spatraster of sensitivity scores mapped to species distribution
#'
#'@export
#'
make_sensitivity_raster <- function(
  species,
  sensitivity,
  stocks,
  species_col,
  stock_col,
  stock_key,
  stock_order,
  total_sens_col,
  r_template,
  bathymetry,
  bathymetry_range
) {
  #build sensitivity raster
  #isolate sensitivity data
  spSens <- sensitivity[sensitivity[, species_col] == species, ] #this will be a vector if there are no stocks, but a data.frame if there are

  #failsafe in case the shape files don't exist but spSens has multiple rows
  if (nrow(spSens) > 1 & is.null(stocks)) {
    spSens <- spSens[which(spSens[, stock_col] == 'global'), ] #subset to just global
  }

  #failsafe for if stocks exist but spSens only has one row
  if (nrow(spSens) == 1 & !is.null(stocks)) {
    stocks <- NULL #set stocks to NULL
  }

  #create raster with the same extent as the template raster (should be an exposure raster)
  sensRast <- r_template
  #fill in values accordingly
  if (!is.null(stocks)) {
    #if stocks is not a NULL object, then spSens is a data.frame and we need to find the 'global' row
    #start by making all of the raster the global value
    sensRast[] <- spSens[which(spSens[, stock_col] == 'global'), total_sens_col]
    #match stocks
    stocks_proj <- terra::project(stocks, sensRast)
    for (x in 1:length(stocks_proj)) {
      stock_full_name <- stock_key[
        names(stock_key) %in% stocks_proj$stock_area[x]
      ]
      sensRast <- terra::mask(
        sensRast,
        stocks_proj[x],
        inverse = TRUE,
        updatevalue = spSens[
          which(spSens[, stock_col] == stock_full_name),
          total_sens_col
        ]
      ) #replace value within stock polygon with corresponding sensitivity value
    }
  } else {
    sensRast[] <- spSens[1, total_sens_col] #make all values equal to the global value which will be the only value in the total sensitivity column if stocks is null
  }

  #mask with bathy object
  if (inherits(bathymetry, "PackedSpatRaster")) {
    bathymetry <- terra::unwrap(bathymetry)
  }
  #reproject because there are slight differences in resolution/extent for some reason, especially with the forecasts
  bathy_aligned <- terra::resample(bathymetry, sensRast, method = "bilinear")

  sensRast <- terra::ifel(
    bathy_aligned <= bathymetry_range[1] | bathy_aligned > bathymetry_range[2],
    NA,
    sensRast
  )

  #mask  with template
  sensRast <- terra::mask(sensRast, r_template)

  return(sensRast)
}
