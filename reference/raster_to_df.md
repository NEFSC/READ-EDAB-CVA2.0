# Convert Environmental Rasters to a Data.Frame

Make a data frame out of environmental data to be used for predictions.
Helps speed up sdmTMB predictions.

## Usage

``` r
raster_to_df(rasts, static_variables, bathy_raster, bathy_max, mask)
```

## Arguments

- rasts:

  list of rasterStacks corresponding to the environmental covariates
  used to build the models. The number of layers in each rasterStack
  should be the same and correspond to the length of the timeseries for
  the models to be predicted on

- static_variables:

  list of rasters containing the static variables used in model. Should
  be the same object used in `merge_spp_env`.

- bathy_raster, bathy_max:

  a raster of bathymetry with the same extent and resolution as `rasts`
  and the maximum depth you want included. For example, if you want to
  mask off waters deeper than 1000 m, `bathy_max` would be set to 1000.
  The value should be positive regardless of the sign of your bathymetry
  data.

- mask:

  TRUE/FALSE indicating whether to mask off certain depths (i.e. waters
  deeper than 1000 m)

## Value

a data.frame containing all of the environmental data in `rasts` for all
timesteps. It will be big depending on the length of your time series.
This is meant to be used to help predict sdmTMB models.

## Examples

``` r
if (FALSE) { # \dontrun{
#create data frame of environmental data
allDF <- raster_to_df(rasts = rasts, static_variables = staticVars,
bathy_raster = bathyR, bathy_max = 1000, mask = T)

#predict sdmTMB model for all data; not for each timestep as in \code{make_sdm_predictions}
#this generates the \code{est} column
pred <- stats::predict(mod, newdata = allDF, type = 'response')

#add appropriate month.year (my) column to create rasters
pred$my <- paste(pred$month, pred$year, sep = '.')
abund <- predict_to_raster(df = pred, staticData = staticVars) #make into rasters
} # }
```
