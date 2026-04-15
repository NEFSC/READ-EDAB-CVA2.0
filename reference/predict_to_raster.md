# Convert Predicted Values in Data.Frame to Rasters

Make rasters from predicted values in a data.frame. Helps speed up
sdmTMB predictions.

## Usage

``` r
predict_to_raster(df, static_variables)
```

## Arguments

- df:

  the output from `raster_to_df`, with a column `est` containing the
  predicted values. See example.

- static_variables:

  list of rasters containing the static variables used in model. Should
  be the same object used in `merge_spp_env`. Used to set resolution of
  output rasters

## Value

a rasterStack with the same number of layers as in `rasts` with
probability of occurance predicted for each day with available data.
Values will range from 0 to 1.

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
