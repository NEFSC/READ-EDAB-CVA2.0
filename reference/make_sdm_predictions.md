# Predict Component or Ensemble SDM

Make model predictions for either a component model
(GAM/MAXENT/RF/BRT/SDMTMB) or the Ensemble model

## Usage

``` r
make_sdm_predictions(
  mod,
  model,
  rasts,
  static_variables,
  bathy_raster,
  mask = T,
  bathy_max,
  se = NULL,
  month_col,
  year_col,
  xy_col,
  weights = NULL
)
```

## Arguments

- mod:

  the output from `make_sdm` - only used for ensemble

- model:

  one of the following indicating the desired model to calculate
  variable importance for: gam, maxent, brt, rf, sdmtmb, or ens

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

- se:

  data frame containing species presence/absence data and desired
  environmental covariate data.

- month_col, year_col:

  column names for month and year columns respectively

- xy_col:

  a vector with a length of 2 indicating the longitude and latitude
  column names

- weights:

  a vector of model weights - used for building the ensemble model

## Value

returns a rasterStack of predicted habitat suitability. The number of
layers will be equal to the number of layers in `rasts`
