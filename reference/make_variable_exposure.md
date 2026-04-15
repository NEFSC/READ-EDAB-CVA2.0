# Make Maps or Timeseries of Variable-Specific Exposure

Combines raw variable exposures and SDM results, using the SDM results
as weights to average the exposure across time to produce maps of
species-specific exposure for each variable or across space to produce
timeseries of species-specific exposure for each variable

## Usage

``` r
make_variable_exposure(type, ranked_exposure, sdm_raster)
```

## Arguments

- type:

  designates desired output, must equal 'map' or 'timeseries'

- ranked_exposure:

  output from `rank_exposure`

- sdm_raster:

  rasterStack with monthly averages of SDM results to use as weights for
  weighted average

## Value

If `type == 'map'`, a rasterStack with the number of layers equal to the
number of variables supplied, containing the weighted average exposure,
weighted by SDM results, for each variable. If `type == 'timeseries'`, a
matrix with the number of rows equal to the number of variables
supplied, and 12 columns (1 for each month) containing the weighted
average exposure, weighted by SDM results, for each variable.
