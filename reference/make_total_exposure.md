# Make Map or Timeseries of Total Species-Specific Exposure

Uses the logic rule from CVA1.0 to combine variable-specific exposures
to create a total exposure map or timeseries

## Usage

``` r
make_total_exposure(
  type,
  variable_exposure,
  count_all,
  weights,
  weight_threshold
)
```

## Arguments

- type:

  designates desired output, must equal 'map' or 'timeseries'

- variable_exposure:

  If `type == 'map'`, a rasterStack output from
  `make_variable_exposure(type == 'map')`. If `type == 'timeseries'`,
  the matrix output from `make_variable_exposure(type == 'timeseries')`.

- count_all:

  TRUE/FALSE to use `weights` and `wThreshold` to subset variables to
  only important variables

- weights:

  output from `combine_weights` - a vector of variable weights in
  ensemble SDM

- weight_threshold:

  numeric value used to subset weights, variables with weights less to
  or equal to this value will be excluded from total exposure
  calculation

## Value

If `type == 'map'`, the output is a raster representing total exposure
across space. If `type == 'timeseries'`, the output is a vector
representing total exposure across time.
