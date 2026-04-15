# Calculate Standard Deviation on MOM6 Data

Calculate monthly standard deviations from raw model data across the
entire provided timeseries. This is built specifically for MOM6 output,
but would work on any rasterStack of gridded data.

## Usage

``` r
sd_model_data(raw_list)
```

## Arguments

- raw_list:

  List of output rasters from `pull_mom6_hindcast` or
  `pull_mom6_forecast`

## Value

a list whose length is equal to the number of variables supplied, where
each item in the list is a rasterStack of data associated with that
variable
