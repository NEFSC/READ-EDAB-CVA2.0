# Normalize Model Data

Normalize model data across the entire provided timeseries using a
z-score

## Usage

``` r
normalize_model_data(raw_list, avg_list, sd_list, short_names)
```

## Arguments

- raw_list:

  List of output rasters from `pull_hind` or `pull_forecast`.
  Alternatively, a list of rasters containing gridded data to be
  normalized, with each entry in the list containing a rasterStack of a
  timeseries of environmental data to be normalized.

- avg_list, sd_list:

  List of output rasters from `avg_model_data` and `sd_env`,
  respectively. Alternatively, averages and standard deviations of
  desired gridded datasets. Should be a list of rasterStacks, with each
  entry in the list containing the rasterStack for a different variable.

- short_names:

  vector of simplified variable names to help name resulting raster
  files.

## Value

a list whose length is equal to the length of `raw_list`, where each
item in the list is a rasterStack of normalized data associated with
that variable
