# Merge Fisheries Presence/Absence/Effort Rasters

Combines rasters from `create_rast` across multiple sources

## Usage

``` r
merge_fisheries_rasters(raster_list)
```

## Arguments

- raster_list:

  list of outputs from `create_rast` to combine

## Value

returns a single rasterBrick with the same spatial and temporal
resolution as provided rasters. Maximum value should be 2 to ensure that
presences are included.
