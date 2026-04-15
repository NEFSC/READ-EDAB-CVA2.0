# Convert Standardized Fisheries Data into Presence/Absence/Effort Raster

Takes output from `standardize_data` and converts it into a raster with
the spatial and temporal resolution of the provided grid ranging from
0-2 where:

- 0 - No fishing effort in the cell

- 1 - Fishing effort present for the target species in the cell

- 2 - Target species caught in the cell

## Usage

``` r
build_fisheries_raster(
  data,
  is_obs = FALSE,
  grid,
  tm_multiplier = 24 * 60 * 60,
  origin = "1993-01-01",
  all_names
)
```

## Arguments

- data:

  data frame to convert to presence/absence raster

- is_obs:

  TRUE/FALSE indicating whether or not the data is observer or similar
  fisheries-dependent data. Will force function to check if species is
  observed at least 30 times throughout timeseries before creating
  raster

- grid:

  static link to a ncdcf object with the variables lon, lat, time - can
  be link to remote data - must be able to be read with nc_open

- tm_multiplier:

  multiplier to help convert timestep to POSIX (seconds since origin),
  defaults to 86400 (number of seconds in a day)

- origin:

  Origin of time series to be used by POSIXct

- all_names:

  a vector containing all possible names for the target species. Must
  have a length \>= 1

## Value

a rasterBrick with the same extent as the provided grid, and a number of
layers equal to the timeseries associated with the provided model data
