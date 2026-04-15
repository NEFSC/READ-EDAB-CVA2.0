# Match Species and Model Rasters and Create a Data Frame

Match species and environmental rasters to build data frame for modeling

## Usage

``` r
match_pa_model_rasters(
  pa_rasters,
  model_data,
  add_static = TRUE,
  static_variables
)
```

## Arguments

- pa_rasters:

  rasterStack representing fisheries presence, absence, and effort with
  a layer for each timestep. Layers should be named by timestamp, for
  example, month.year

- model_data:

  A list containing rasterStacks of each environmental variable to be
  matched, or the output of `norm_env`. Each member of the list should
  have a name which will correspond to the environmental variable's
  column name in the resulting dataframe. Each rasterStack should have
  the same number of layers as `pa_rasters`, and each layer should have
  a name that can be matched to the names in `pa_rasters`

- add_static:

  TRUE/FALSE to add static variables such as bathymetry. If TRUE,
  staticData needs to be supplied.

- static_variables:

  list of rasters containing static variables such as bathymetry. Names
  of the listed rasters will become the column names.

## Value

`merge_pa_model_rasters` returns a data frame with the following
columns:

- `x, y, month, year` indicating longitude, latitude, month, and year

- `value` - the column indicating presence, absence, and effort from the
  fisheries raster

- A number of columns equal to the length of `model_data` with the same
  names containing the corresponding environmental data at each location

- A number of columns equal to the length of `static_variables`, if
  `add_static = T`, with column names matching the names of the rasters
  in the list
