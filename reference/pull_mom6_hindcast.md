# Pull MOM6 Hindcast Data

Pull hindcast data from the MOM6 model output based on provided URL from
the CEFI portal

## Usage

``` r
pull_mom6_hindcast(
  var_url,
  req_vars,
  short_names,
  gt = "regrid",
  of = "monthly",
  bounds = c(-78, -65, 35, 45),
  static,
  release
)
```

## Arguments

- var_url:

  URL pointing to JSON table variable lists for desired MOM6 hindcast
  and domain

- req_vars:

  vector of variable names to pull. Must match names in the
  'cefi_long_name' column provided JSON table

- short_names:

  vector of simplified variable names to help name resulting raster
  files. Must be the same length as req_vars.

- gt:

  desired grid type. Must match one of the options in the
  'cefi_grid_type' column in provided JSON table

- of:

  desired output frequency. Must match one of the options in the
  'cefi_output_frequency' column in provided JSON table

- bounds:

  xmin, xmax, ymin, ymax of desired output raster

- static:

  URL to static grid for MOM6

- release:

  release code. Must match one of the options in the 'cefi_release'
  column in provided JSON table

## Value

a list whose length is equal to the number of variables supplied, where
each item in the list is a rasterStack of data associated with that
variable
