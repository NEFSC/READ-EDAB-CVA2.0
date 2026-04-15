# Make Plots of Exposure Results

Produces exposure plots. Requires directory to be set up per directions
in the package documentation/manual.

## Usage

``` r
make_exposure_plots(
  species,
  type = c("variable", "total", "important", "radar"),
  present_time,
  future_time,
  variable_df,
  coastline
)
```

## Arguments

- species:

  names of the species to plot. Must match folder name to pull correct
  data and save figures correctly.

- type:

  at least one of the following: 'variable', 'total', 'important',
  'radar. Used to determine what to plot. Defaults to all.

- present_time, future_time:

  character strings indicating the present and future time series to
  compare. Example: '1993-2019'. Used to pull correct calculations and
  save the data properly

- variable_df:

  a data.frame containing all possible environmental variables, such as
  from the MOM6 model. Must contain columns `Long.Name` and
  `Short.Name`, containing the full names and abbreviated names of the
  variables. Abbreviated names should correspond to those in the weights
  vector produced by `combineWeights`

- coastline:

  shapefile used to plot land in model prediction plots

## Value

Function does not return anything. Figures are saved to species-specific
`figures` folder.
