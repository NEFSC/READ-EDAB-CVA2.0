# Make Plots for Species Distribution Models

Produces plots necessary for SDM reports. Requires directory to be set
up per directions in the package documentation/manual.

## Usage

``` r
make_sdm_plots(
  species,
  type = c("ensemble", "weights", "importance", "residuals"),
  yr_min,
  yr_max,
  coastline,
  bathy,
  model_metrics
)
```

## Arguments

- species:

  name of the species to plot. Must match folder name to pull correct
  data and save figures correctly.

- type:

  at least one of the following: 'ensemble', 'weights', 'importance',
  'residuals'. Used to determine what to plot.

- yr_min, yr_max:

  start and end of timeseries to plot, used to pick correct model
  predictions

- coastline:

  shapefile used to plot land in model prediction plots

- bathy:

  raster file of bathymetry data, used to help mask off deeper
  predictions

- model_metrics:

  a data.table containing the model metrics for each species. Generated
  by `make_evaluation_csv`.

## Value

Function does not return anything. Figures are saved to species-specific
`figures` folder.
