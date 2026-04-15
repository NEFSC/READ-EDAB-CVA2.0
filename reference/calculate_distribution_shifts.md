# Calculate Change in Distribution Metrics

Calculate two metrics to help describe changes in species
distributions - Center of Gravity and Area of Preferred Habitat

## Usage

``` r
calculate_distribution_shifts(abund, area_threshold = 0.75, cell_area = 8 * 8)
```

## Arguments

- abund:

  model predictions, output from `make_predictions` as a rasterStack.
  Layers need names that correspond to the time stamps represented.

- area_threshold:

  a value between 0 and 1, defaults to 0.75. The total area of grid
  cells with the probability of occurrence of a species equal to or
  above this threshold will be calculated.

- cell_area:

  a number representing the area of a single model grid cell on which
  the model is predicted.

## Value

a data.frame with three columns: 1) timestamp, equal to the names of the
layers in `abund`; 2) COG; 3) Area of probabilities greater than the
`area.threshold`
