# Calculate Variable Importance

Calculate variable importance for one of the component models

## Usage

``` r
calculate_sdm_variable_importance(
  mod,
  se,
  pa_col,
  xy_col,
  month_col,
  year_col,
  model
)
```

## Arguments

- mod:

  the output from `make_sdm` - only used for ensemble

- se:

  data frame containing species presence/absence data and desired
  environmental covariate data.

- pa_col:

  column name for presence/absence column

- xy_col:

  a vector with a length of 2 indicating the longitude and latitude
  column names

- month_col, year_col:

  column names for month and year columns respectively

- model:

  one of the following indicating the desired model to calculate
  variable importance for: gam, maxent, brt, rf, or sdmtmb

## Value

a vector of the variable importance for the given model. Each model
calculates these differently, so the values should be normalized in
order to compare across models.
