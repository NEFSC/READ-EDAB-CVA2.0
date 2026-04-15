# Component Model Cross-Validation

Perform cross-validation on one of the five component model types. k is
5 for all models.

## Usage

``` r
cross_validate_sdm(mod, se, pa_col, xy_col, month_col, year_col, model)
```

## Arguments

- mod:

  the output from `build_sdm`

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

  one of the following indicating the desired model to perform cross
  validation on: gam, maxent, brt, rf, or sdmtmb

## Value

the cross-validation results from the desired model. Objects will differ
slightly depending on model type.
