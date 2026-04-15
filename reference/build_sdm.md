# Build Component or Ensemble Species Distribution Models

Build one of the five component models or the final ensemble model

## Usage

``` r
build_sdm(
  se,
  pa_col,
  xy_col,
  month_col,
  year_col,
  model,
  ensemble_weights = NULL,
  ensemble_preds = NULL
)
```

## Arguments

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

  one of the following indicating the desired model to build: gam,
  maxent, brt, rf, sdmtmb, or ens

- ensemble_weights:

  vector of model weights used to build ensemble model. The vector must
  have the same length and be in the same order as the corresponding
  predictions list.

- ensemble_preds:

  a list of prediction values from component models used to build
  ensemble model. The list of predictions must have the same length and
  be in the same order as the corresponding weight vector.

## Value

the model object from the desired model. Model objects will differ
depending on the model type.
