# Remove Correlated Environmental Covariates from Species/Environmental Data Frame

This function removed correlated covariates to prepare data for modeling

## Usage

``` r
remove_corr(se, pa_col, xy_col, month_col, year_col)
```

## Arguments

- se:

  data.frame of fisheries and environmental data

- pa_col:

  column name for presence/absence column

- xy_col:

  a vector with a length of 2 indicating the longitude and latitude
  column names

- month_col, year_col:

  column names for month and year columns respectively

## Value

a data frame with correlated covariates removed
