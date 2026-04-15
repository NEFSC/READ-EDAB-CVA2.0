# Bhattacharyya's coefficient test

test Bhattacharyya's coefficient ranges from 0 to 1, with 0 being no
overlap in distributions and 1 being perfect overlap. Used within
`eval_brt`. Included to ensure functionality.

## Usage

``` r
bhattacharyya_stat_brt(data, response, vars)
```

## Source

From Camrin Brawn (WHOI): https://zenodo.org/records/7971532 and
http://tguillerme.github.io/R/bhatt.coef.R

## Arguments

- data:

  data used to test the model. Must have a presence absence column
  called "presabs" and column names that match the names in vars

- response:

  name of the response variable

- vars:

  character vector of the desired variables

## Value

vector of Bhattacharyya's coefficients for each variable
