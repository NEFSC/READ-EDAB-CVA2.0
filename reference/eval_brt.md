# Calculate BRT model evaluation statistics

Calculate evaluation statistics for BRT models. From Camrin Brawn
(WHOI): https://zenodo.org/records/7971532. Included to ensure
functionality.

## Usage

``` r
eval_brt(model, test_data, response, plot = TRUE)
```

## Source

Much of this code is derived from https://github.com/elhazen/PA-paper

## Arguments

- model:

  a BRT model object

- test_data:

  is data used to fit the model. Must contain a column name that matches
  response.

- response:

  is character providing the name of the response variable in test_data

- plot:

  TRUE/FALSE to turn on/off plotting

## Value

model evaluation statistics
