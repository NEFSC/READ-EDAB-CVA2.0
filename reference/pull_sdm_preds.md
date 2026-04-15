# Pull Predicted Values from SDMs CV

Extract the predictions from the cross-validation of one of the five
component models

## Usage

``` r
pull_sdm_preds(cv, model)
```

## Arguments

- cv:

  output from `cross_validate_sdm`

- model:

  one of the following indicating the desired model to extract predicted
  values from: gam, maxent, brt, rf, or sdmtmb

## Value

a data.frame containing the prediction outputs from the cross-validation
necessary to calculate evaluation metric
