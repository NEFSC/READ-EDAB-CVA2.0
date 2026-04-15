# Calculate Performance Metric

Calculate RMSE or AUC for one of the component models or the final
ensemble

## Usage

``` r
evaluate_sdm(preds, model, metric, mod = NULL)
```

## Arguments

- preds:

  output from `pull_sdm_preds`

- model:

  one of the following indicating the desired model to calculate the
  metric for from: gam, maxent, brt, rf, sdmtmb, or ens

- metric:

  either rmse or auc - determines which metric is calculated

- mod:

  the output from `build_sdm` - only used for ensemble

## Value

the desired evaluation metric for the given model
