# Calculate Area under the Curve (AUC) for BRT models

Calculate AUC for BRT k-folds. Used within `eval_brt`. Included with
CVA2.0 package to ensure functionality.

## Usage

``` r
save_auc_brt(truth, predicted, plot_roc = FALSE, ...)
```

## Source

From Camrin Brawn (WHOI): https://zenodo.org/records/7971532.

## Arguments

- truth:

  observations

- predicted:

  predicted values

- plot_roc:

  TRUE/FALSE to turn on/off plotting

- ...:

  arguments passed to `ROCR::plot`

## Value

an objected of class `performance`. Same as
[`ROCR::performance`](https://ipa-tys.github.io/ROCR/reference/performance.html).
