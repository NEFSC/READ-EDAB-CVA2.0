# K-fold BRT fit

k-fold cross validation code from Camrin Brawn (WHOI):
https://zenodo.org/records/7971532. Function originally by B. Abrahms
from
https://github.com/briana-abrahms/DynamicEnsembleSDM/blob/master/model_evaluation.R.
Modified to return information necessary to calculate evaluation
metrics. Follows many of the same inputs as `gbm.step`

## Usage

``` r
eval_kfold_brt(
  data_input,
  gbm_x,
  gbm_y,
  learning_rate = 0.05,
  k_folds = 5,
  tree_complexity = 3,
  bag_fraction = 0.6,
  is_fixed = TRUE,
  max_trees = 2000
)
```

## Arguments

- data_input:

  input data frame

- gbm_x:

  names of predictor variables in `data_input`

- gbm_y:

  name of response variable in `data_input`

- learning_rate:

  sets the weight applied to individual trees input to dismo::gbm.step

- k_folds:

  number of folds to perform

- tree_complexity:

  is tree complexity input to dismo::gbm.step - sets complexity of
  individual trees

- bag_fraction:

  is bag fraction input to dismo::gbm.step - sets the proportion of
  observations used in selecting variables

- is_fixed:

  TRUE/FALSE to determine if gbm.step or gbm.fixed should be used.

- max_trees:

  maximum number of trees to fit before stopping

## Value

a list containing 1) the output from `eval_brt` and 2) a data.frame of
the observed and predicted values used to calculate RMSE/AUC by `sdm_cv`

## Examples

``` r
if (FALSE) { # \dontrun{
eval_kfold_brt(data_input_Fit, gbm_x=c("curl","ild", "ssh", "sst","sst_sd"), "presabs")
} # }
```
