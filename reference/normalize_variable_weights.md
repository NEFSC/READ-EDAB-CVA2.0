# Normalize Dynamic Variable Weights

Normalizes the component model variable importance, as well as the model
weights in the final ensemble, to create a weighted average vector of
variable importance for the ensemble model. Options allow for static
variables, such as time and space, or non-dynamic environmental
variables such as bathymetry, to be excluded from the normalization and
weighted average.

## Usage

``` r
normalize_variable_weights(vars, ens_weights, imp_flist)
```

## Arguments

- vars:

  character vector of variables to generate normalized and averaged
  ensemble importance for

- ens_weights:

  vector of component model weights in the ensemble

- imp_flist:

  character vector of file paths to variable importance for component
  models

## Value

A vector representing the weighted average of normalized variable
importance, representing variable importance in the final ensemble SDM.
