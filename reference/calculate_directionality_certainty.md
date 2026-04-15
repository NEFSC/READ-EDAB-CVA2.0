# Calculate Certainty on Bootstrapped Directionality Scores

Calculates certainty - the percentage of bootstrapped values that match
expert scores - on directionality

## Usage

``` r
calculate_directionality_certainty(bootstrap_scores, directionality_scores)
```

## Arguments

- bootstrap_scores:

  output from `calculate_directionality(bootstrap = TRUE)`

- directionality_scores:

  output from `calculate_directionality(bootstrap = FALSE)` to append
  bootstrap certainty to.

## Value

returns a vector with a length of two: 1) the output from
`calculate_directionality(bootstrap = FALSE)` and 2) the certainty
value, which is the percentage of bootstrapped sensitivities that
matched the weighted average final directionality
