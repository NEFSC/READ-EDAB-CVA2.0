# Calculate Certainty on Bootstrapped Sensitivity Scores

Calculates certainty - the percentage of bootstrapped values that match
expert scores - on sensitivity

## Usage

``` r
calculate_sensitivity_certainty(bootstrap_scores, expert_scores)
```

## Arguments

- bootstrap_scores:

  output from `calculate_sensitivity(bootstrap = TRUE)`

- expert_scores:

  output from `calculate_sensitivity(bootstrap = FALSE)` to append
  bootstrap certainty to.

## Value

returns the `expert_scores` with a new column called 'Certainty', which
includes the percentage of bootstrapped sensitivities that matched the
weighted average final sensitivity.
