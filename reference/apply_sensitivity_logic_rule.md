# Calculate Total Sensitivity Score with the Logic Rule

Calculates total sensitivity from expert counts or bootstrapped values.
Based on code from the Highly Migratory Species (HMS) CVA team.

## Usage

``` r
apply_sensitivity_logic_rule(attribute_vec, bootstrap)
```

## Arguments

- attribute_vec:

  output from `attribute_score`

- bootstrap:

  binary TRUE/FALSE to turn on/off bootstrapping of sensitivity scores.
  Will impact outputs.

## Value

returns the final sensitivity scores based on the logic rule. The output
is similar to `attribute_score`, where if `bootstrap = TRUE`, the
function returns a vector containing the final sensitivity score for
each sample, and if `bootstrap = FALSE`, it returns a single final
value.
