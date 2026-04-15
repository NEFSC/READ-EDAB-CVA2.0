# Calculate Directionality Score

Calculates weighted average directionality scores for individual
sensitivity attributes, with or without bootstrapping. Based on code
from the Highly Migratory Species (HMS) CVA team.

## Usage

``` r
calculate_directionality(species, bootstrap = TRUE, samples = 10000)
```

## Arguments

- species:

  a data frame consisting of all expert directionality scores for a
  single species.

- bootstrap:

  binary TRUE/FALSE to turn on/off bootstrapping of sensitivity scores.
  Will impact outputs.

- samples:

  number of samples to run for bootstrapping. Default is 10,000.

## Value

weighted averaged directionality score. If `bootstrap = TRUE`, this is a
vector containing bootstrapped scores, with a length equal to the number
of samples. If `bootstrap = FALSE`, this is a single score.
