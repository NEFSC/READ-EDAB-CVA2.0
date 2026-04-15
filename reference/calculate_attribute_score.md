# Calculate Attribute Scores

Calculates weighted average attribute scores for individual sensitivity
attributes, with or without bootstrapping. Based on code from the Highly
Migratory Species (HMS) CVA team.

## Usage

``` r
calculate_attribute_score(attribute, bootstrap = TRUE, samples = 10000)
```

## Arguments

- attribute:

  a data frame consisting of all expert scores for a single sensitivity
  attribute and species.

- bootstrap:

  binary TRUE/FALSE to turn on/off bootstrapping of sensitivity scores.
  Will impact outputs.

- samples:

  number of samples to run for bootstrapping. Default is 10,000.

## Value

weighted averaged scores for each attribute. If `bootstrap = TRUE`, this
is a data frame with the number of rows equal to the number of samples,
and the number of columns equal to the number of attributes. If
`bootstrap = FALSE`, this is a vector equal to the length of the number
of attributes.
