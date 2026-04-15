# Rank Raw Environmental Exposure

Converts raw exposure from `calculate_raw_exposure` to a rank between
1-4 based on thresholds used in CVA1.0. This function also flips
variables as needed such that positive exposure is always negative (i.e.
both increasing temperature and decreasing oxygen concentrations result
in positive exposure).

## Usage

``` r
rank_exposure(
  exposure,
  flip = T,
  no_flip_vars = c("bottomT", "surfaceT", "bottomArg", "MLD")
)
```

## Arguments

- exposure:

  output from `calculate_raw_exposure`

- flip:

  TRUE/FALSE option to multiply exposure by -1 to ensure that positive
  exposure values represent negative variable change (increasing
  temperatures, decreasing oxygen concentrations, etc).

- no_flip_vars:

  character vector naming which vectors should NOT be flipped -
  counterintuitive, yes, but this list was actually shorter than the
  variables that needed to be flipped for NECVA2.0. Names should match
  the names in `exposure`

## Value

A list of rasterStacks with each layer containing ranked values between
1 - 4. The length of the list is equal to the length of the lists
supplied as `exposure`
