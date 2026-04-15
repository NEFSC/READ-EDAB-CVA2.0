# Convert Species/Environmental Data Frame to Presence/Absence

This function converts the presence/absence/effort (0-2) data to
presence/absence data (0-1)

## Usage

``` r
set_binary_pa(se, pa_col)
```

## Arguments

- se:

  data.frame of fisheries and environmental data

- pa_col:

  column name for presence/absence column

## Value

a data frame where presence/absence has been set to 0 for absent and 1
for present
