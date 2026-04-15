# Subset Species/Environmental Dataset to Ecologically-Relevant Dynamic Variables

This function reduces the number of environmental covariates in the
species/environmental data frame based on defined feeding and habitat
guilds

## Usage

``` r
match_guilds(
  spp_env,
  spp,
  spp_col = "name",
  spp_guild,
  feeding_key,
  feeding_col = "Feeding.Guild",
  habitat_key,
  habitat_col = "Habitat.Guild",
  static_vars = c("x", "y", "month", "year", "bathy", "rugosity", "dist2coast"),
  pa_col = "value"
)
```

## Arguments

- spp_env:

  data.frame of fisheries and environmental data

- spp:

  Species name to add to log files and save data to correct directory
  (see vignette for recommended directory set up)

- spp_col:

  name of column with species name in spp_guild key

- spp_guild:

  file path to csv file containing a list of species and corresponding
  feeding and habitat guilds. See vignette for recommended set up of
  this file

- feeding_key, habitat_key:

  file path to csv for feeding and habitat guild keys respectively. See
  vignette for recommended set up of these files.

- feeding_col, habitat_col:

  columns names indicating the guild column in the feeding and habitat
  guild csvs

- static_vars:

  a vector with the column names to be retained

- pa_col:

  column name for presence/absence column

## Value

a data frame that contains the static variables listed and only
environmental covariates associated with the species' feeding and
habitat guilds
