# Calculate Data Quality Scores for each Attribute & Species

Calculates mean data quality score

## Usage

``` r
calculate_data_quality(species_attributes)
```

## Arguments

- species_attributes:

  a data frame containing all of the expert scores for one species. See
  example for how to generate this list from the FCVA output.

## Value

a vector of mean data quality scores for each attribute with the length
of the number of attributes.

## Examples

``` r
if (FALSE) { # \dontrun{
data <-utils::read.csv('expert_scores.csv') #sensitivity data from FSCVA Portal
species.data.list <- split(data, data$Stock.Name)

species.dqs <- lapply(species.data.list, calculate_data_quality)

speciesDQ <- do.call(rbind, species.dqs)
#this combines the list of vectors into a data.frame
} # }
```
