# Calculate Model Confidence Score for each Species

Calculates mean and standard deviation of the model confidence score

## Usage

``` r
calculate_model_confidence(species)
```

## Arguments

- species:

  a data frame containing all of the expert scores for one species in a
  column called `Score`. See example for how to generate this list from
  the combined model confidence score spreadsheets.

## Value

a data frame including mean and standard deviation of data quality
scores for each attribute with the length of the number of attributes.

## Examples

``` r
if (FALSE) { # \dontrun{
species.data.list <- split(data, data$Species)
#data is the combined data frame of
#all model confidence spreadsheets from scorers

species.conf <- lapply(species.data.list, model.confidence)

speciesMC <- do.call(rbind, species.conf)
#this combines the list of vectors into a data.frame
} # }
```
