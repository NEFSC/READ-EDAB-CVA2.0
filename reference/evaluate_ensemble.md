# Test the Ensemble Model with Observations

Calculate Area under the Curve (AUC) for the ensemble species
distribution model on an external set of data (not used for building the
model). This function builds a dataset of presence/absence with the new
data across a variety of sources, extracts the predicted values from the
ensemble, and calculates the resulting AUC.

## Usage

``` r
evaluate_ensemble(spp, spp_names, sources, yr_min, yr_max)
```

## Arguments

- spp:

  species name. Used to help pull correct ensemble model.

- spp_names:

  a vector containing all possible names for the target species. Must
  have a length \>= 1

- sources:

  a vector containing all of the source csvs containing the observation
  data to use in the test. These are loaded in directly by `read.csv` so
  they should either be full file paths or in the current working
  directory.

- yr_min:

  start of year range to specify which predictions to use

- yr_max:

  end of year range to specify which predictions to use

## Value

returns the AUC value for the given model on the new data
