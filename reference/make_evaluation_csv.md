# Create a spreadsheet containing all SDM performance metrics

Loads existing performance metrics for component models, and calculates
Area under the Curve (AUC) for the ensemble model 2 ways. This function
also has the option to test the model on an external dataset using
`test_ens`. The final product is a saved csv file containing all of the
performance metrics.

## Usage

``` r
make_evaluation_csv(spp_list, test_ens = T, yr_min, yr_max)
```

## Arguments

- spp_list:

  the data.frame containing species names and alternative names for
  `test_ens`. Must contain the column `Name`, which matches the names of
  the species folders to help pull correct data.

- test_ens:

  a TRUE/FALSE indicating whether or not to test the ensemble model on
  an external dataset

- yr_min:

  start of year range to specify which predictions to use. Only used if
  `test_ens == TRUE`

- yr_max:

  end of year range to specify which predictions to use. Only used if
  `test_ens == TRUE`

## Value

nothing is returned. The resulting CSV file called
'species_evaluation_metrics.csv' is saved in the working directory.

## Details

The resulting CSV file contains the following columns:

- Common.Name:

  name of species

- Managing.Body:

  name of group responsible for species management, if exists

- Feeding.Guild:

  assigned feeding guild

- Habitat.Guild:

  assigned habitat guild

- N.PRESENCE, N.ABSENCE:

  number of presences/absences in the data set used to build the model

- BRT, GAM, MAXENT, RF, SDMTMB:

  AUC for each component model

- BRT.WT, GAM.WT, MAXENT.WT, RF.WT, SDMTMB.WT:

  weights for each component model in the final ensemble

- ENS.AUC:

  Ensemble model AUC based on weighted sum of CV predictions from each
  component model

- WAVG.ENS.AUC:

  Ensemble model AUC calculated as a weighted average of AUCs, using the
  corresponding weights

- AUC.yr_min.yr_max:

  Ensemble model AUC calculated on the external data, if
  `test_ens == TRUE`. Column name will reflect the timeseries used
