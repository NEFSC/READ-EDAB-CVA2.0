# Render Species Distribution Model (SDM) reports

Builds SDM report using a provided template. Requires directory to be
set up per directions in the package documentation/manual.

## Usage

``` r
make_sdm_reports(
  species_list,
  yr_min,
  yr_max,
  model_metrics,
  feeding_key,
  habitat_key,
  variable_key,
  template,
  report_path
)
```

## Arguments

- species_list:

  a data.table containing the species names. Must contain the column
  `Common.name`

- yr_min, yr_max:

  start and end of timeseries to plot, used to pick correct model
  predictions to include

- model_metrics:

  a data.table containing the model metrics for each species. Generated
  by `model_metrics`.

- feeding_key, habitat_key:

  file path to csv for feeding and habitat guild keys respectively. Used
  for matching variable names. See vignette for recommended set up of
  these files.

- variable_key:

  a data.frame containing all possible modeled variables, such as from
  the MOM6 model. Must contain columns `Long.Name` and `Short.Name`,
  containing the full names and abbreviated names of the variables.
  Abbreviated names should correspond to those in the variable
  importance plots

- template:

  a character string designating the .qmd file to use as a template

- report_path:

  character string designating where the reports should be saved

## Value

Function does not return anything. Reports are saved to `report_path`
