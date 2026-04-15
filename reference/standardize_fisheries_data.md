# Standardize Fisheries Dependent and Independent Datasets

Pulls and standardizes fisheries independent and dependent datasets from
the NEFSC survey and observer programs using ROracle. The function also
can open a local CSV file and standardize the output.

## Usage

``` r
standardize_fisheries_data(
  data_type,
  channel = channel,
  csv,
  csv_columns,
  yr_range
)
```

## Arguments

- data_type:

  type of data to be standardize. Must be one of the following:
  'NESurveys', 'NEObserver', or 'CSV'

- channel:

  connection to remote databases. Only required for 'NESurveys' or
  'NEObserver'

- csv:

  path to local CSV file

- csv_columns:

  Column names in csv file in the following order: 'towid', 'longitude',
  'latitude', 'date', 'count (can be count/abundance/density, etc)',
  'name'. Date must be in a format that can be converted to POSIX with
  as.POSIXct

- yr_range:

  A vector with length of 2 indicating the start and end, inclusive,
  year of the desired time series

## Value

a data frame. It is recommended to save this as a csv file as pulling
survey and observer datasets does take time
