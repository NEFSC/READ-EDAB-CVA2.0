# Make Sensitivity Barplot Reports

Generates a set of PDFs, with each page containing stacked barplots
illustrating the distribution of tallies for each sensitivity attribute
for each species scored. The option to generate scorer-specific reports
is available.

## Usage

``` r
make_sensitivity_barplots(
  data,
  species,
  plots_name,
  sensitivity = TRUE,
  plot_data_quality = FALSE,
  plot_legend = TRUE,
  preliminary = TRUE,
  scorer = TRUE,
  scorer_dir
)
```

## Source

This function was originally written as two seperate functions to
generate 1) a PDF report of all scores and 2) individual reports for
each scorer for their assigned species. These were originally written by
Megan Stachura in December 2014 for the original CVA analysis in the
Northeast US. The functions were combined in 2026 by the package author.

## Arguments

- data:

  all data exported from database

- species:

  list of all the species in the order you want them to appear in the
  plots

- plots_name:

  desired name of plot file. file extension not needed.

- sensitivity:

  TRUE indicates sensitivity attribute scores should be plotted. FALSE
  indicates exposure factor scores should be plotted

- plot_data_quality:

  TRUE indicates the average data quality for each attribute should be
  printed on the plot, FALSE indicates this should not be printed on the
  plot

- plot_legend:

  TRUE indicates a legend of the scorers should be printed on the plots,
  FALSE indicates this legend should not be printed

- preliminary:

  TRUE indicates that the "Preliminary Results" label should be printed
  at the top of the plots, FALSE indicates this should not be printed on
  the plots

- scorer:

  TRUE indicates that individual PDF reports for each scorer should also
  be produced.

- scorer_dir:

  the directory to save individual scorer PDF reports to

## Value

No returns. Saves pdf file of the stacked barplots, with one page of
barplots per species scored by the corresponding experts. PDFs will have
naming convention FIRSTNAME_LASTNAME.pdf and be saved in the designated
folder name
