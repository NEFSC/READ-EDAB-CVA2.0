#' @title Make Sensitivity Barplot Reports
#' @description
#' Generates a set of PDFs, with each page containing stacked barplots illustrating the distribution of tallies for each sensitivity attribute for each species scored. The option to generate scorer-specific reports is available.
#'
#' @source This function was originally written as two seperate functions to generate 1) a PDF report of all scores and 2) individual reports for each scorer for their assigned species. These were originally written by Megan Stachura in December 2014 for the original CVA analysis in the Northeast US. The functions were combined in 2026 by the package author.
#'
#' @param data all data exported from database
#' @param species list of all the species in the order you want them to appear in the plots
#' @param plots_name desired name of plot file. file extension not needed.
#' @param sensitivity TRUE indicates sensitivity attribute scores should be plotted. FALSE indicates exposure factor scores should be plotted
#' @param plot_data_quality TRUE indicates the average data quality for each attribute should be printed on the plot, FALSE indicates this should not be printed on the plot
#' @param plot_legend TRUE indicates a legend of the scorers should be printed on the plots, FALSE indicates this legend should not be printed
#' @param preliminary TRUE indicates that the "Preliminary Results" label should be printed at the top of the plots, FALSE indicates this should not be printed on the plots
#' @param scorer TRUE indicates that individual PDF reports for each scorer should also be produced.
#' @param scorer_dir the directory to save individual scorer PDF reports to
#'
#' @return No returns. Saves pdf file of the stacked barplots, with one page of barplots per species scored by the corresponding experts. PDFs will have naming convention FIRSTNAME_LASTNAME.pdf and be saved in the designated folder name

make_sensitivity_barplots <- function(
  data,
  species,
  plots_name,
  sensitivity = TRUE,
  plot_data_quality = FALSE,
  plot_legend = TRUE,
  preliminary = TRUE,
  scorer = TRUE,
  scorer_dir
) {
  #########################################################
  ### Get information we need from the data file

  # get a list of the attributes, either sensitivity or exposure
  if (sensitivity == TRUE) {
    attributes <- levels(as.factor(data$Attribute.Name[which(
      data$Attributes == 'Sensitivity Attribute'
    )]))
  } else {
    attributes <- levels(as.factor(data$Attribute.Name[which(
      data$Attributes == 'Exposure Factor'
    )]))
  }

  # Get a list of the sensitivity attribute/exposure factor scorers
  if (sensitivity == TRUE) {
    scorers <- levels(as.factor(data$Scorer[which(
      data$Attributes == 'Sensitivity Attribute' & !is.na(data$Scoring.Rank1)
    )]))
  } else {
    scorers <- levels(as.factor(data$Scorer[which(
      data$Attributes == 'Exposure Factor' & !is.na(data$Scoring.Rank1)
    )]))
  }

  # Get the total number of scorers
  num.scorers <- length(scorers)

  #########################################################
  ### Set the colors to choose from for each scorer
  all.colors <- c(
    'red',
    'blue',
    'green3',
    'orange',
    'yellow',
    'gray',
    'cyan',
    'purple',
    'brown',
    'black'
  )

  #########################################################
  ### Set up plotting window

  # open window of specified size
  # windows(w=8.5, h=11)

  # save plots as pdf with specified name
  pdf(
    paste(plots_name, '.pdf', sep = ''),
    pointsize = 12,
    width = 8.5,
    height = 11
  )

  #########################################################
  ### Plot the barcharts for each species, with each page containing barcharts of all the
  ### sensitivity or exposure attributes for a species

  for (i in 1:length(species)) {
    # for each species

    # Set a seed value so that we always get the same result from random sampling for this species
    set.seed(i)

    # start new plotting page with specified margins and plotting parameters
    par(mfrow = c(4, 3), oma = c(3, 3, 6, 1), mar = c(2, 3, 2, 1), las = 1)

    # get all the rows with scores for that species
    if (sensitivity == TRUE) {
      rows <- which(
        data$Stock.Name == species[i] &
          !is.na(data$Scoring.Rank1) &
          data$Attributes == 'Sensitivity Attribute'
      )
    } else {
      rows <- which(
        data$Stock.Name == species[i] &
          !is.na(data$Scoring.Rank1) &
          data$Attributes == 'Exposure Factor'
      )
    }
    species.data <- data[rows, ]

    # get all the scorers who scored for this species, across all sensitivity attributes/exposure factors
    species.scorers <- levels(factor(species.data$Scorer))

    # Randomply assign bar colors for all the scorers for this species
    species.colors <- sample(
      all.colors[1:length(species.scorers)],
      length(species.scorers),
      replace = FALSE
    )

    for (j in 1:length(attributes)) {
      # draw a barchart on the page for each sensitivity attribute/exposure factors

      # get the data for this species and attribute
      rows <- which(
        data$Stock.Name == species[i] &
          data$Attribute.Name == attributes[j] &
          !is.na(data$Scoring.Rank1)
      )
      species.attribute.data <- as.matrix(data[rows, 7:10])
      average.data.quality <- formatC(mean(data[rows, 6]), digits = 2)

      # if there is not data for this attribute, end this interation of the loop and go on to the
      # next attribute
      if (length(rows) == 0) {
        next
      }

      # identify the scorers for this attribute and find matching colors
      species.attribute.scorers <- as.matrix(data[rows, 1])
      bar.colors <- species.colors[match(
        species.attribute.scorers,
        species.scorers
      )]

      # calculate maximum y value for barchart
      y.max <- length(species.attribute.scorers) * 5

      # plot barchart
      barplot <- barplot(
        species.attribute.data,
        main = '',
        xlab = '',
        col = bar.colors,
        ylim = c(0, y.max),
        names.arg = c('Low', 'Medium', 'High', 'Very High')
      )

      # print the name of the attribute/factor on the plot
      # also print the data quality score if that is specific as TRUE by the plot_data_quality parameter
      if (plot_data_quality == TRUE) {
        title(
          main = paste(
            attributes[j],
            "\nData quality mean=",
            average.data.quality,
            sep = ''
          ),
          adj = 0.05,
          font.main = 1,
          cex.main = 1,
          line = 0.3
        )
      } else {
        title(
          main = attributes[j],
          adj = 0.05,
          font.main = 1,
          cex.main = 1,
          line = 0.2
        )
      }
    }

    # Add the plot labels
    mtext('Ranking', 1, outer = T, line = 0.75, cex = 1.5) # x-axis label
    mtext('Score', 2, outer = T, line = 0.75, las = 0, cex = 1.5) # y-axis label
    if (preliminary == TRUE) {
      mtext('Preliminary Results', 3, outer = T, line = 4, cex = 1, col = 'red')
    }
    mtext(species[i], 3, outer = T, line = 2, cex = 1.5) # species title label

    # Add a legend of scorers along the top of plot if the plot_legend parameter equals TRUE
    if (plot_legend == TRUE) {
      legend(
        -4,
        y.max * 4.96,
        legend = species.scorers,
        fill = species.colors,
        horiz = T,
        bty = 'n',
        xpd = NA,
        cex = 1.1,
        x.intersp = 0.05,
        xjust = 0.5,
        yjust = 0
      )
    }
  } # end species barchart loop

  #########################################################
  ### Close the plotting window

  dev.off()
  graphics.off()

  if (scorer) {
    ###############################################################
    ### Make the plots for each scorer

    for (i in 1:num.scorers) {
      # got through each scorer to make plots for

      #########################################################
      ### Set up plotting window

      # open window of specified size
      #windows(w=8.5, h=11)

      # save plots as pdf with specified name
      pdf(
        paste(scorer_dir, '/', scorers[i], '.pdf', sep = ''),
        pointsize = 12,
        width = 8.5,
        height = 11
      )

      #########################################################
      ### Find all the species scored by this scorer

      scorer.rows <- which(
        data$Scorer == scorers[i] & !is.na(data$Scoring.Rank1)
      )
      scorer.data <- data[scorer.rows, ]
      species.scored <- levels(factor(scorer.data$Stock.Name))
      species.scored <- species.scored[order(match(species.scored, species))]

      #########################################################
      ### Make plots for each species scored by this scorer

      for (j in 1:length(species.scored)) {
        # Set a seed value so that we always get the same result from random sampling for this species
        set.seed(match(species.scored[j], species))

        # start new plotting page with specified margins and plotting parameters
        par(mfrow = c(4, 3), oma = c(3, 3, 6, 1), mar = c(2, 3, 2, 1), las = 1)

        # get all the rows with scores for that species
        if (sensitivity == TRUE) {
          rows <- which(
            data$Stock.Name == species.scored[j] &
              !is.na(data$Scoring.Rank1) &
              data$Attributes == 'Sensitivity Attribute'
          )
        } else {
          rows <- which(
            data$Stock.Name == species.scored[j] &
              !is.na(data$Scoring.Rank1) &
              data$Attributes == 'Exposure Factor'
          )
        }
        species.data <- data[rows, ]

        # get all the scorers who scored for this species, across sensitivity attributes
        species.scorers <- levels(factor(species.data$Scorer))

        # Randomply assign bar colors for all the scorers for this species
        species.colors <- sample(
          all.colors[1:length(species.scorers)],
          length(species.scorers),
          replace = FALSE
        )

        for (k in 1:length(attributes)) {
          # draw a barchart on the page for each sensitivity attribute/exposure factor

          # get the data for this species and attribute
          rows <- which(
            data$Stock.Name == species.scored[j] &
              data$Attribute.Name == attributes[k] &
              !is.na(data$Scoring.Rank1)
          )
          species.attribute.data <- as.matrix(data[rows, 7:10])
          average.data.quality <- formatC(mean(data[rows, 6]), digits = 2)

          # if there is not data for this attribute, end this interation of the loop and go on to the
          # next attribute
          if (length(rows) == 0) {
            next
          }

          # identify the scorers for this attribute and find matching colors
          species.attribute.scorers <- as.matrix(data[rows, 1])
          bar.colors <- species.colors[match(
            species.attribute.scorers,
            species.scorers
          )]

          # calculate maximum y value for barchart
          y.max <- length(species.attribute.scorers) * 5

          # plot barchart
          barplot <- barplot(
            species.attribute.data,
            main = '',
            xlab = '',
            col = bar.colors,
            ylim = c(0, y.max),
            names.arg = c('Low', 'Medium', 'High', 'Very High')
          )

          # print the name of the sensitivity attribute/exposure factor
          # also print the data quality score if that is specific as TRUE by the plot_data_quality parameter
          if (plot_data_quality == TRUE) {
            title(
              main = paste(
                attributes[k],
                "\nData quality mean=",
                average.data.quality,
                sep = ''
              ),
              adj = 0.05,
              font.main = 1,
              cex.main = 1,
              line = 0.3
            )
          } else {
            title(
              main = attributes[k],
              adj = 0.05,
              font.main = 1,
              cex.main = 1,
              line = 0.2
            )
          }
        }

        # Labels
        mtext('Ranking', 1, outer = T, line = 0.75, cex = 1.5) # x-axis label
        mtext('Score', 2, outer = T, line = 0.75, las = 0, cex = 1.5) # y-axis label
        if (preliminary == TRUE) {
          mtext(
            'Preliminary Results',
            3,
            outer = T,
            line = 4,
            cex = 1,
            col = 'red'
          )
        }
        mtext(species.scored[j], 3, outer = T, line = 2, cex = 1.5) # species title label

        # plot legend of scorers along top of plot
        # only plot for scorer we are giving this to
        color <- species.colors[match(scorers[i], species.scorers)]
        legend(
          -4,
          y.max * 4.96,
          legend = scorers[i],
          fill = color,
          horiz = T,
          bty = 'n',
          xpd = NA,
          cex = 1.1,
          x.intersp = 0.1,
          xjust = 0.5,
          yjust = 0
        )
      } # end species loop
    } # end scorer loop

    #########################################################
    ### Close the plotting window

    dev.off()
    graphics.off()
  } #end scorer
} # end function
