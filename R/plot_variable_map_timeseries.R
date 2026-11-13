#' @title Make Variable Vulnerability Plots
#' @description
#' Produces variable Vulnerability maps and time series
#'
#' @param data either a multi-layer spatRaster or a matrix containing the map data or timeseries data to plot
#' @param stocks spatVector of stock polygons to add to plot. Defaults to NULL meaning no stocks are available.
#' @param stock_key named vector of stock long names and abbreviations. Used to label timeseries.
#' @param variable_df a data.frame containing Short.Names and Long.Names. Used to isolate variables and label plots.
#' @param fig_name name to save figure as. Figure will save to current working directory if desired directory is not included in the figure name
#' @param coastline shapefile used to plot land in model prediction plots
#' @param bathymetry spatRaster file of bathymetry data; used to plot bathymetry in stock boundary plots

#' @return Function does not return anything. Figures are saved to species-specific \code{figures} folder.
#'
#'@export
#'
plot_variable_map_timeseries <- function(
  data,
  stocks = NULL,
  stock_key,
  variable_df,
  fig_name,
  coastline,
  bathymetry
) {
  # Failsafe: If the function crashes, forcefully close any open PDFs
  on.exit(
    while (grDevices::dev.cur() > 1) {
      grDevices::dev.off()
    },
    add = TRUE
  )

  plot_variable_exp <- function() {
    if (inherits(data, 'SpatRaster')) {
      message('Plotting Variable Exposure Maps')
      # Use dynamic rows and set global outer margins (oma)
      graphics::par(
        mfrow = c(n_rows, 3),
        oma = c(1, 1, 1, 1)
      )

      for (y in 1:terra::nlyr(data)) {
        #get full name of variable
        i <- variable_df$Short.Name %in% names(data)[y]

        #draw_legend <- y == terra::nlyr(data)

        #map
        #graphics::par(plt = c(0.1, 0.98, 0.1, 0.95))
        terra::plot(
          data[[y]],
          type = 'continuous',
          range = c(1, 4),
          col = cmocean::cmocean('deep')(64),
          mar = c(2.5, 2.5, 1.5, 0.5), # Explicitly set margins inside terra::plot
          legend = TRUE,
          xlab = expression('Longitude (' * degree * ')'),
          ylab = expression('Latitude (' * degree * ')'),
          main = variable_df$Long.Name[i],
          pax = list(cex.axis = 1.5, xat = seq(-80, -60, by = 2)),
          cex.lab = 1.25,
          plg = list(
            title = "Exposure",
            title.cex = 1.5,
            cex = 1.5,
            horizontal = TRUE,
            x = -73.5,
            y = 37,
            at = 1:4,
            n = 4,
            # 1. Scale the size of the color bar itself (width, height)
            size = c(1, 2.5)
          )
        )

        #add bathy contours, coastline, and stocks if necessary
        terra::contour(
          bathymetry,
          filled = F,
          levels = c(-1000, -100, -50),
          add = T
        )
        terra::plot(coastline['id'], col = 'grey', add = T)
        if (!is.null(stocks)) {
          terra::plot(stocks, add = T, lwd = 2)
        }
      } #end y
      # 2. Draw the legend independently if it is the last panel
      #   if (draw_legend) {
      #     graphics::par(mgp = c(3, 0.1, 0))
      #     terra::plot(
      #      data[[y]],
      #     type = 'continuous',
      #    range = c(1, 4),
      #   col = cmocean::cmocean('deep')(64),
      #  mar = c(2.5, 2.5, 1.5, 0.5), # Explicitly set margins inside terra::plot
      # legend.only = TRUE, # <-- Draws only the legend elements
      #        plg = list(
      #         title = "Exposure",
      #        title.cex = 1.5,
      #       cex = 1.5,
      #      horizontal = TRUE,
      #     x = -73.5,
      #    y = 37,
      #   at = 1:4,
      #  n = 4,
      # 1. Scale the size of the color bar itself (width, height)
      #  size = c(1, 2.5)
      #)
      #)
      #}
    } else {
      #end if spatRaster
      #if data is not a spatraster, then its a matrix or list of matrices
      message('Plotting Variable Exposure Timeseries')
      if (!inherits(data, 'list')) {
        #if vecExp is NOT a list and is just a single matrix, just plot a single line

        #set up panels according to the number of variables
        if (nrow(data) <= 6) {
          graphics::par(mfrow = c(2, 3), mar = c(4, 3, 2, 2))
        } else {
          graphics::par(mfrow = c(3, 3), mar = c(4, 3, 2, 2))
        }

        for (y in 1:nrow(data)) {
          #get full name of variable
          i <- variable_df$Short.Name %in% rownames(data)[y]

          plot(
            data[y, ],
            t = 'b',
            lty = 1,
            lwd = 1,
            cex = 1,
            pch = 1,
            ylim = c(1, 4),
            ylab = "Exposure",
            xlab = "Month",
            yaxt = 'n',
            xaxt = 'n',
            main = variable_df$Long.Name[i]
          )
          graphics::axis(
            1,
            at = 1:12,
            labels = month.abb,
            las = 2,
            cex.lab = 0.5
          )
          graphics::axis(
            2,
            at = 1:4,
            labels = c('L', "M", "H", "VH"),
            las = 2,
            cex.lab = 1.25
          )
        }
      } else {
        #if it is a list, then exposure is calculated within multiple stocks

        #set up panels according to the number of variables
        if (nrow(data[[1]]) <= 6) {
          graphics::par(mfrow = c(2, 3), mar = c(4, 3, 2, 2))
        } else {
          graphics::par(mfrow = c(3, 3), mar = c(4, 3, 2, 2))
        }

        for (y in 1:nrow(data[[1]])) {
          i <- variable_df$Short.Name %in% rownames(data[[1]])[y]
          vecSub <- do.call(rbind, lapply(data, function(s) s[y, ]))
          plot(
            #explicitly call the first row, which is the global value and then add the additional ones
            vecSub[1, ],
            t = 'b',
            lty = 1,
            lwd = 1,
            cex = 1,
            pch = 1,
            ylim = c(1, 4),
            ylab = "",
            xlab = "",
            yaxt = 'n',
            xaxt = 'n',
            main = variable_df$Long.Name[i]
          )
          for (m in 2:nrow(vecSub)) {
            graphics::lines(
              vecSub[m, ],
              t = 'b',
              lty = m,
              lwd = 1,
              cex = 1,
              pch = m
            )
          } #end m

          graphics::axis(
            1,
            at = 1:12,
            labels = month.abb,
            las = 2,
            cex.lab = 0.5
          )
          graphics::axis(
            2,
            at = 1:4,
            labels = c('L', "M", "H", "VH"),
            las = 2,
            cex.lab = 1.25
          )
        } #end y
        #add legend
        graphics::legend(
          'top',
          bty = 'n',
          legend = c('All', names(data)[-1]),
          pch = 1:length(data),
          lty = 1:length(data),
          cex = 1,
          title = 'Stocks',
          ncol = 2
        )
      } #end if vecExp is a list
    } #end if not spatraster
  }

  # Determine number of panels and rows
  if (inherits(data, 'SpatRaster')) {
    n_panels <- terra::nlyr(data)
  } else {
    n_panels <- if (inherits(data, 'list')) nrow(data[[1]]) else nrow(data)
  }

  n_rows <- ceiling(n_panels / 3)

  # Calculate dynamic height (approx 3.5 inches per row to match width)
  dynamic_height <- n_rows * 3.5

  grDevices::pdf(
    paste0(fig_name, '.pdf'),
    width = 8,
    height = dynamic_height # Updated from 11
  )
  plot_variable_exp()
  grDevices::dev.off()

  grDevices::png(
    paste0(fig_name, '.png'),
    width = 8,
    height = dynamic_height, # Updated from 11
    units = 'in',
    res = 300
  )
  plot_variable_exp()
  grDevices::dev.off()
}
