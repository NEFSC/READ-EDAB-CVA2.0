#' @title Make Plots of Exposure Results
#' @description
#' Produces exposure plots. Requires directory to be set up per directions in the package documentation/manual.
#'
#' @param species names of the species to plot. Must match folder name to pull correct data and save figures correctly.
#' @param type at least one of the following: 'variable', 'total', 'important', 'radar'. Used to determine what to plot. Defaults to all.
#' @param forecast_release,hindcast_release MOM6 release codes for the (f)orecast and (h)indcasts used. Used to pull correct variable exposures
#' @param forecast_init forecast_initialization code corresponding to the forecast_initalization date of the desired forecast data. Used to pull correct variable exposures
#' @param hindcast_yr_range character string corresponding to the years in the hindcast data used. Used to pull correct data and save the data properly
#' @param variable_df a data.frame containing all possible environmental variables, such as from the MOM6 model. Must contain columns \code{Long.Name} and \code{Short.Name}, containing the full names and abbreviated names of the variables. Abbreviated names should correspond to those in the weights vector produced by \code{combineWeights}
#' @param coastline shapefile used to plot land in model prediction plots
#'
#' @return Function does not return anything. Figures are saved to species-specific \code{figures} folder.
#'
#'@export

make_exposure_plots <- function(
  species,
  type = c('variable', 'total', 'important', 'radar'),
  forecast_release, 
  forecast_init,
  hindcast_release, 
  hindcast_yr_range,
  variable_df,
  coastline
) {
  #determine what to plot
  ind <- c('variable', 'total', 'important', 'radar') %in% type

  for (x in species) {
    message(paste0('Plotting ', x, ' Exposure...'))

    if (ind[1]) {
      message(paste("Plotting Variable-Specific Exposure..."))
      ###VARIABLE-LEVEL EXPOSURE
      #load variable weights
      imp <- readRDS(paste0(
        file.path(getwd(), x, 'Data'),
        '/normalized_dynamic_variable_importance.rds'
      )) 

      #load maps
      varMaps <- terra::rast(file.path(getwd(), x, 'Data',
      paste0('variable_exposure_maps_', forecast_release, '_', forecast_init, '_', hindcast_release,'_',hindcast_yr_range, '.tif')))

      #load timeseries
      vecExp <- readRDS(
        file.path(getwd(), x, 'Data',
        paste0('variable_exposure_timeseries_', forecast_release, '_', forecast_init, '_', hindcast_release,'_',hindcast_yr_range,'.rds')
      ))

      #plot
      grDevices::pdf(
        paste0(
          file.path(getwd(), x, 'Figures',
          paste0('variable_exposure_maps_inset_timeseries',forecast_release, '_', forecast_init, '_', hindcast_release,'_',hindcast_yr_range,'.pdf'))
        ),
        width = 8,
        height = 11
      )
      #set up panels according to the number of variables
      if (terra::nlyr(varMaps) < 6) {
        graphics::par(mfrow = c(2, 3), mar = c(2, 2.5, 1.5, 0.5), mgp = c(1.2, 0.5, 0))
      } else {
        graphics::par(mfrow = c(3, 3), mar = c(2, 2.5, 1.5, 0.5), mgp = c(1.2, 0.5, 0))
      }

      for (y in 1:terra::nlyr(varMaps)) {
        #get full name of variable
        i <- variable_df$Short.Name %in% names(varMaps)[y]

        #map
        graphics::par(plt = c(0.1, 0.98, 0.1, 0.95))
       terra::plot(
          varMaps[[y]],
          zlim = c(1, 4),
          col = cmocean::cmocean('matter')(4),
          legend = F,
          legend.mar = 0,
          xlab = expression('Longitude (' * degree * ')'),
          ylab = expression('Latitude (' * degree * ')'),
          pax = list(cex = 1.3, xat = seq(-80, -60, by = 2)),
          main = variable_df$Long.Name[i],
          cex.main = 1.5
        )
        plot(coastline['id'], col = 'grey', add = T)
        plot(stocks, add = T)

        #inset timeseries
        graphics::par(plt = c(0.55, 0.9, 0.3, 0.5), new = TRUE)
        if(!inherits(vecExp, 'list')){ #if vecExp is NOT a list and is just a single matrix, just plot a single line 
          plot(
            vecExp[y, ],
            t = 'b',
            lty = 1,
            lwd = 0.8,
            cex = 0.8,
            pch = 1,
            ylim = c(1, 4),
            ylab = "",
            xlab = "",
            yaxt = 'n',
            xaxt = 'n'
          )
        } else { #if vecExp is a list, then exposure is calculated within multiple stocks 
          vecSub <- do.call(rbind, lapply(vecExp, function(s) s[y,]))
          plot( #explicitly call the first row, which is the global value and then add the additional ones 
            vecSub[1, ],
            t = 'b',
            lty = 1,
            lwd = 0.8,
            cex = 0.8,
            pch = 1,
            ylim = c(1, 4),
            ylab = "",
            xlab = "",
            yaxt = 'n',
            xaxt = 'n'
          )
          for(m in 2:nrow(vecSub)){
            lines(vecSub[m,], 
                  t = 'b', 
                  lty = m,
                  lwd = 0.8,
                  cex = 0.8,
                  pch = m)
          } #end m
        } #end if vecExp is a list
        graphics::axis(1, at = 1:12, labels = month.abb, las = 2, cex.lab = 0.5)
        graphics::axis(
          2,
          at = 1:4,
          labels = c('L', "M", "H", "VH"),
          las = 2,
          cex.lab = 0.5
        )
      } #end y 

      if (terra::nlyr(varMaps) != 6) {
        #add legend on the last one if the number of variables is not 6
        plot(
          1:10,
          t = 'n',
          axes = F,
          xaxt = 'n',
          yaxt = 'n',
          xlab = '',
          ylab = ''
        )
        fields::image.plot(
          matrix(seq(1, 4, length.out = 16), 4, 4),
          legend.only = T,
          horizontal = F,
          legend.shrink = 0.7,
          smallplot = c(0.4, 0.6, 0.2, 0.8),
          legend.args = list(
            text = 'Exposure',
            cex = 1.25,
            side = 3,
            line = 0.1
          ),
          axis.args = list(
            cex.axis = 1,
            at = 1:4,
            labels = c('Low (L)', "Moderate (M)", "High (H)", "Very High (VH)"),
            mgp = c(3, 0.5, 0)
          ),
          col = cmocean::cmocean('matter')(4)
        )
      } else {
        #if the number of variables is 6, it will still be a 3x3 grid, so put legend in the middle by adding an extra plot
        plot(
          1:10,
          t = 'n',
          axes = F,
          xaxt = 'n',
          yaxt = 'n',
          xlab = '',
          ylab = ''
        )
        plot(
          1:10,
          t = 'n',
          axes = F,
          xaxt = 'n',
          yaxt = 'n',
          xlab = '',
          ylab = ''
        )
        fields::image.plot(
          matrix(seq(1, 4, length.out = 16), 4, 4),
          legend.only = T,
          horizontal = F,
          legend.shrink = 0.7,
          smallplot = c(0.4, 0.6, 0.2, 0.8),
          legend.args = list(
            text = 'Exposure',
            cex = 1.25,
            side = 3,
            line = 0.1
          ),
          axis.args = list(
            cex.axis = 1,
            at = 1:4,
            labels = c('Low (L)', "Moderate (M)", "High (H)", "Very High (VH)"),
            mgp = c(3, 0.5, 0)
          ),
          col = cmocean::cmocean('matter')(4)
        )
      }

      grDevices::dev.off()
    }

    if (ind[2]) {
      message(paste("Plotting Total Exposure with All Variables..."))
      ##TOTAL EXPOSURE - ALL VARIABLES
      #load total map
      load(paste0(
        file.path(getwd(), x, 'Data'),
        '/',
        present_time,
        ' vs ',
        future_time,
        '/total_exposure_maps_all.RData'
      )) #totalM

      #load total timeseries
      load(paste0(
        file.path(getwd(), x, 'Data'),
        '/',
        present_time,
        ' vs ',
        future_time,
        '/total_exposure_timeseries_all.RData'
      )) #totalT

      #plot
      grDevices::pdf(
        paste0(
          file.path(getwd(), x, 'Figures'),
          '/',
          present_time,
          ' vs ',
          future_time,
          '/total_exposure_maps_inset_timeseries_allvars.pdf'
        ),
        width = 8,
        height = 11
      )
      #map
      graphics::par(fig = c(0, 1, 0, 1))
      plot(
        totalM,
        zlim = c(1, 4),
        col = cmocean::cmocean('matter')(4),
        ylim = c(35, 45),
        legend = F,
        xlab = expression('Longitude (' * degree * ')'),
        ylab = expression('Latitude (' * degree * ')'),
        xaxt = 'n',
        yaxt = 'n',
        legend.mar = 0
      )
      graphics::axis(
        2,
        at = seq(30, 50, by = 1),
        labels = seq(30, 50, by = 1),
        las = 2
      )
      graphics::axis(
        1,
        at = seq(-85, -65, by = 1),
        labels = seq(-85, -65, by = 1)
      )
      plot(coastline['id'], col = 'grey', add = T)
      fields::image.plot(
        matrix(seq(1, 4, length.out = 16), 4, 4),
        legend.only = T,
        horizontal = T,
        legend.shrink = 0.7,
        smallplot = c(0.5, 0.9, 0.15, 0.2),
        legend.args = list(text = 'Exposure', cex = 1.5, side = 3, line = 0.1),
        axis.args = list(
          cex.axis = 1,
          at = 1:4,
          labels = c('Low', "Moderate", "High", "Very High"),
          mgp = c(3, 0.5, 0)
        ),
        col = cmocean::cmocean('matter')(4)
      )

      graphics::par(fig = c(0.125, 0.6, 0.65, 0.95), new = TRUE)
      plot(
        totalT,
        t = 'b',
        lty = 8,
        lwd = 1.5,
        pch = 19,
        ylim = c(1, 4),
        ylab = "",
        xlab = "Month",
        yaxt = 'n',
        xaxt = 'n'
      )
      graphics::axis(1, at = 1:12, labels = month.abb, las = 2)
      graphics::axis(
        2,
        at = 1:4,
        labels = c('Low', "Moderate", "High", "Very High"),
        las = 2,
        cex.lab = 0.75
      )
      grDevices::dev.off()
    }

    if (ind[3]) {
      message(paste("Plotting Total Exposure with Important Variables..."))
      ### ONLY IMPORTANT VARS
      #load total map
      load(paste0(
        file.path(getwd(), x, 'Data'),
        '/',
        present_time,
        ' vs ',
        future_time,
        '/total_exposure_maps_subset.RData'
      )) #totalM

      #load total timeseries
      load(paste0(
        file.path(getwd(), x, 'Data'),
        '/',
        present_time,
        ' vs ',
        future_time,
        '/total_exposure_timeseries_subset.RData'
      )) #totalT

      #plot
      grDevices::pdf(
        paste0(
          file.path(getwd(), x, 'Figures'),
          '/',
          present_time,
          ' vs ',
          future_time,
          '/total_exposure_maps_inset_timeseries_impvars.pdf'
        ),
        width = 8,
        height = 11
      )
      #map
      graphics::par(fig = c(0, 1, 0, 1))
      plot(
        totalM,
        zlim = c(1, 4),
        col = cmocean::cmocean('matter')(4),
        ylim = c(35, 45),
        legend = F,
        xlab = expression('Longitude (' * degree * ')'),
        ylab = expression('Latitude (' * degree * ')'),
        xaxt = 'n',
        yaxt = 'n',
        legend.mar = 0
      )
      graphics::axis(
        2,
        at = seq(30, 50, by = 1),
        labels = seq(30, 50, by = 1),
        las = 2
      )
      graphics::axis(
        1,
        at = seq(-85, -65, by = 1),
        labels = seq(-85, -65, by = 1)
      )
      plot(coastline['id'], col = 'grey', add = T)
      fields::image.plot(
        matrix(seq(1, 4, length.out = 16), 4, 4),
        legend.only = T,
        horizontal = T,
        legend.shrink = 0.7,
        smallplot = c(0.5, 0.9, 0.15, 0.2),
        legend.args = list(text = 'Exposure', cex = 1.5, side = 3, line = 0.1),
        axis.args = list(
          cex.axis = 1,
          at = 1:4,
          labels = c('Low', "Moderate", "High", "Very High"),
          mgp = c(3, 0.5, 0)
        ),
        col = cmocean::cmocean('matter')(4)
      )

      graphics::par(fig = c(0.125, 0.6, 0.65, 0.95), new = TRUE)
      plot(
        totalT,
        t = 'b',
        lty = 8,
        lwd = 1.5,
        pch = 19,
        ylim = c(1, 4),
        ylab = "",
        xlab = "Month",
        yaxt = 'n',
        xaxt = 'n'
      )
      graphics::axis(1, at = 1:12, labels = month.abb, las = 2)
      graphics::axis(
        2,
        at = 1:4,
        labels = c('Low', "Moderate", "High", "Very High"),
        las = 2,
        cex.lab = 0.75
      )
      grDevices::dev.off()
    }

    if (ind[4]) {
      message(paste("Plotting Radar Plot of Relative Variable Importance..."))
      #load variable weights
      load(paste0(
        file.path(getwd(), x, 'Data'),
        '/combined_variable_weights.RData'
      )) #cW

      cW <- rbind(rep(1, length(cW)), rep(0, length(cW)), cW)

      #plot
      grDevices::pdf(
        paste0(
          file.path(getwd(), x, 'Figures'),
          '/dynamic_variable_weights.pdf'
        ),
        width = 8,
        height = 8
      )
      fmsb::radarchart(as.data.frame(cW), pfcol = scales::alpha('grey', 0.5))
      grDevices::dev.off()
    }
  } #end x
}
