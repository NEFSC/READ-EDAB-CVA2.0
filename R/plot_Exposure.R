#' @title Make Plots of Exposure Results
#' @description
#' Produces exposure plots. Requires directory to be set up per directions in the package documentation/manual.
#'
#' @param species names of the species to plot. Must match folder name to pull correct data and save figures correctly.
#' @param type at least one of the following: 'variable', 'total', 'important', 'radar. Used to determine what to plot. Defaults to all.
#' @param presentTime,futureTime character strings indicating the present and future time series to compare. Example: '1993-2019'. Used to pull correct calculations and save the data properly
#' @param variableDF a data.frame containing all possible environmental variables, such as from the MOM6 model. Must contain columns \code{Long.Name} and \code{Short.Name}, containing the full names and abbreviated names of the variables. Abbreviated names should correspond to those in the weights vector produced by \code{combineWeights}
#' @param coastline shapefile used to plot land in model prediction plots
#'
#' @return Function does not return anything. Figures are saved to species-specific \code{figures} folder.


plot_Exposure <- function(species, type = c('variable', 'total', 'important', 'radar'), presentTime, futureTime, variableDF, coastline){

  #determine what to plot
  ind <- c('variable', 'total', 'important', 'radar') %in% type

  for(x in species){
    message(paste0('Plotting ', x , ' Exposure...'))

    if(ind[1]){
      message(paste("Plotting Variable-Specific Exposure..."))
        ###VARIABLE-LEVEL EXPOSURE
      #load variable weights
      load(paste0(file.path(getwd(),x, 'Data'), '/combined_variable_weights.RData')) #cW

      #load maps
      load(paste0(file.path(getwd(),x, 'Data'), '/', presentTime, ' vs ', futureTime, '/variable_exposure_maps.RData')) #mapExp
      #subset timeseries matrix by rownames
      i <- names(mapExp) %in% names(cW)
      mapSub <- raster::subset(mapExp, which(i == T))

      #load timeseries
      load(paste0(file.path(getwd(),x, 'Data'), '/', presentTime, ' vs ', futureTime, '/variable_exposure_timeseries.RData')) #vecExp
      #subset timeseries matrix by rownames
      i <- rownames(vecExp) %in% names(cW)
      vecSub <- vecExp[i,]

        #plot
        pdf(paste0(file.path(getwd(),x, 'Figures'), '/', presentTime, ' vs ', futureTime, '/variable_exposure_maps_inset_timeseries.pdf'), width = 8, height = 11)
        #set up panels according to the number of variables
        if(raster::nlayers(mapSub) < 6){
          par(mfrow=c(2,3))
        } else {
          par(mfrow=c(3,3))
        }

        for(y in 1:raster::nlayers(mapSub)){
          #get full name of variable
          i <- variableDF$Short.Name %in% names(mapSub)[y]

          #map
          par(plt = c(0.2, 0.9, 0.15, 0.875))
          plot(raster::subset(mapSub, y), zlim = c(1,4), col = cmocean::cmocean('matter')(4), legend = F, legend.mar = 0, xlab = expression('Longitude ('*degree*')'), ylab = expression('Latitude ('*degree*')'), xaxt = 'n', yaxt = 'n', main = variableDF$Long.Name[i])
          axis(2, at = seq(30, 50, by = 1), labels = seq(30, 50, by = 1), las = 2)
          axis(1, at = seq(-85, -65, by = 1), labels = seq(-85, -65, by = 1))
          plot(coastline['id'], col = 'grey', add = T)

          #inset timeseries
          par(plt = c(0.55, 0.9, 0.25, 0.45), new = TRUE)
          plot(vecSub[y,], t = 'b', lty = 8, lwd = 0.8, cex = 0.8, pch = y, ylim = c(1, 4), ylab = "", xlab = "", yaxt = 'n', xaxt = 'n')
          axis(1, at = 1:12, labels = month.abb, las = 2)
          axis(2, at = 1:4, labels = c('L', "M", "H", "VH"), las = 2, cex.lab = 0.75)
        }

        if(raster::nlayers(mapSub) != 6){ #add legend on the last one if the number of variables is not 6
          plot(1:10, t = 'n', axes = F, xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
          fields::image.plot(matrix(seq(1,4,length.out = 16), 4,4), legend.only = T, horizontal = F, legend.shrink = 0.7,
                             smallplot = c(0.4, 0.6, 0.2, 0.8),
                             legend.args = list(text = 'Exposure', cex = 1.25, side = 3, line = 0.1),
                             axis.args = list(cex.axis =1, at = 1:4, labels = c('Low (L)', "Moderate (M)", "High (H)", "Very High (VH)"), mgp = c(3, 0.5, 0)), col = cmocean::cmocean('matter')(4))
        } else { #if the number of variables is 6, it will still be a 3x3 grid, so put legend in the middle by adding an extra plot
          plot(1:10, t = 'n', axes = F, xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
          plot(1:10, t = 'n', axes = F, xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
          fields::image.plot(matrix(seq(1,4,length.out = 16), 4,4), legend.only = T, horizontal = F, legend.shrink = 0.7,
                             smallplot = c(0.4, 0.6, 0.2, 0.8),
                             legend.args = list(text = 'Exposure', cex = 1.25, side = 3, line = 0.1),
                             axis.args = list(cex.axis =1, at = 1:4, labels = c('Low (L)', "Moderate (M)", "High (H)", "Very High (VH)"), mgp = c(3, 0.5, 0)), col = cmocean::cmocean('matter')(4))
        }

        dev.off()
  }

    if(ind[2]){
      message(paste("Plotting Total Exposure with All Variables..."))
      ##TOTAL EXPOSURE - ALL VARIABLES
      #load total map
      load(paste0(file.path(getwd(),x, 'Data'), '/', presentTime, ' vs ', futureTime, '/total_exposure_maps_all.RData')) #totalM

      #load total timeseries
      load(paste0(file.path(getwd(),x, 'Data'), '/', presentTime, ' vs ', futureTime, '/total_exposure_timeseries_all.RData')) #totalT

      #plot
      pdf(paste0(file.path(getwd(),x, 'Figures'), '/', presentTime, ' vs ', futureTime, '/total_exposure_maps_inset_timeseries_allvars.pdf'), width = 8, height = 11)
      #map
      par(fig = c(0, 1, 0, 1))
      plot(totalM, zlim = c(1,4), col = cmocean::cmocean('matter')(4), ylim = c(35,45), legend = F, xlab = expression('Longitude ('*degree*')'), ylab = expression('Latitude ('*degree*')'), xaxt = 'n', yaxt = 'n', legend.mar = 0)
      axis(2, at = seq(30, 50, by = 1), labels = seq(30, 50, by = 1), las = 2)
      axis(1, at = seq(-85, -65, by = 1), labels = seq(-85, -65, by = 1))
      plot(coastline['id'], col = 'grey', add = T)
      fields::image.plot(matrix(seq(1,4,length.out = 16), 4,4), legend.only = T, horizontal = T, legend.shrink = 0.7,
                         smallplot = c(0.5, 0.9, 0.15, 0.2),
                         legend.args = list(text = 'Exposure', cex = 1.5, side = 3, line = 0.1),
                         axis.args = list(cex.axis =1, at = 1:4, labels = c('Low', "Moderate", "High", "Very High"), mgp = c(3, 0.5, 0)), col = cmocean::cmocean('matter')(4))

      par(fig = c(0.125, 0.6, 0.65, 0.95), new = TRUE)
      plot(totalT, t = 'b', lty = 8, lwd = 1.5, pch = 19, ylim = c(1, 4), ylab = "", xlab = "Month", yaxt = 'n', xaxt = 'n')
      axis(1, at = 1:12, labels = month.abb, las = 2)
      axis(2, at = 1:4, labels = c('Low', "Moderate", "High", "Very High"), las = 2, cex.lab = 0.75)
      dev.off()
  }

    if(ind[3]){
      message(paste("Plotting Total Exposure with Important Variables..."))
      ### ONLY IMPORTANT VARS
      #load total map
      load(paste0(file.path(getwd(),x, 'Data'), '/', presentTime, ' vs ', futureTime, '/total_exposure_maps_subset.RData')) #totalM

      #load total timeseries
      load(paste0(file.path(getwd(),x, 'Data'), '/', presentTime, ' vs ', futureTime, '/total_exposure_timeseries_subset.RData')) #totalT

      #plot
      pdf(paste0(file.path(getwd(),x, 'Figures'), '/', presentTime, ' vs ', futureTime, '/total_exposure_maps_inset_timeseries_impvars.pdf'), width = 8, height = 11)
      #map
      par(fig = c(0, 1, 0, 1))
      plot(totalM, zlim = c(1,4), col = cmocean::cmocean('matter')(4), ylim = c(35,45), legend = F, xlab = expression('Longitude ('*degree*')'), ylab = expression('Latitude ('*degree*')'), xaxt = 'n', yaxt = 'n', legend.mar = 0)
      axis(2, at = seq(30, 50, by = 1), labels = seq(30, 50, by = 1), las = 2)
      axis(1, at = seq(-85, -65, by = 1), labels = seq(-85, -65, by = 1))
      plot(coastline['id'], col = 'grey', add = T)
      fields::image.plot(matrix(seq(1,4,length.out = 16), 4,4), legend.only = T, horizontal = T, legend.shrink = 0.7,
                         smallplot = c(0.5, 0.9, 0.15, 0.2),
                         legend.args = list(text = 'Exposure', cex = 1.5, side = 3, line = 0.1),
                         axis.args = list(cex.axis =1, at = 1:4, labels = c('Low', "Moderate", "High", "Very High"), mgp = c(3, 0.5, 0)), col = cmocean::cmocean('matter')(4))

      par(fig = c(0.125, 0.6, 0.65, 0.95), new = TRUE)
      plot(totalT, t = 'b', lty = 8, lwd = 1.5, pch = 19, ylim = c(1, 4), ylab = "", xlab = "Month", yaxt = 'n', xaxt = 'n')
      axis(1, at = 1:12, labels = month.abb, las = 2)
      axis(2, at = 1:4, labels = c('Low', "Moderate", "High", "Very High"), las = 2, cex.lab = 0.75)
      dev.off()
    }

    if(ind[4]){
      message(paste("Plotting Radar Plot of Relative Variable Importance..."))
      #load variable weights
      load(paste0(file.path(getwd(),x, 'Data'), '/combined_variable_weights.RData')) #cW

      cW <- rbind(rep(1, length(cW)), rep(0, length(cW)), cW)

      #plot
      pdf(paste0(file.path(getwd(),x, 'Figures'), '/dynamic_variable_weights.pdf'), width = 8, height = 8)
      fmsb::radarchart(as.data.frame(cW), pfcol = scales::alpha('grey', 0.5))
      dev.off()

    }

  } #end x

}
