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
  coastline, 
  bathymetry
) {
  # Failsafe: If the function crashes, forcefully close any open PDFs
  on.exit(while (grDevices::dev.cur() > 1) grDevices::dev.off(), add = TRUE)
  
  #determine what to plot
  ind <- c('variable', 'total', 'important', 'radar') %in% type

  for (x in species) {
    message(paste0('Plotting ', x, ' Exposure...'))
    
    if(file.exists(paste0('../shpfiles/species_stock_areas/', x, '.shp'))){
      stocks <- terra::vect(paste0('../shpfiles/species_stock_areas/', x, '.shp'))
    } else {
      stocks <- NULL
    }
    
    if (ind[1]) {
      message(paste("Plotting Variable-Specific Exposure..."))
      ###VARIABLE-LEVEL EXPOSURE

      #load maps
      varMaps <- terra::rast(file.path(getwd(), x, 'Data',
      paste0('variable_exposure_maps_', forecast_release, '_', forecast_init, '_', hindcast_release,'_',hindcast_yr_range, '.tif')))

      #load timeseries
      vecExp <- readRDS(
        file.path(getwd(), x, 'Data',
        paste0('variable_exposure_timeseries_', forecast_release, '_', forecast_init, '_', hindcast_release,'_',hindcast_yr_range,'.rds')
      ))

      #plot maps
      grDevices::pdf(
        paste0(
          file.path(getwd(), x, 'Figures',
          paste0('variable_exposure_maps_',forecast_release, '_', forecast_init, '_', hindcast_release,'_',hindcast_yr_range,'.pdf'))
        ),
        width = 11,
        height = 8
      )
      #set up panels according to the number of variables
      if (terra::nlyr(varMaps) <= 6) {
        graphics::par(mfrow = c(2, 3), mar = c(2, 3.5, 1.5, 0.5), mgp = c(1.2, 0.5, 0))
      } else {
        graphics::par(mfrow = c(3, 3), mar = c(2, 3.5, 1.5, 0.5), mgp = c(1.2, 0.5, 0))
      }

      for (y in 1:terra::nlyr(varMaps)) {
        #get full name of variable
        i <- variable_df$Short.Name %in% names(varMaps)[y]
        
        draw_legend <- y == terra::nlyr(varMaps)

        #map
        #graphics::par(plt = c(0.1, 0.98, 0.1, 0.95))
       terra::plot(
          varMaps[[y]],
          type = 'continuous',
          range = c(1, 4),
          col = cmocean::cmocean('matter')(64),
          legend = FALSE,
          xlab = expression('Longitude (' * degree * ')'),
          ylab = expression('Latitude (' * degree * ')'),
          main = variable_df$Long.Name[i],
          pax = list(cex.axis = 1.5, xat = seq(-80, -60, by = 2)), cex.lab = 1.25
       )
        
       #add bathy contours, coastline, and stocks if necessary
       terra::contour(bathymetry, filled = F, levels = c(-1000, -100, -50), add = T)
       plot(coastline['id'], col = 'grey', add = T)
       if(!is.null(stocks)){
         terra::plot(stocks, add = T, lwd = 2)
       }
        
      } #end y 
      # 2. Draw the legend independently if it is the last panel
      if (draw_legend) {
        terra::plot(
          varMaps[[y]],
          type = 'continuous',
          range = c(1, 4),
          col = cmocean::cmocean('matter')(64),
          legend.only = TRUE, # <-- Draws only the legend elements
          plg = list(
            title = "Exposure",
            title.cex = 1.5,
            cex = 1.5,
            horizontal = TRUE,
            x = -73.5,
            y = 36.5,
            at = 1:4, n = 4,
            # 1. Scale the size of the color bar itself (width, height)
            size = c(1, 2.5), 
            # 2. Control the distance of the labels from the bar
            pax = list(
              mgp = c(3, 2, 0) 
            )
          )
        )
      }

      grDevices::dev.off()
      
      #timeseries
      grDevices::pdf(
        paste0(
          file.path(getwd(), x, 'Figures',
                    paste0('variable_exposure_timeseries_',forecast_release, '_', forecast_init, '_', hindcast_release,'_',hindcast_yr_range,'.pdf'))
        ),
        width = 8,
        height = 11
      )
      #set up panels according to the number of variables
      if (terra::nlyr(varMaps) <= 6) {
        graphics::par(mfrow = c(2, 3), mar = c(4,3,2,2))
      } else {
        graphics::par(mfrow = c(3, 3), mar = c(4,3,2,2))
      }
      
      for (y in 1:terra::nlyr(varMaps)) {
        #get full name of variable
        i <- variable_df$Short.Name %in% names(varMaps)[y]
        
        if(!inherits(vecExp, 'list')){ #if vecExp is NOT a list and is just a single matrix, just plot a single line 
          plot(
            vecExp[y, ],
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
        } else { #if vecExp is a list, then exposure is calculated within multiple stocks 
          vecSub <- do.call(rbind, lapply(vecExp, function(s) s[y,]))
          plot( #explicitly call the first row, which is the global value and then add the additional ones 
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
          for(m in 2:nrow(vecSub)){
            lines(vecSub[m,], 
                  t = 'b', 
                  lty = m,
                  lwd = 1,
                  cex = 1,
                  pch = m)
          } #end m
        } #end if vecExp is a list
        graphics::axis(1, at = 1:12, labels = month.abb, las = 2, cex.lab = 0.5)
        graphics::axis(
          2,
          at = 1:4,
          labels = c('L', "M", "H", "VH"),
          las = 2,
          cex.lab = 1.25
        )
      }
      
        if(inherits(vecExp, 'list')){ #add legend if necessary
          legend('top', bty = 'n', legend = c('All', names(vecExp)[-1]), pch = 1:length(vecExp), lty = 1:length(vecExp), cex = 1, title = 'Stocks', ncol = 2)
        }
      
      grDevices::dev.off()
      
    }

    if (ind[2]) {
      message(paste("Plotting Total Exposure with All Variables..."))
      ##TOTAL EXPOSURE - ALL VARIABLES
      #load maps
      varMaps <- terra::rast(file.path(getwd(), x, 'Data',
                                       paste0('total_exposure_map_all_var_', forecast_release, '_', forecast_init, '_', hindcast_release,'_',hindcast_yr_range, '.tif')))
      
      #load timeseries
      vecExp <- readRDS(
        file.path(getwd(), x, 'Data',
                  paste0('total_exposure_timeseries_all_var_', forecast_release, '_', forecast_init, '_', hindcast_release,'_',hindcast_yr_range,'.rds')
        ))
      if(!is.null(nrow(vecExp))){
        rownames(vecExp)[1] <- 'All Stocks'
      }
      
      #plot map
      grDevices::pdf(
        paste0(
          file.path(getwd(), x, 'Figures',
                    paste0('total_exposure_all_var_map_',forecast_release, '_', forecast_init, '_', hindcast_release,'_',hindcast_yr_range,'.pdf'))
        ),
        width = 8,
        height = 8
      )
      #map
      terra::plot(
        varMaps,
        type = 'classes',
        levels = c("1", "2", "3", "4"),
        range = c(1, 4),
        col = cmocean::cmocean('matter')(4),
        ylim = c(35, 45),
        legend = F,
        pax = list(cex.axis = 1.5), cex.lab = 1.25,
        xlab = expression('Longitude (' * degree * ')'),
        ylab = expression('Latitude (' * degree * ')'),
        mar = c(3,3,1.5,0.5)
      )

      #add bathy contours, coastline, and stocks if necessary
      terra::contour(bathymetry, filled = F, levels = c(-1000, -100, -50), add = T)
      plot(coastline['id'], col = 'grey', add = T)
      if(!is.null(stocks)){
        terra::plot(stocks, add = T, lwd = 2)
      }
      
    #legend
      terra::plot(
        varMaps,
        type = 'classes',
        levels = c("1", "2", "3", "4"),
        range = c(1, 4),
        col = cmocean::cmocean('matter')(4),
        legend.only = TRUE, # <-- Draws only the legend elements
        plg = list(
          title = "Exposure",
          title.cex = 1.5,
          cex = 1.5,
          x = -68,
          y = 38, 
          legend = c("Low", "Moderate", "High", "Very High")
        )
      )
      grDevices::dev.off()
      
      #plot timeseries
      grDevices::pdf(
        paste0(
          file.path(getwd(), x, 'Figures',
                    paste0('total_exposure_all_vars_timeseries_',forecast_release, '_', forecast_init, '_', hindcast_release,'_',hindcast_yr_range,'.pdf'))
        ),
        width = 11,
        height = 8
      )
      graphics::par(mar = c(5,5,2,2))
      
      if(!inherits(vecExp, 'matrix')){
        plot(
          vecExp,
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
          labels = c('Low', "Moderate", "High", "Very\nHigh"),
          las = 2,
          cex.lab = 0.75
        )
      } else {
        if (nrow(vecExp) == 3) {
          graphics::par(mfrow = c(2, 2))
        } 
        if (nrow(vecExp) <= 6 & nrow(vecExp) >= 4) {
          graphics::par(mfrow = c(2, 3))
        } 
        if (nrow(vecExp) > 6) {
          graphics::par(mfrow = c(3, 3))
        }
        for(s in 1:nrow(vecExp)){
          plot(
            vecExp[s,],
            t = 'b',
            lty = s,
            lwd = 1.5,
            pch = s,
            ylim = c(1, 4),
            ylab = "",
            xlab = "Month",
            yaxt = 'n',
            xaxt = 'n',
            main = rownames(vecExp)[s]
          )
          
          graphics::axis(1, at = 1:12, labels = month.abb, las = 2)
          graphics::axis(
            2,
            at = 1:4,
            labels = c('Low', "Moderate", "High", "Very\nHigh"),
            las = 2,
            cex.lab = 0.75
          )
          
        }
      }

      grDevices::dev.off()
    }

    if (ind[3]) {
      message(paste("Plotting Total Exposure with Important Variables..."))
      ### ONLY IMPORTANT VARS
      #load maps
      varMaps <- terra::rast(file.path(getwd(), x, 'Data',
                                       paste0('total_exposure_map_imp_var_', forecast_release, '_', forecast_init, '_', hindcast_release,'_',hindcast_yr_range, '.tif')))
      
      #load timeseries
      vecExp <- readRDS(
        file.path(getwd(), x, 'Data',
                  paste0('total_exposure_timeseries_imp_var_', forecast_release, '_', forecast_init, '_', hindcast_release,'_',hindcast_yr_range,'.rds')
        ))
      if(!is.null(nrow(vecExp))){
        rownames(vecExp)[1] <- 'All Stocks'
      }
      
      #plot map
      grDevices::pdf(
        paste0(
          file.path(getwd(), x, 'Figures',
                    paste0('total_exposure_imp_var_map_',forecast_release, '_', forecast_init, '_', hindcast_release,'_',hindcast_yr_range,'.pdf'))
        ),
        width = 8,
        height = 8
      )
      #graphics::par(mar = c(5, 5, 1.5, 0.5))
      #map
      terra::plot(
        varMaps,
        type = 'classes',
        levels = c("1", "2", "3", "4"),
        range = c(1, 4),
        col = cmocean::cmocean('matter')(4),
        ylim = c(35, 45),
        legend = F,
        xlab = expression('Longitude (' * degree * ')'),
        ylab = expression('Latitude (' * degree * ')'),
        pax = list(cex.axis = 1.5), cex.lab = 1.25,
        mar = c(3,3,1.5,0.5)
      )
      
      #add bathy contours, coastline, and stocks if necessary
      terra::contour(bathymetry, filled = F, levels = c(-1000, -100, -50), add = T)
      plot(coastline['id'], col = 'grey', add = T)
      if(!is.null(stocks)){
        terra::plot(stocks, add = T, lwd = 2)
      }
      
      #legend
      terra::plot(
        varMaps,
        type = 'classes',
        levels = c("1", "2", "3", "4"),
        range = c(1, 4),
        col = cmocean::cmocean('matter')(4),
        legend.only = TRUE, # <-- Draws only the legend elements
        plg = list(
          title = "Exposure",
          title.cex = 1.5,
          cex = 1.5,
          x = -68,
          y = 38, 
          legend = c("Low", "Moderate", "High", "Very High")
        )
      )
      grDevices::dev.off()
      
      #plot timeseries
      grDevices::pdf(
        paste0(
          file.path(getwd(), x, 'Figures',
                    paste0('total_exposure_imp_vars_timeseries_',forecast_release, '_', forecast_init, '_', hindcast_release,'_',hindcast_yr_range,'.pdf'))
        ),
        width = 11,
        height = 8
      )
      graphics::par(mar = c(5,5,2,2))
      
      if(!inherits(vecExp, 'matrix')){
        plot(
          vecExp,
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
          labels = c('Low', "Moderate", "High", "Very\nHigh"),
          las = 2,
          cex.lab = 0.75
        )
      } else {
        if (nrow(vecExp) == 3) {
          graphics::par(mfrow = c(2, 2))
        } 
        if (nrow(vecExp) <= 6 & nrow(vecExp) >= 4) {
          graphics::par(mfrow = c(2, 3))
        } 
        if (nrow(vecExp) > 6) {
          graphics::par(mfrow = c(3, 3))
        }
        for(s in 1:nrow(vecExp)){
          plot(
            vecExp[s,],
            t = 'b',
            lty = s,
            lwd = 1.5,
            pch = s,
            ylim = c(1, 4),
            ylab = "",
            xlab = "Month",
            yaxt = 'n',
            xaxt = 'n',
            main = rownames(vecExp)[s]
          )
          
          graphics::axis(1, at = 1:12, labels = month.abb, las = 2)
          graphics::axis(
            2,
            at = 1:4,
            labels = c('Low', "Moderate", "High", "Very\nHigh"),
            las = 2,
            cex.lab = 0.75
          )
          
        }
      }
      
      grDevices::dev.off()
    }

    if (ind[4]) {
      message(paste("Plotting Radar Plot of Relative Variable Importance..."))
      #load variable weights
      imp <- readRDS(
        file.path(getwd(), x, 'Data',
        'normalized_dynamic_variable_importance.rds')
      ) 
      cW <- imp[nrow(imp),]
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
      fmsb::radarchart(as.data.frame(cW), pfcol = scales::alpha('grey', 0.5), seg = 10)
      grDevices::dev.off()
    }
  } #end x
}