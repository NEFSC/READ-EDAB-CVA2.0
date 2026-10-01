#' @title Make Variable Vulnerability Plots
#' @description
#' Produces variable Vulnerability maps and time series 
#'
#' @param data either a multi-layer spatRaster or a matrix containing the map data or timeseries data to plot
#' @param timeseries a vector or data.frame of timeseries data to plot if desired. Defaults to NULL, which just plots the map  
#' @param stocks spatVector of stock polygons to add to plot. Defaults to NULL meaning no stocks are available. 
#' @param stock_key named vector of stock long names and abbreviations. Used to label timeseries. 
#' @param fig_name name to save figure as. Figure will save to current working directory if desired directory is not included in the figure name
#' @param coastline shapefile used to plot land in model prediction plots
#' @param bathymetry spatRaster file of bathymetry data; used to plot bathymetry in stock boundary plots

#' @return Function does not return anything. Figures are saved to species-specific \code{figures} folder.
#'
#'@export
#'
plot_variable_map_timeseries <- function(data, 
                                      stocks = NULL, 
                                      stock_key, 
                                      variable_df,
                                      fig_name,
                                      coastline, 
                                      bathymetry){
  
  # Failsafe: If the function crashes, forcefully close any open PDFs
  on.exit(
    while (grDevices::dev.cur() > 1) {
      grDevices::dev.off()
    },
    add = TRUE
  )
  
  
  grDevices::pdf(
    fig_name,
    width = 8,
    height = 11
  )     
  
  if(inherits(data, 'SpatRaster')){
    message('Plotting Variable Exposure Maps')
    #set up panels according to the number of variables
    if (terra::nlyr(data) <= 6) {
      graphics::par(
        mfrow = c(2, 3)
      )
    } else {
      graphics::par(
        mfrow = c(3, 3)
      )
    }
    
    for (y in 1:terra::nlyr(data)) {
      #get full name of variable
      i <- variable_df$Short.Name %in% names(data)[y]
      
      draw_legend <- y == terra::nlyr(data)
      
      #map
      #graphics::par(plt = c(0.1, 0.98, 0.1, 0.95))
      terra::plot(
        data[[y]],
        type = 'continuous',
        range = c(1, 4),
        col = cmocean::cmocean('deep')(64),
        mar = c(2.5, 2.5, 1.5, 0.5), # Explicitly set margins inside terra::plot
        legend = FALSE,
        xlab = expression('Longitude (' * degree * ')'),
        ylab = expression('Latitude (' * degree * ')'),
        main = variable_df$Long.Name[i],
        pax = list(cex.axis = 1.5, xat = seq(-80, -60, by = 2)),
        cex.lab = 1.25
      )
      
      #add bathy contours, coastline, and stocks if necessary
      terra::contour(
        bathymetry,
        filled = F,
        levels = c(-1000, -100, -50),
        add = T
      )
      plot(coastline['id'], col = 'grey', add = T)
      if (!is.null(stocks)) {
        terra::plot(stocks, add = T, lwd = 2)
      }
    } #end y
    # 2. Draw the legend independently if it is the last panel
    if (draw_legend) {
      graphics::par(mgp = c(3, 0.1, 0))
      terra::plot(
        data[[y]],
        type = 'continuous',
        range = c(1, 4),
        col = cmocean::cmocean('deep')(64),
        mar = c(2.5, 2.5, 1.5, 0.5), # Explicitly set margins inside terra::plot
        legend.only = TRUE, # <-- Draws only the legend elements
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
    }
    
  } else { #end if spatRaster
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
      } 
      
    } else {
        #if data is not a spatraster, then its a matrix or list of matrices
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
      
      graphics::axis(1, at = 1:12, labels = month.abb, las = 2, cex.lab = 0.5)
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

  }#end if not spatraster

  grDevices::dev.off()
}