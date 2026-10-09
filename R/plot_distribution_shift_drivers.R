#' @title Plot drivers of distribution shifts
#' @description
#' make gif of maps of species distribution and timeseries of range edges
#'
#' @param plot_name full file path to save gif to 
#' @param avg_core_polys,avg_range_polys polygons representing the core and range distributions to map; output from \code{calculate_distribution_shift}. 
#' @param poly_pals a two column data.frame of colors to plot for core and range polygons
#' @param avg_metrics data.frame output from \code{calculate_distribution_shift} containing x/y coordinates of the weighted centroids to add to the plot. should contain columns \code{centroid_x/y} and \code{timestamp}, which should match the desired legend names
#' @param hindcast_timeseries,forecast_timeseries data.frame output from \code{calculate_distribution_shift} at the desired timestamps for timeseries plots of core and total range area and latitude/longitude ranges. \code{forecast_timeseries} is optional. 
#' @param bathymetry,coastline spatRaster objects to add to maps
#'
#' @return a list containing a data.frame of metrics where the number of rows is equal to the number of timestamps, the core kde polygons, and the range kde polygons
#'
#'@export
#'
#'
plot_distribution_shift_drivers <- function(plot_name, sdm, core_polys, range_polys, metrics, bathymetry, coastline){
  gifski::save_gif(expr = 
             for(i in 1:nrow(metrics)){
               layout(matrix(c(1,1,1,1,2,4,3,5), byrow=F, nrow = 2, ncol = 4), heights = c(3,3), widths = c(1,1,1,1))
               par(mar=c(4,4,1,1), oma = c(0,0,3,0))
               #map of habitat shifts
               terra::plot(sdm[[i]], range = c(0,1),
                    pax = list(cex.axis = 2, 
                               yat = seq(35,45,by=1),
                               retro = T),
                    mar = c(3.1, 3.1, 2.1, 2.1), # Remove outer right margin space
                    xlab = '', 
                    ylab = '',
                    col = cmocean::cmocean('rain')(30),
                    plg = list(
                      x = c(-73, -65.5, 35.25, 36),
                      horiz = TRUE, 
                      
                      # Tick formatting - none added below because terra is difficult 
                      tick = 'none',
                      labels = '',
                      box.col = "black"
                    ))
               #add ticks/legend manually because terra is difficult
               tick_x <- seq(-73, -65.5, length.out = 5)
               segments(x0 = tick_x, y0 = 36, x1 = tick_x, y1 = 36 + 0.1, col = "black")
               
               #Add tick text on top of bar
               tick_vals <- seq(0,1, length.out = 5)
               text(x = tick_x, y = 36 + 0.25, labels = tick_vals, cex = 2)
               
               #Add title above tick text
               text(x = mean(c(-73, -65.5)), y = 36 + 0.5, 
                    labels = "Probability of Occurrence", 
                    cex = 2.5, font = 2)
               
               #abline(v = seq(-76,66,by=2), h = seq(36,44,by = 2), lwd = 0.5, lty = 3, col = 'grey')
               
               #add core and total area polygons + centroid
               terra::polys(core_polys[i], border = 'black', lwd = 3)
               terra::polys(range_polys[i], border = 'black', lwd = 3, lty = 3)
               
               points(centroid_y ~ centroid_x, data = metrics[i,], pch = 23, col = 'black', bg = 'green4', cex = 5)
               
               terra::polys(coastline, col = 'grey25')
               
               text(-76, 44, paste(month.abb[lubridate::month(metrics$timestamp[i])], lubridate::year(metrics$timestamp[i]), sep = ' '), cex = 4)
               
               add_legend(x = -71, y = 39, pch = c(NA, NA, 18), lty = c(1,3,NA), col = c('black', 'black', 'green4'), legend = c("Core Habitat (50% KDE)", 'Total Range (95% KDE)', "Weighted Centroid"), bty = 'n', cex = 2.5)
               
               #change in leading edges
               #latitude
               plot(leading_edge_y  ~ timestamp, data = metrics, t = 'n', pch = 19, ylab = 'Latitude', cex.axis = 1.5, cex.lab = 1.75, xlab = '', main = 'Northern Edge', cex.main = 1.75)
               lines(leading_edge_y  ~ timestamp, data = metrics[1:i,], t = 'b', pch = 19, col = 'grey')
               points(leading_edge_y  ~ timestamp, data = metrics[i,], pch = 19, cex = 2, col = 'red4')
               
               #longitude
               plot(leading_edge_x ~ timestamp, data = metrics, t = 'n', pch = 19, ylab = 'Longitude', cex.axis = 1.5, cex.lab = 1.75, xlab = '', main = 'Eastern Edge', cex.main = 1.75)
               lines(leading_edge_x  ~ timestamp, data = metrics[1:i,], t = 'b', pch = 19, col = 'grey')
               points(leading_edge_x  ~ timestamp, data = metrics[i,], pch = 19, cex = 2, col = 'red4')
               
               #change in trailing edges
               #latitude
               plot(trailing_edge_y  ~ timestamp, data = metrics, t = 'n', pch = 19, ylab = 'Latitude', cex.axis = 1.5, cex.lab = 1.75, xlab = 'Year', main = 'Southern Edge', cex.main = 1.75)
               lines(trailing_edge_y  ~ timestamp, data = metrics[1:i,], t = 'b', pch = 19, col = 'grey')
               points(trailing_edge_y  ~ timestamp, data = metrics[i,], pch = 19, cex = 2, col = 'red4')
               
               #longitude
               plot(trailing_edge_x ~ timestamp, data = metrics, t = 'n', pch = 19, ylab = 'Longitude', cex.axis = 1.5, cex.lab = 1.75, xlab = 'Year', main = 'Western Edge', cex.main = 1.75)
               lines(trailing_edge_x  ~ timestamp, data = metrics[1:i,], t = 'b', pch = 19, col = 'grey')
               points(trailing_edge_x  ~ timestamp, data = metrics[i,], pch = 19, cex = 2, col = 'red4')
               
             },
           width = 1440, height = 720, delay = 0.25, loop = T, progress = T,
           gif_file = plot_name
  )
}