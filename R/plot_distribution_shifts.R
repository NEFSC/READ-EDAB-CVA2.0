#' @title Plot distribution shift 
#' @description
#' plot map of distribution shifts and timeseries of area and range metrics
#'
#' @param plot_name full file path to save plot to 
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
plot_distribution_shifts <- function(plot_name, avg_core_polys, avg_range_polys, poly_pals, avg_metrics, hindcast_timeseries, forecast_timeseries,bathymetry, coastline){
  
  if(!is.null(forecast_timeseries)){
    tm_range <- c(min(hindcast_timeseries$timestamp,na.rm = T), max(forecast_timeseries$timestamp, na.rm = T))
  } else {
    tm_range <- range(hindcast_timeseries$timestamp,na.rm = T)
  }
  
  #plot static
  pdf(plot_name,
  width = 12, 
  height = 6)
  
  ###set up layout
  layout(matrix(c(1,1,1,1,2,4,3,5), byrow=F, nrow = 2, ncol = 4), heights = c(3,3), widths = c(1,1,1,1))
  par(mar=c(4,4,1,1), oma = c(0,0,3,0))
  
  #map of habitat shifts - using KDE polygons 
  #bathy contours as base layer
  terra::contour(bathymetry, filled = F, levels = c(-1000, -100, -50), 
                 cex.axis = 2, cex.lab = 1.25,
                 xlab = '', 
                 ylab = '')
  # 1. Plot coastline first so habitat polygons overlay it 
  # (add = T is implied/not needed for terra::polys)
  terra::polys(coastline, col = 'grey25')
  
  # 2. Add habitat polygons using polys() instead of plot()
  for(x in 1:length(avg_core_polys)){
    terra::polys(avg_core_polys[x], border = scales::alpha(poly_pals$core[x], 0.75), lwd = 3, col = NA)
    terra::polys(avg_range_polys[x], border = scales::alpha(poly_pals$range[x], 0.5), lwd = 3, lty = 3, col = NA)
  }
  
  points(centroid_y ~ centroid_x, data = avg_metrics, pch = 21, bg = scales::alpha(poly_pals$core, 0.75), col = 'black', cex = 2)
  
  terra::add_legend(x = -70, y = 39.5,  
                    fill = poly_pals$core,
                    legend = avg_metrics$timestep,
                    cex = 1.5, bty = 'n')
  
  terra::add_legend(x = -70.2, y = 38, 
                    pch = c(NA, NA, 21),
                    lty = c(1, 3, NA),
                    lwd = c(2,2,NA),
                    bg = c(NA, NA, 'grey90'),
                    legend = c('Core Habitat', 'Species Range', 'Weighted Centroid'),
                    cex = 1.5, bty = 'n')
  
  #time series 
  #change in areas through time
  #core
  if(!is.null(forecast_timeseries)){
    y_range <- range(c(hindcast_timeseries$kde_core_area/1000, forecast_timeseries$kde_core_area/1000), na.rm = T)
  } else {
    y_range <- range(hindcast_timeseries$kde_core_area/1000,na.rm = T)
  }
  
  #annual hindcast 
  plot(kde_core_area/1000 ~ timestamp, data = hindcast_timeseries, t = 'b', pch = 19, cex = 1.5, lwd = 1.5, col = 'black', ylab = 'Area (10^3 km2)', cex.axis = 1.25, xlab = '', main = 'Core Habitat Area', xlim = tm_range, ylim = y_range)

  #linear model
  m <- lm(kde_core_area ~ timestep, data = hindcast_timeseries)
  # Extract key statistics & plot if significant
  p_val     <- summary(m)$coefficients["timestep", "Pr(>|t|)"]
  if(p_val <= 0.05){
    y_pred <- predict(m)
    lines(hindcast_timeseries$timestamp, y_pred/1000, col = "firebrick", lty = 2, lwd = 2)
  }
  
  if(!is.null(forecast_timeseries)){
    #repeat for forecast
    #annual 
    lines(kde_core_area/1000 ~ timestamp, data = forecast_timeseries, t = 'b', pch = 17, cex = 1.5, lwd = 1.5, col = 'goldenrod4')
    #linear model
    m <- lm(kde_core_area ~ timestep, data = forecast_timeseries)
    # Extract key statistics & plot if significant
    p_val     <- summary(m)$coefficients["timestep", "Pr(>|t|)"]
    if(p_val <= 0.05){
      y_pred <- predict(m)
      lines(forecast_timeseries$timestamp, y_pred/1000, col = "blue4", lty = 2, lwd = 2)
    }
  }
  
  #all (95% kde)  
  if(!is.null(forecast_timeseries)){
    y_range <- range(c(hindcast_timeseries$kde_all_area/1000, forecast_timeseries$kde_all_area/1000), na.rm = T)
  } else {
    y_range <- range(hindcast_timeseries$kde_all_area/1000,na.rm = T)
  }
  #annual hindcast 
  plot(kde_all_area/1000 ~ timestamp, data = hindcast_timeseries, t = 'b', pch = 19, cex = 1.5, lwd = 1.5, col = 'black', ylab = 'Area (10^3 km2)', cex.axis = 1.25, xlab = '', main = 'Habitat Range Area', xlim = tm_range, ylim = y_range)

  #linear model
  m <- lm(kde_all_area ~ timestep, data = hindcast_timeseries)
  # Extract key statistics & plot if significant
  p_val     <- summary(m)$coefficients["timestep", "Pr(>|t|)"]
  if(p_val <= 0.05){
    y_pred <- predict(m)
    lines(hindcast_timeseries$timestamp, y_pred/1000, col = "firebrick", lty = 2, lwd = 2)
  }
  
  if(!is.null(forecast_timeseries)){
    #repeat for forecast
    #annual 
    lines(kde_all_area/1000 ~ timestamp, data = forecast_timeseries, t = 'b', pch = 17, cex = 1.5, lwd = 1.5, col = 'goldenrod4', ylab = 'Area (10^3 km2)', cex.axis = 1.25, xlab = '')
  
    #linear model
    m <- lm(kde_all_area ~ timestep, data = forecast_timeseries)
    # Extract key statistics & plot if significant
    p_val     <- summary(m)$coefficients["timestep", "Pr(>|t|)"]
    if(p_val <= 0.05){
      y_pred <- predict(m)
      lines(forecast_timeseries$timestamp, y_pred/1000, col = "blue4", lty = 2, lwd = 2)
    }
  }
  
  #change in ranges
  #latitude
  if(!is.null(forecast_timeseries)){
    y_range <- range(c(hindcast_timeseries$range_y_deg, forecast_timeseries$range_y_deg), na.rm = T)
  } else {
    y_range <- range(hindcast_timeseries$range_y_deg,na.rm = T)
  }
  #hindcast
  plot(range_y_deg ~ timestamp, data = hindcast_timeseries, t = 'b', pch = 19, cex = 1.5, lwd = 1.5, col = 'black', ylab = 'Degrees', cex.axis = 1.25, xlab = 'Year', main = 'Latitude Range', xlim = tm_range, ylim = y_range)

  #linear model
  m <- lm(range_y_deg ~ timestep, data = hindcast_timeseries)
  # Extract key statistics & plot if significant
  p_val     <- summary(m)$coefficients["timestep", "Pr(>|t|)"]
  if(p_val <= 0.05){
    y_pred <- predict(m)
    lines(hindcast_timeseries$timestamp, y_pred, col = "firebrick", lty = 2, lwd = 2)
  }
  if(!is.null(forecast_timeseries)){
    #repeat for forecast
    #annual 
    lines(range_y_deg~ timestamp, data = forecast_timeseries, t = 'b', pch = 17, cex = 1.5, lwd = 1.5, col = 'goldenrod4', ylab = 'Degrees', cex.axis = 1.25, xlab = '')
    #linear model
    m <- lm(range_y_deg ~ timestep, data = forecast_timeseries)
    # Extract key statistics & plot if significant
    p_val     <- summary(m)$coefficients["timestep", "Pr(>|t|)"]
    if(p_val <= 0.05){
      y_pred <- predict(m)
      lines(forecast_timeseries$timestamp, y_pred, col = "blue4", lty = 2, lwd = 2)
    }
  }
  
  #longitude
  if(!is.null(forecast_timeseries)){
    y_range <- range(c(hindcast_timeseries$range_x_deg, forecast_timeseries$range_x_deg), na.rm = T)
  } else {
    y_range <- range(hindcast_timeseries$range_x_deg,na.rm = T)
  }
  #hindcast
  plot(range_x_deg ~ timestamp, data = hindcast_timeseries, t = 'b', pch = 19, cex = 1.5, lwd = 1.5, col = 'black', ylab = '', cex.axis = 1.25, xlab = 'Year', main = 'Longitude Range', xlim = tm_range, ylim = y_range)
  #linear model
  m <- lm(range_x_deg ~ timestep, data = hindcast_timeseries)
  # Extract key statistics & plot if significant
  p_val     <- summary(m)$coefficients["timestep", "Pr(>|t|)"]
  if(p_val <= 0.05){
    y_pred <- predict(m)
    lines(hindcast_timeseries$timestamp, y_pred, col = "firebrick", lty = 2, lwd = 2)
  }
  if(!is.null(forecast_timeseries)){
    #repeat for forecast
    #annual 
    lines(range_x_deg~ timestamp, data = forecast_timeseries, t = 'b', pch = 17, cex = 1.5, lwd = 1.5, col = 'goldenrod4', ylab = 'Degrees', cex.axis = 1.25, xlab = '')
    #linear model
    m <- lm(range_x_deg ~ timestep, data = forecast_timeseries)
    # Extract key statistics & plot if significant
    p_val     <- summary(m)$coefficients["timestep", "Pr(>|t|)"]
    if(p_val <= 0.05){
      y_pred <- predict(m)
      lines(forecast_timeseries$timestamp, y_pred, col = "blue4", lty = 2, lwd = 2)
    }
  }
  
  dev.off()
}