#' @title Make Vulnerability Plots
#' @description
#' Produces total Vulnerability maps and time series 
#'
#' @param map spatRaster to plot 
#' @param timeseries a vector or data.frame of timeseries data to plot if desired. Defaults to NULL, which just plots the map  
#' @param stocks spatVector of stock polygons to add to plot. Defaults to NULL meaning no stocks are available. 
#' @param stock_key named vector of stock long names and abbreviations. Used to label timeseries. 
#' @param fig_name name to save figure as. Figure will save to current working directory if desired directory is not included in the figure name
#' @param metric changes figure range depending on which metric is being plotted. Must be 'exposure', 'sensitivity', or 'vulnerability' 
#' @param coastline shapefile used to plot land in model prediction plots
#' @param bathymetry spatRaster file of bathymetry data; used to plot bathymetry in stock boundary plots

#' @return Function does not return anything. Figures are saved to species-specific \code{figures} folder.
#'
#'@export
#'
plot_total_map_timeseries <- function(map, 
                               timeseries = NULL, 
                               stocks = NULL, 
                               stock_key, 
                               fig_name,
                               metric,
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
    width = ifelse(!is.null(timeseries), 11, 6),
    height = ifelse(!is.null(timeseries), 6, 11)
  )
  
  #if multi-panel plot is desired
  if(!is.null(timeseries)){
    #dynamic layout 
    if(is.null(stocks)){
      layout(matrix(c(1,1,2,2,1,1,3,3,1,1,3,3), byrow = T, nrow = 3, ncol = 4), height = c(2,2,2), width = c(2,2,2,2))
    } else {
      if(nrow(timeseries) <= 3){
        layout(matrix(c(1,1,2,2,1,1,3,3,1,1,4,4), byrow = T, nrow = 3, ncol = 4), height = c(2,2,2), width = c(2,2,2,2))
      }
      if(nrow(timeseries) > 3){
        layout(matrix(c(1,1,2,2,1,1,3,4,1,1,5,6), byrow = T, nrow = 3, ncol = 4), height = c(2,2,2), width = c(2,2,2,2))
      }
    } #end if stocks 
  } #end if timeseries
  
  #map
  terra::plot(
    map,
    type = 'classes',
    levels = c("1", "2", "3", "4"),
    range = c(1, 4),
    col = cmocean::cmocean('matter')(4),
    ylim = c(35, 45),
    legend = F,
    pax = list(cex.axis = 1.5),
    cex.lab = 1.25,
    xlab = expression('Longitude (' * degree * ')'),
    ylab = expression('Latitude (' * degree * ')'),
    mar = c(3, 3, 1.5, 0.5)
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
  
  #legend
  terra::plot(
    map,
    type = 'classes',
    levels = c("1", "2", "3", "4"),
    range = c(1, 4),
    col = cmocean::cmocean('matter')(4),
    legend.only = TRUE, # <-- Draws only the legend elements
    plg = list(
      title = "Vulnerability",
      title.cex = 1.5,
      cex = 1.5,
      x = -68,
      y = 38,
      legend = c("Low", "Moderate", "High", "Very High")
    )
  )
  
  #timeseries
  if(!is.null(timeseries)){
    if(is.null(stocks)){
      #plot vector
      plot(
        timeseries,
        t = 'b',
        lty = 8,
        lwd = 1.5,
        pch = 19,
        ylim = c(1, 4),
        ylab = "",
        xlab = "Month",
        yaxt = 'n',
        xaxt = 'n',
        main = 'Range'
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
      for(x in 1:nrow(timeseries)){
        #plot vector
        plot(
          timeseries[x,],
          t = 'b',
          lty = 7+x,
          lwd = 1.5,
          pch = 18+x,
          ylim = c(1, 4),
          ylab = "",
          xlab = "Month",
          yaxt = 'n',
          xaxt = 'n',
          main = stock_key[names(stock_key) %in% rownames(totV)[x]]
        )
        
        graphics::axis(1, at = 1:12, labels = month.abb, las = 2)
        graphics::axis(
          2,
          at = 1:4,
          labels = c('Low', "Moderate", "High", "Very\nHigh"),
          las = 2,
          cex.lab = 0.75
        )
      } #end x
    } #end if stocks 
  } #end if timeseries 
  
  grDevices::dev.off()
}