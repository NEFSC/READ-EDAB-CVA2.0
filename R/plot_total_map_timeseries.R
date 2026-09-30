#' @title Make Vulnerability Plots
#' @description
#' Produces total Vulnerability maps and time series 
#'
#' @param map spatRaster to plot 
#' @param timeseries a vector or data.frame of timeseries data to plot if desired. Defaults to NULL, which just plots the map  
#' @param stocks spatVector of stock polygons to add to plot. Defaults to NULL meaning no stocks are available. 
#' @param stock_key named vector of stock long names and abbreviations. Used to label timeseries. 
#' @param fig_name name to save figure as. Figure will save to current working directory if desired directory is not included in the figure name
#' @param metric changes color pallete depending on which metric is being plotted. Must be 'exposure', 'sensitivity', or 'vulnerability' 
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
    width = ifelse(!is.null(timeseries), 11, 8),
    height = ifelse(!is.null(timeseries), 8, 8)
  )
  
  #if multi-panel plot is desired
  if(!is.null(timeseries)){
    #dynamic layout 
    if(is.null(stocks)){
      layout(matrix(c(1,1,2,2,1,1,3,3,1,1,3,3), byrow = T, nrow = 3, ncol = 4), height = c(2,2,2), width = c(3,3,2,2))
    } else {
      if(nrow(timeseries) <= 3){
        layout(matrix(c(1,1,2,2,1,1,3,3,1,1,4,4), byrow = T, nrow = 3, ncol = 4), height = c(2,2,2), width = c(3,3,2,2))
      }
      if(nrow(timeseries) > 3){
        layout(matrix(c(1,1,2,2,1,1,3,4,1,1,5,6), byrow = T, nrow = 3, ncol = 4), height = c(2,2,2), width = c(3,3,2,2))
      }
    } #end if stocks 
  } #end if timeseries
  
  pal <- switch(metric,
                'exposure' = 'dense',
                'sensitivity' = 'amp',
                'vulnerability' = 'matter')
  
  # 1. Hardcode the color mapping directly into the raster object
  terra::coltab(map) <- data.frame(
    value = 1:4,
    col = cmocean::cmocean(pal)(4)
  )
  
  # 2. Plot the map (remove col, breaks, range, and type arguments)
  terra::plot(
    map,
    ylim = c(35, 45),
    legend = FALSE, # coltab handles the map colors; we will build the legend separately
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
  legend(
    x = -69,
    y = 38,
    title = "Vulnerability",
    legend = c("Low", "Moderate", "High", "Very High"),
    fill = cmocean::cmocean(pal)(4),
    cex = 1.5,
    bty = "n" # Removes the box around the legend (optional)
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
          main = stock_key[names(stock_key) %in% rownames(timeseries)[x]]
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