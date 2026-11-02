#' @title Make Total Vulnerability Plots
#' @description
#' Produces total Vulnerability maps and time series 
#'
#' @param map spatRaster to plot. Can be NULL if plot_map is FALSE.
#' @param timeseries a vector or data.frame of timeseries data to plot if desired. Defaults to NULL.  
#' @param stocks spatVector of stock polygons to add to plot. Defaults to NULL. 
#' @param stock_key named vector of stock long names and abbreviations. Used to label timeseries. 
#' @param fig_name name to save figure as. 
#' @param metric changes color pallete depending on which metric is being plotted. Must be 'exposure', 'sensitivity', or 'vulnerability' 
#' @param coastline shapefile used to plot land in model prediction plots
#' @param bathymetry spatRaster file of bathymetry data; used to plot bathymetry in stock boundary plots
#'
#' @return Function does not return anything. Figures are saved to species-specific \code{figures} folder.
#'
#' @export
#'
plot_total_map_timeseries <- function(map = NULL, 
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
  
  if (is.null(map) && is.null(timeseries)) {
    stop("Both plot_map is FALSE and timeseries is NULL. Nothing to plot.")
  }
  
  plot_total <- function(){
    
    # --- 1. Layout Setup ---
    if (!is.null(map)) {
      # If multi-panel plot with map is desired
      if (!is.null(timeseries)) {
        if (is.null(stocks)) {
          layout(matrix(c(1,1,2,2,1,1,3,3,1,1,3,3), byrow = T, nrow = 3, ncol = 4), height = c(2,2,2), width = c(3,3,2,2))
        } else {
          if (nrow(timeseries) <= 3) {
            layout(matrix(c(1,1,2,2,1,1,3,3,1,1,4,4), byrow = T, nrow = 3, ncol = 4), height = c(2,2,2), width = c(3,3,2,2))
          }
          if (nrow(timeseries) > 3) {
            layout(matrix(c(1,1,2,2,1,1,3,4,1,1,5,6), byrow = T, nrow = 3, ncol = 4), height = c(2,2,2), width = c(3,3,2,2))
          }
        } 
      }
    } else {
      # Timeseries only layout
      if (is.null(stocks)) {
        par(mfrow = c(1, 1), mar = c(4, 4, 2, 2))
      } else {
        n_plots <- nrow(timeseries)
        if (n_plots <= 2) {
          par(mfrow = c(1, n_plots), mar = c(4, 4, 2, 2))
        } else if (n_plots <= 4) {
          par(mfrow = c(2, 2), mar = c(4, 4, 2, 2))
        } else {
          par(mfrow = c(2, ceiling(n_plots / 2)), mar = c(4, 4, 2, 2))
        }
      }
    }
    
    # --- 2. Map Plotting ---
    if (!is.null(map)) {
      pal <- switch(metric,
                    'exposure' = 'dense',
                    'sensitivity' = 'amp',
                    'vulnerability' = 'matter')
      
      terra::coltab(map) <- data.frame(
        value = 1:4,
        col = cmocean::cmocean(pal)(4)
      )
      
      terra::plot(
        map,
        ylim = c(35, 45),
        legend = FALSE, 
        pax = list(cex.axis = 1.5),
        cex.lab = 1.25,
        xlab = expression('Longitude (' * degree * ')'),
        ylab = expression('Latitude (' * degree * ')'),
        mar = c(3, 3, 1.5, 0.5)
      )
      
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
      
      legend(
        x = -69,
        y = 38,
        title = switch(metric,
                       'exposure' = 'Exposure',
                       'sensitivity' = 'Sensitivity',
                       'vulnerability' = 'Vulnerability'),
        legend = c("Low", "Moderate", "High", "Very High"),
        fill = cmocean::cmocean(pal)(4),
        cex = 1.5,
        bty = "n" 
      )
    }
    
    # --- 3. Timeseries Plotting ---
    if (!is.null(timeseries)) {
      if (is.null(stocks)) {
        plot(
          timeseries,
          t = 'b',
          lty = 8,
          lwd = 1.5,
          pch = 19,
          ylim = c(1, 4),
          ylab = ifelse(!plot_map, "Vulnerability", ""),
          xlab = "Month",
          yaxt = 'n',
          xaxt = 'n',
          main = 'Range'
        )
        
        graphics::axis(1, at = 1:12, labels = month.abb, las = 2)
        graphics::axis(
          2,
          at = 1:4,
          labels = c('L', "M", "H", "VH"),
          las = 2,
          cex.lab = 0.75
        )
      } else {
        for(x in 1:nrow(timeseries)){
          
          # Set title to 'Range' for the first panel, otherwise use the stock key
          panel_title <- if(x == 1) "Range" else stock_key[names(stock_key) %in% rownames(timeseries)[x]]
          
          plot(
            timeseries[x,],
            t = 'b',
            lty = 7+x,
            lwd = 1.5,
            pch = 18+x,
            ylim = c(1, 4),
            ylab = ifelse(!plot_map, "Vulnerability", ""),
            xlab = "Month",
            yaxt = 'n',
            xaxt = 'n',
            main = panel_title
          )
          
          graphics::axis(1, at = 1:12, labels = month.abb, las = 2)
          graphics::axis(
            2,
            at = 1:4,
            labels = c('L', "M", "H", "VH"),
            las = 2,
            cex.lab = 0.75
          )
        } 
      } 
    } 
    
  } 
  
  # --- 4. Figure Export ---
  
  # Adjust dimensions based on what is being plotted
  if (is.null(map)) {
    w <- 8
    h <- 6
  } else if (!is.null(timeseries)) {
    w <- 11
    h <- 8
  } else {
    w <- 8
    h <- 8
  }
  
  grDevices::pdf(
    paste0(fig_name, '.pdf'),
    width = w,
    height = h
  )
  plot_total()
  grDevices::dev.off()
  
  grDevices::png(
    paste0(fig_name, '.png'),
    width = w,
    height = h,
    units = 'in', 
    res = 300
  )
  plot_total()
  grDevices::dev.off()
}