#' @title Make SDm Plots
#' @description
#' Produces SDM result and residual plots 
#'
#' @param sdm,obs spatRaster(s) to plot 
#' @param stocks spatVector of stock polygons to add to plot. Defaults to NULL meaning no stocks are available. 
#' @param fig_name,hist_name name(s) to save maps and histograms (if plotting) as without extension. 
#' @param type character string describing what to plot. Options are 'model' or 'residuals'; if \code{type = 'residuals'}, \code{obs} is necessary
#' @param coastline shapefile used to plot land in model prediction plots
#' @param bathymetry spatRaster file of bathymetry data to add to plots 

#' @return Function does not return anything. PDFs and PNGs of figures are saved to path outlined in fig_name/hist_name. If no path is provided, figures will save in working directory.
#'
#'@export
#'
plot_sdms <- function(sdm, obs = NULL, 
                                      xy_col = NULL,
                                      time_col = NULL,
                      panel_names = month.abb,
                                      stocks = NULL,  
                                      fig_name,
                                      hist_name = NULL,
                                      type,
                                      coastline, 
                                      bathymetry){
  
  # Failsafe: If the function crashes, forcefully close any open PDFs
  on.exit(
    while (grDevices::dev.cur() > 1) {
      grDevices::dev.off()
    },
    add = TRUE
  )
  
  if(type == 'model'){
    rast <- sdm
    pal = 'rain'
    rng = c(0,1)
    tks = c(0, 0.25, 0.5, 0.75, 1)
  } else { #if type == 'residuals'
    if(inherits(obs, 'SpatRaster')){ #if obs is already a spatraster, just substract
      rast <- obs - sdm
    } else { #if obs is a data.frame, convert to spatraster using sdm as template
      #append sdm results to data frame
      preds <- build_preds_df(obs, xy_col = xy_col, sdm)
      #calculate residuals 
      preds$residuals <- preds$pa - preds$predicted
      
      #convert to spatRaster and average
      template_r <- sdm[[1]]
      #avg residuals
      avgR <- vector(mode = 'list', length = length(unique(preds[,time_col]))) #going by the number of unique values in the desired time_col
      for (y in 1:length(avgR)) {
        sub <- preds[preds[,time_col] == y, ]
        pts <- terra::vect(
          sub,
          geom = xy_col,
          crs = 'EPSG:4326' #assumption here
        )
        avgR[[y]] <- terra::rasterize(
          pts,
          template_r,
          field = 'residuals',
          fun = mean
        )
      } #end for
      avgR <- terra::rast(avgR)
      names(avgR) <- panel_names
      
      rast <- avgR
    
      #set color pallete & range
      pal = 'balance'
      rng = c(-1,1)
      tks = c(-1, -0.5, 0, 0.5, 1)
      
      ## 1. Bonus residual histogram inner plotting function
      draw_hist <- function() {
        graphics::hist(
          preds$residuals,
          main = '',
          xlab = 'Residuals',
          xlim = c(-1, 1)
        )
        graphics::abline(
          v = mean(preds$residuals, na.rm = T),
          lty = 2,
          col = 'red4'
        )
        graphics::legend(
          'topleft',
          legend = paste0(
            'Mean (+/- SD) =\n',
            round(mean(preds$residuals, na.rm = T), 2),
            ' +/- ',
            round(stats::sd(preds$residuals, na.rm = T), 2)
          ),
          bty = 'n'
        )
      }
      
      # Save Histogram as PDF
      grDevices::pdf(paste0(hist_name, ".pdf"), width = 6, height = 6)
      draw_hist()
      grDevices::dev.off()
      
      # Save Histogram as PNG
      grDevices::png(paste0(hist_name, ".png"), width = 6, height = 6, units = "in", res = 300)
      draw_hist()
      grDevices::dev.off()
    }
  }
  
  ## 2. Main maps inner plotting function
  draw_maps <- function() {
    graphics::par(
      mfrow = c(4, 3),
      mar = c(2.2, 2.2, 1, 0.5),
      oma = c(0, 0, 0, 0)
    )
    for (y in 1:terra::nlyr(rast)) {
      
      draw_legend <- y == terra::nlyr(rast)
      
      terra::plot(
        rast[[y]],
        type = 'continuous',
        range = rng,
        col = cmocean::cmocean(pal)(64),
        mar = c(2.5, 2.5, 1.5, 0.5), # Explicitly set margins inside terra::plot
        legend = FALSE,
        xlab = expression('Longitude (' * degree * ')'),
        ylab = expression('Latitude (' * degree * ')'),
        main = panel_names[y],
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
      terra::plot(coastline['id'], col = 'grey', add = T)
      if (!is.null(stocks)) {
        terra::plot(stocks, add = T, lwd = 2)
      }
    } #end y
    
    # Draw the legend independently if it is the last panel
    if (draw_legend) {
      graphics::par(mgp = c(3, 0.1, 0))
      terra::plot(
        rast[[terra::nlyr(rast)]],
        type = 'continuous',
        range = rng,
        col = cmocean::cmocean(pal)(64),
        mar = c(2.5, 2.5, 1.5, 0.5),
        legend.only = TRUE, 
        plg = list(
          title = ifelse(type == 'model', 'Mean Probability\nof Occurance', 'Mean Residuals'),
          title.cex = 1.25,
          cex = 1.25,
          horizontal = TRUE,
          x = -73.5,
          y = 36.5,
          at = tks,
          n = 4,
          size = c(1, 2.5)
        )
      )
    }
  }
  
  # Save Maps as PDF
  grDevices::pdf(paste0(fig_name, ".pdf"), width = 8, height = 11)
  draw_maps()
  grDevices::dev.off()
  
  # Save Maps as PNG
  grDevices::png(paste0(fig_name, ".png"), width = 8, height = 11, units = "in", res = 300)
  draw_maps()
  grDevices::dev.off()
}