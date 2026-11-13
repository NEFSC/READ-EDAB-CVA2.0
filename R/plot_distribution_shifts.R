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
plot_distribution_shifts <- function(
  plot_name,
  avg_core_polys,
  avg_range_polys,
  poly_pals,
  avg_metrics,
  hindcast_timeseries,
  forecast_timeseries,
  bathymetry,
  coastline
) {
  if (!is.null(forecast_timeseries)) {
    tm_range <- c(
      min(hindcast_timeseries$timestamp, na.rm = T),
      max(forecast_timeseries$timestamp, na.rm = T)
    )
  } else {
    tm_range <- range(hindcast_timeseries$timestamp, na.rm = T)
  }

  # Failsafe: If the function crashes, forcefully close any open PDFs
  on.exit(
    while (grDevices::dev.cur() > 1) {
      grDevices::dev.off()
    },
    add = TRUE
  )

  plot_shift <- function() {
    # Set global outer margin to zero out excess space
    par(oma = c(0, 0, 0, 0))

    # ==========================================
    # MAP 1 (Left Half: x = 0 to 0.5, y = 0 to 1)
    # ==========================================
    par(fig = c(0, 0.5, 0, 1))
    # Update this for BOTH Map 1 and Map 2 sections:
    par(mar = c(3.5, 4.5, 1.5, 1), mgp = c(2.5, 0.75, 0))

    #map of habitat shifts - using KDE polygons
    #bathy contours as base layer
    terra::contour(
      bathymetry,
      filled = F,
      levels = c(-1000, -100, -50),
      cex.axis = 1.2,
      cex.lab = 1.25,
      xlab = expression('Longitude (' * degree * ')'),
      ylab = expression('Latitude (' * degree * ')'),
      xaxs = "i",
      yaxs = "i",
      asp = NA
    ) # Removes the 4% whitespace buffer
    # 1. Plot coastline first so habitat polygons overlay it
    # (add = T is implied/not needed for terra::polys)
    terra::polys(coastline, col = 'grey25')

    # 2. Add habitat polygons using polys() instead of plot()
    for (x in 1:length(avg_core_polys)) {
      terra::polys(
        avg_core_polys[x],
        border = scales::alpha(poly_pals$core[x], 0.5),
        lwd = 3,
        col = scales::alpha(poly_pals$core[x], 0.1)
      )
      #terra::polys(avg_range_polys[x], border = scales::alpha(poly_pals$range[x], 0.5), lwd = 3, lty = 3, col = NA)
    }

    points(
      centroid_y ~ centroid_x,
      data = avg_metrics,
      pch = 21,
      bg = scales::alpha(poly_pals$core, 0.75),
      col = 'black',
      cex = 2
    )

    #legend
    terra::add_legend(
      x = -77.5,
      y = 44.5,
      fill = poly_pals$core,
      legend = avg_metrics$timestep,
      cex = 1,
      bty = 'o',
      bg = 'white',
      box.col = 'black'
    )

    #time series
    # --- INSET 1 (Time Series inside Map 1) ---
    # Bottom right corner of Map 1: x = 0.25 to 0.48, y = 0.05 to 0.40
    par(fig = c(0.20, 0.48, 0.12, 0.45), new = TRUE)
    par(mar = c(2, 4.5, 1, 0.5), mgp = c(1.5, 0.5, 0)) # Widened left margin

    #change in areas through time
    #core
    if (!is.null(forecast_timeseries)) {
      y_range <- range(
        c(
          hindcast_timeseries$kde_core_area / 1000,
          forecast_timeseries$kde_core_area / 1000
        ),
        na.rm = T
      )
    } else {
      y_range <- range(hindcast_timeseries$kde_core_area / 1000, na.rm = T)
    }

    #annual hindcast
    plot(
      kde_core_area / 1000 ~ timestamp,
      data = hindcast_timeseries,
      t = 'b',
      pch = 19,
      cex = 1.25,
      lwd = 1.25,
      col = 'grey50',
      ylab = expression("Area (10"^3 * " km"^2 * ")"),
      cex.axis = 1,
      xlab = '',
      main = '',
      xlim = tm_range,
      ylim = y_range
    )

    # 2. Add this to color the first 5 points (e.g., dodgerblue)
    points(
      kde_core_area / 1000 ~ timestamp,
      data = hindcast_timeseries[1:5, ],
      pch = 19,
      cex = 1.5,
      col = poly_pals$range[1]
    )

    # 3. Add this to color the last 5 points (e.g., firebrick1)
    n <- nrow(hindcast_timeseries)
    points(
      kde_core_area / 1000 ~ timestamp,
      data = hindcast_timeseries[(n - 4):n, ],
      pch = 19,
      cex = 1.5,
      col = poly_pals$range[2]
    )

    #linear model
    m <- lm(kde_core_area ~ timestep, data = hindcast_timeseries)
    # Extract key statistics & plot if significant
    p_val <- summary(m)$coefficients["timestep", "Pr(>|t|)"]
    if (p_val <= 0.05 & !is.na(p_val)) {
      y_pred <- predict(m)
      lines(
        hindcast_timeseries$timestamp,
        y_pred / 1000,
        col = "black",
        lty = 2,
        lwd = 2
      )
    }

    if (!is.null(forecast_timeseries)) {
      #repeat for forecast
      #annual
      lines(
        kde_core_area / 1000 ~ timestamp,
        data = forecast_timeseries,
        t = 'b',
        pch = 17,
        cex = 1.25,
        lwd = 1.25,
        col = 'grey50'
      )
      n <- nrow(forecast_timeseries)
      points(
        kde_core_area / 1000 ~ timestamp,
        data = forecast_timeseries[(n - 4):n, ],
        pch = 17,
        cex = 1.5,
        col = poly_pals$range[3]
      )
      #linear model
      m <- lm(kde_core_area ~ timestep, data = forecast_timeseries)
      # Extract key statistics & plot if significant
      p_val <- summary(m)$coefficients["timestep", "Pr(>|t|)"]
      if (p_val <= 0.05 & !is.na(p_val)) {
        y_pred <- predict(m)
        lines(
          forecast_timeseries$timestamp,
          y_pred / 1000,
          col = "black",
          lty = 2,
          lwd = 2
        )
      }
    }

    ## ==========================================
    # MAP 2 (Right Half: x = 0.5 to 1, y = 0 to 1)
    # ==========================================
    par(fig = c(0.5, 1, 0, 1), new = TRUE)
    # Update this for BOTH Map 1 and Map 2 sections:
    par(mar = c(3.5, 4.5, 1.5, 1), mgp = c(2.5, 0.75, 0))
    #map of habitat shifts - using KDE polygons
    #bathy contours as base layer
    terra::contour(
      bathymetry,
      filled = F,
      levels = c(-1000, -100, -50),
      cex.axis = 1.2,
      cex.lab = 1.25,
      xlab = expression('Longitude (' * degree * ')'),
      ylab = expression('Latitude (' * degree * ')'),
      xaxs = "i",
      yaxs = "i",
      asp = NA
    ) # Removes the 4% whitespace buffer
    # 1. Plot coastline first so habitat polygons overlay it
    # (add = T is implied/not needed for terra::polys)
    terra::polys(coastline, col = 'grey25')

    # 2. Add habitat polygons using polys() instead of plot()
    for (x in 1:length(avg_core_polys)) {
      #terra::polys(avg_core_polys[x], border = scales::alpha(poly_pals$core[x], 0.75), lwd = 3, col = NA)
      terra::polys(
        avg_range_polys[x],
        border = scales::alpha(poly_pals$range[x], 0.5),
        lwd = 3,
        col = scales::alpha(poly_pals$core[x], 0.1)
      )
    }

    points(
      centroid_y ~ centroid_x,
      data = avg_metrics,
      pch = 21,
      bg = scales::alpha(poly_pals$core, 0.75),
      col = 'black',
      cex = 2
    )

    #legend
    terra::add_legend(
      x = -77.5,
      y = 44.5,
      fill = poly_pals$core,
      legend = avg_metrics$timestep,
      cex = 1,
      bty = 'o',
      bg = 'white',
      box.col = 'black'
    )

    # --- INSET 2 (Time Series inside Map 2) ---
    # Shifted left (0.68 to 0.95) to clear the right-side margin
    par(fig = c(0.70, 0.98, 0.12, 0.45), new = TRUE)
    par(mar = c(2, 4.5, 1, 0.5), mgp = c(1.5, 0.5, 0)) # Widened left margin
    #all (95% kde)
    if (!is.null(forecast_timeseries)) {
      y_range <- range(
        c(
          hindcast_timeseries$kde_all_area / 1000,
          forecast_timeseries$kde_all_area / 1000
        ),
        na.rm = T
      )
    } else {
      y_range <- range(hindcast_timeseries$kde_all_area / 1000, na.rm = T)
    }
    #annual hindcast
    plot(
      kde_all_area / 1000 ~ timestamp,
      data = hindcast_timeseries,
      t = 'b',
      pch = 19,
      cex = 1.25,
      lwd = 1.25,
      col = 'grey50',
      ylab = expression("Area (10"^3 * " km"^2 * ")"),
      cex.axis = 1,
      xlab = '',
      main = '',
      xlim = tm_range,
      ylim = y_range
    )

    # 2. Add this to color the first 5 points (e.g., dodgerblue)
    points(
      kde_all_area / 1000 ~ timestamp,
      data = hindcast_timeseries[1:5, ],
      pch = 19,
      cex = 1.5,
      col = poly_pals$range[1]
    )

    # 3. Add this to color the last 5 points (e.g., firebrick1)
    n <- nrow(hindcast_timeseries)
    points(
      kde_all_area / 1000 ~ timestamp,
      data = hindcast_timeseries[(n - 4):n, ],
      pch = 19,
      cex = 1.5,
      col = poly_pals$range[2]
    )

    #linear model
    m <- lm(kde_all_area ~ timestep, data = hindcast_timeseries)
    # Extract key statistics & plot if significant
    p_val <- summary(m)$coefficients["timestep", "Pr(>|t|)"]
    if (p_val <= 0.05 & !is.na(p_val)) {
      y_pred <- predict(m)
      lines(
        hindcast_timeseries$timestamp,
        y_pred / 1000,
        col = "black",
        lty = 2,
        lwd = 2
      )
    }

    if (!is.null(forecast_timeseries)) {
      #repeat for forecast
      #annual
      lines(
        kde_all_area / 1000 ~ timestamp,
        data = forecast_timeseries,
        t = 'b',
        pch = 17,
        cex = 1.25,
        lwd = 1.25,
        col = 'grey50'
      )
      n <- nrow(forecast_timeseries)
      points(
        kde_all_area / 1000 ~ timestamp,
        data = forecast_timeseries[(n - 4):n, ],
        pch = 17,
        cex = 1.5,
        col = poly_pals$range[3]
      )

      #linear model
      m <- lm(kde_all_area ~ timestep, data = forecast_timeseries)
      # Extract key statistics & plot if significant
      p_val <- summary(m)$coefficients["timestep", "Pr(>|t|)"]
      if (p_val <= 0.05 & !is.na(p_val)) {
        y_pred <- predict(m)
        lines(
          forecast_timeseries$timestamp,
          y_pred / 1000,
          col = "black",
          lty = 2,
          lwd = 2
        )
      }
    }
  } #end function

  #plot static
  #pdf
  grDevices::pdf(paste0(plot_name, '.pdf'), width = 10, height = 4.8)

  plot_shift()

  grDevices::dev.off()

  grDevices::png(
    paste0(plot_name, '.png'),
    width = 10,
    height = 4.8,
    units = 'in',
    res = 300
  )

  plot_shift()

  grDevices::dev.off()
}
