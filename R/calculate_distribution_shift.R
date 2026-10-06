#' @title Calculate distribution shift metrics
#' @description
#' calculate trailing and leading edges of species distributions, desired kdes, kde areas, and weighted centroids
#'
#' @param r_stack a spatraster of sdm results
#' @param trailing_p,leading_p probabilities for leading and trailing edges; defaults to 5% and 95%
#' @param core_p,range_p probabilities associated with desired core and range KDEs; defaults to 50% and 95% respectively
#' @param poly_dir location to save kde polygons
#'
#' @return a list containing a data.frame of metrics where the number of rows is equal to the number of timestamps, the core kde polygons, and the range kde polygons
#'
#'@export
#'

analyze_distribution_dynamics <- function(
  r_stack,
  trailing_p = 0.05,
  leading_p = 0.95,
  core_p = 0.5,
  all_p = 0.95,
  poly_dir = NULL
) {
  # 1. Extract coordinates and cell density values
  xy <- terra::xyFromCell(r_stack, 1:ncell(r_stack))
  vals_mat <- terra::values(r_stack)
  vals_mat[is.na(vals_mat)] <- 0

  n_layers <- terra::nlyr(r_stack)
  layer_names <- names(r_stack)
  if (is.null(layer_names)) {
    layer_names <- paste0("Time_", 1:n_layers)
  }

  results_list <- core_poly_list <- all_poly_list <- vector("list", n_layers)

  # 2. Iterate through each timestep layer
  for (i in 1:n_layers) {
    r_single <- r_stack[[i]]
    w <- vals_mat[, i]
    valid <- w > 0

    if (!any(valid)) {
      next
    }

    x_val <- xy[valid, 1]
    y_val <- xy[valid, 2]
    w_val <- w[valid]

    # Weighted Centroid (Center of Gravity)
    tot_w <- sum(w_val)
    cx <- sum(x_val * w_val) / tot_w
    cy <- sum(y_val * w_val) / tot_w

    # Weighted Quantiles (Trailing and Leading Edges)
    y_quantiles <- get_edges(y_val, w_val, probs = c(trailing_p, leading_p))
    x_quantiles <- get_edges(x_val, w_val, probs = c(trailing_p, leading_p))

    # Calculate exact KDEs
    kde_core <- get_kde(r_single, level = core_p)
    kde_all <- get_kde(r_single, level = all_p)

    #Calculate KDE areas
    kde_core_area <- sum(expanse(kde_core, unit = 'km'))
    kde_all_area <- sum(expanse(kde_all, unit = 'km'))

    results_list[[i]] <- data.frame(
      layer = layer_names[i],
      timestep = i,
      centroid_x = cx,
      centroid_y = cy,
      trailing_edge_y = y_quantiles[1],
      leading_edge_y = y_quantiles[2],
      trailing_edge_x = x_quantiles[1],
      leading_edge_x = x_quantiles[2],
      kde_core_area = kde_core_area,
      kde_all_area = kde_all_area
    )

    core_poly_list[[i]] <- kde_core
    all_poly_list[[i]] <- kde_all
  }

  df <- do.call(rbind, results_list)

  # 3. Calculate step-by-step centroid displacement (t to t+1)
  df$centroid_dx <- c(diff(df$centroid_x), NA)
  df$centroid_dy <- c(diff(df$centroid_y), NA)

  # Convert coordinate centroids to a SpatVector point object
  centroid_pts <- vect(
    as.matrix(df[, c("centroid_x", "centroid_y")]),
    type = "points",
    crs = crs(r_stack)
  )

  # Calculate true geodesic distance (in meters) between consecutive points
  step_dists_m <- sapply(1:(nrow(df) - 1), function(i) {
    distance(centroid_pts[i], centroid_pts[i + 1])
  })

  # Add to data frame in kilometers
  df$centroid_step_dist_km <- c(step_dists_m / 1000, NA)

  # Compass bearing (0° N, 90° E, 180° S, 270° W)
  bearing_rad <- c(atan2(diff(df$centroid_x), diff(df$centroid_y)), NA)
  df$centroid_bearing_deg <- (bearing_rad * 180 / pi) %% 360

  # 4. Cumulative & Net Displacement relative to baseline (Time 1)
  total_dists_m <- sapply(1:(nrow(df) - 1), function(i) {
    distance(centroid_pts[1], centroid_pts[i])
  })
  df$total_displacement_km <- c(total_dists_m / 1000, NA)

  step_dists <- df$centroid_step_dist_km[!is.na(df$centroid_step_dist_km)]
  df$cum_distance_traveled <- c(0, cumsum(step_dists_m / 1000))

  # 5. Boundary & Area shifts
  df$leading_y_shift <- c(diff(df$leading_edge_y), NA)
  df$trailing_y_shift <- c(diff(df$trailing_edge_y), NA)
  df$area_core_change <- c(diff(df$kde_core_area), NA)
  df$area_all_change <- c(diff(df$kde_all_area), NA)

  # Combine polygons into a single SpatVector
  core_polygons <- do.call(rbind, core_poly_list)
  all_polygons <- do.call(rbind, all_poly_list)

  # Save polygons to disk if path is provided
  if (!is.null(poly_dir) && !is.null(all_polygons)) {
    terra::writeVector(
      all_polygons,
      filename = paste0(poly_dir, '/kde_all_', all_p * 100, '.gpkg'),
      overwrite = TRUE
    )
    terra::writeVector(
      core_polygons,
      filename = paste0(poly_dir, '/kde_core_', core_p * 100, '.gpkg'),
      overwrite = TRUE
    )
    cat(sprintf(
      "Successfully saved 50%% and 95%% KDE polygons to: %s\n",
      poly_dir
    ))
  }

  # Return BOTH the summary dataframe and the SpatVector polygon object
  return(list(
    metrics = df,
    core_polygons = core_polygons,
    all_polygons = all_polygons
  ))
}
