#' @title Calculate leading and trailing edges of species distributions
#' @description
#' calculate trailing and leading edges of species distributions using cumulative density function
#'
#' @param r_single a spatRaster with one layer
#' @param level the probability to calculate KDE for
#'
#' @return a polygon representing the KDE at the given probability
#'
#'@export
#'

get_kde <- function(r_single, level) {
  # 1. Calculate sum over non-NA cells
  total_val <- terra::global(r_single, "sum", na.rm = TRUE)$sum
  if (is.na(total_val) || total_val == 0) {
    return(NA)
  }

  # 2. Extract values and isolate non-NA indices
  vals <- terra::values(r_single, mat = FALSE)
  valid_idx <- which(!is.na(vals))

  if (length(valid_idx) == 0) {
    return(0)
  }

  # 3. Work ONLY with valid cells to normalize and rank
  valid_vals <- vals[valid_idx] / total_val
  ord <- order(valid_vals, decreasing = TRUE)

  # 4. Calculate cumulative density
  cum_sum <- numeric(length(valid_vals))
  cum_sum[ord] <- cumsum(valid_vals[ord])

  # 5. Build cumulative raster initialized to NA
  r_cum <- terra::rast(r_single)
  terra::values(r_cum) <- NA
  r_cum[valid_idx] <- cum_sum

  # 6. Mask cells that fall within the threshold LEVEL
  r_mask <- r_cum <= level

  # CRITICAL: Mask out original NA regions (e.g. ocean) so they remain NA
  r_mask <- terra::mask(r_mask, r_single)

  # 7. Polygonize keeping ONLY the TRUE (1) foreground cells
  poly <- terra::as.polygons(r_mask, values = TRUE)

  # Keep only the true interior polygon region
  poly <- poly[!is.na(poly[[1]]) & poly[[1]] == 1, ]

  if (nrow(poly) == 0) {
    return(0)
  }

  # 8. Return poly
  return(poly)
}
