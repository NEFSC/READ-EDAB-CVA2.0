#' @title Calculate leading and trailing edges of species distributions
#' @description
#' calculate trailing and leading edges of species distributions using cumulative density function
#'
#' @param coords a vector of x or y coordinates
#' @param weights a vector of sdm values predicted at x/y coordinates.
#' @param probs the desired probabilities to calculate as leading and trailing edges. defaults to 5 and 95%.
#'
#' @return a vector with length equal to the number of probabilities supplied with the x or y values corresponding to the provided probabilities
#'
#'@export

get_edges <- function(coords, weights, probs = c(0.05, 0.95)) {
  # 1. Sort coordinates in ascending order
  ord <- order(coords)
  coords_sorted <- coords[ord]
  weights_sorted <- weights[ord]

  # 2. Compute the cumulative density distribution (normalized from 0 to 1 by the total sum of the weights)
  cum_weights <- cumsum(weights_sorted) / sum(weights_sorted)

  # 3. Find the coordinate values corresponding to the target probabilities
  sapply(probs, function(p) {
    coords_sorted[which(cum_weights >= p)[1]]
  })
}
