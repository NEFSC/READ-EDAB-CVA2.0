#' @title Helper: Repair Column Names
#' @description
#' Ensure column names are formatted correctly
#'
#' @param nms a vector of character strings to convert 
#' @param targets a vector of desired names
#'
#' @return a character vector of corrected names
#'
#'@export
#'

repair_names <- function(nms, targets) {
  res <- as.character(nms)
  for (i in seq_along(res)) {
    if (res[i] %in% targets) {
      next
    }
    clean_target <- make.names(targets)
    if (res[i] %in% clean_target) {
      res[i] <- targets[clean_target == res[i]][1]
    } else if (res[i] %in% gsub(" ", ".", targets)) {
      res[i] <- targets[gsub(" ", ".", targets) == res[i]][1]
    }
  }
  return(res)
}