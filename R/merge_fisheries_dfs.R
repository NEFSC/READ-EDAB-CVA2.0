#' @title Merge Fisheries Presence/Absence/Effort Data Frames
#' @description Combines data.frames from \code{build_fisheries_df} across multiple sources
#'
#' @param df_dir full path to directory containing the data.frames to merge
#'
#' @return returns a single data.frame combining presence/absence data.frames from all sources. Data are aggregated within grid cells for each month and year such that each grid cell reflects if a species was present or absent within that grid cell across data sources
#'
#'
#'@export

merge_fisheries_dfs <- function(df_dir) {
  
  flist <- dir(df_dir, full.names = T)
  data <- NULL
  for(x in 1:length(flist)){
    d <- read.csv(flist[x])
    d$source <- flist[x]
    data <- rbind(data, d)
  }
  
  paMax <- stats::aggregate(data, by = list(data$year, data$month, data$gridID), FUN = max) #combine by grid cell ID to get max presence/absence within each grid cell
  paMax <- paMax[,-grep('Group', names(paMax))]

  return(paMax)
}
