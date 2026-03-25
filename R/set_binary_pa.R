#' @title Convert Species/Environmental Data Frame to Presence/Absence
#' @description
#' This function converts the presence/absence/effort (0-2) data to presence/absence data (0-1)
#'
#' @param se data.frame of fisheries and environmental data
#' @param pa_col column name for presence/absence column
#'
#' @return a data frame where presence/absence has been set to 0 for absent and 1 for present

set_binary_pa <- function(se, pa_col){
  paDF <- se[se[,pa_col] != 0,] #remove unsampled cells (pa_col == 0)
  #switch it back to 0/1 now that we only have grid cells that were sampled
  paDF[,pa_col] <- replace(paDF[,pa_col], paDF[,pa_col] == 1, 0)
  paDF[,pa_col] <- replace(paDF[,pa_col], paDF[,pa_col] == 2, 1)

  return(paDF)
}

