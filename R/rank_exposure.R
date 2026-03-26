#' @title Rank Raw Environmental Exposure
#' @description
#' Converts raw exposure from \code{calculate_raw_exposure} to a rank between 1-4 based on thresholds used in CVA1.0. This function also flips variables as needed such that positive exposure is always negative (i.e. both increasing temperature and decreasing oxygen concentrations result in positive exposure).
#'
#' @param exposure output from \code{calculate_raw_exposure}
#' @param flip TRUE/FALSE option to multiply exposure by -1 to ensure that positive exposure values represent negative variable change (increasing temperatures, decreasing oxygen concentrations, etc).
#' @param no_flip_vars character vector naming which vectors should NOT be flipped - counterintuitive, yes, but this list was actually shorter than the variables that needed to be flipped for NECVA2.0. Names should match the names in \code{exposure}
#'
#' @return A list of rasterStacks with each layer containing ranked values between 1 - 4. The length of the list is equal to the length of the lists supplied as \code{exposure}

rank_exposure <- function(exposure, flip = T, no_flip_vars = c('bottomT', 'surfaceT', 'bottomArg', 'MLD')){
  rankE <- vector(mode = 'list', length = length(exposure)) #list of vectors for
  for(x in 1:length(exposure)){
    s <- exposure[[x]]

    #change sign - negative = exposure to worse habitat?
    if(flip){
      if(names(exposure)[x] %in% no_flip_vars){
        s <- s
      } else {
        s <- s*-1
      }
    }

    mapE <- vector(mode = 'list', length = 12) #list for maps
    for(m in 1:12){
      r <- raster::subset(s, m)
      #rank
      QR <- replace(r, r <= 0.5, 1)
      QR <- replace(QR, r > 0.5 & r <= 1.5, 2)
      QR <- replace(QR, r > 1.5 & r <= 2, 3)
      QR <- replace(QR, r > 2, 4)

      mapE[[m]] <- QR
    }

    rankE[[x]] <- raster::stack(mapE)

  }
  names(rankE) <- names(exposure)
  return(rankE)
}
