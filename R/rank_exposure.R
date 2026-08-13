#' @title Rank Raw Environmental Exposure
#' @description
#' Converts raw exposure from \code{calculate_raw_exposure} to a rank between 1-4 based on thresholds used in CVA1.0. This function also flips variables as needed such that positive exposure is always negative (i.e. both increasing temperature and decreasing oxygen concentrations result in positive exposure).
#'
#' @param exposure output from \code{calculate_raw_exposure}
#' @param flip TRUE/FALSE option to multiply exposure by -1 to ensure that positive exposure values represent negative variable change (increasing temperatures, decreasing oxygen concentrations, etc).
#'
#' @return A spatRaster with each layer containing ranked values between 1 - 4. 
#'
#'@export

rank_exposure <- function(
  exposure,
  flip = T
) {

    #change sign - negative = exposure to worse habitat?
    if (flip) {
      exposure <- exposure * -1
      }

    #rank
    QR <- terra::ifel(!is.na(exposure), 1, NA) #everything starts as 1 and we build from there 
    QR <- terra::ifel(exposure > 0.5 & exposure <= 1.5, QR + 1, QR + 0) #add one if exposure is between 0.5 and 1.5 to bring maximum to 2
    QR <- terra::ifel(exposure > 1.5 & exposure <= 2, QR + 1, QR + 0) #add one if exposure is between 1.5 and 2 to bring maximum to 3
    QR <- terra::ifel(exposure > 2, QR + 1, QR + 0) #add one if exposure is greater than 2 to bring maximum to 4

  names(QR) <- names(exposure)
  
  return(QR)
}
