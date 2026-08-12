#' @title Make Map or Timeseries of Total Species-Specific Exposure
#'
#' @description Uses the logic rule from CVA1.0 to combine variable-specific exposures to create a total exposure map or timeseries
#'
#' @param type designates desired output, must equal 'map' or 'timeseries'
#' @param variable_exposure If \code{type == 'map'}, a spatRaster output from \code{make_variable_exposure(type == 'map')}. If \code{type == 'timeseries'}, the matrix or list output from \code{make_variable_exposure(type == 'timeseries')}.
#' @param count_all TRUE/FALSE to use \code{weights} and \code{wThreshold} to subset variables to only important variables
#' @param variable_weights output from \code{combine_weights} - a vector of variable weights in ensemble SDM
#' @param weight_threshold numeric value used to subset weights, variables with weights less to or equal to this value will be excluded from total exposure calculation
#'
#' @return If \code{type == 'map'}, the output is a raster representing total exposure across space. If \code{type == 'timeseries'}, the output is a vector representing total exposure across time if a single matrix is supplied in variable_exposure, or a matrix with a number of rows equal to the length of the list supplied in variable_exposure.
#'
#'@export

make_total_exposure <- function(
  type,
  variable_exposure,
  count_all,
  variable_weights,
  weight_threshold
) {
  if (type == 'map') {
    if (count_all) {
      mapSub <- variable_exposure
    } else {
      wi <- which(weights >= weight_threshold)
      mapSub <- variable_exposure[wi]
    }

    #count how many layers have each rank within each cell
    hr <- sum(mapSub >= 3.5)
    hh <- sum(mapSub >= 3)
    md <- sum(mapSub >= 2.5)

    #apply logic rule from Hare et al 2015
    expL <- ifel(!is.na(hr), 1, NA) #everything starts as 1 (low)
    expL <- ifel(md >= 2, expL + 1, expL + 0) #add 1 if moderate threshold met, max now = 2
    expL <- ifel(hh >= 2, expL + 1, expL + 0) #add 1 if high threshold met, max now = 3
    expL <- ifel(hr >= 3, expL + 1, expL + 0) #add 1 if very high threshold met, max now = 4
    

    return(expL)
  }

  if (type == 'timeseries') {
    if(inherits(variable_exposure, 'list')){ #if variable_exposure is a list, meaning that the variable-level exposures were averaged within different stock polygons, then repeat the procedure for each matrix in the list
      vMat <- matrix(nrow = length(variable_exposure), ncol = 12)
      for(v in 1:length(variable_exposure)){
        if (count_all) {
          #if counting all included factors and not taking weight into account
          matSub <- variable_exposure[[v]]
        } else {
          wi <- which(weights > weight_threshold)
          matSub <- variable_exposure[[v]][wi, ]
        }
        
        #count each rank in each column
        hr <- hh <- md <- vector(length = ncol(matSub))
        for (x in 1:ncol(matSub)) {
          hr[x] <- length(which(matSub[, x] >= 3.5))
          hh[x] <- length(which(matSub[, x] >= 3))
          md[x] <- length(which(matSub[, x] >= 2.5))
        }
        
        #apply logic rule
        expV <- rep(1, times = ncol(matSub))
        expV <- replace(expV, md >= 2, 2)
        expV <- replace(expV, hh >= 2, 3)
        expV <- replace(expV, hr >= 3, 4)
        
        vMat[v,] <- expV
        return(vMat)
      }
      rownames(vMat) <- names(variable_exposure)
      colnames(vMat) <- month.abb
    } else {
      if (count_all) {
        #if counting all included factors and not taking weight into account
        matSub <- variable_exposure
      } else {
        wi <- which(weights > weight_threshold)
        matSub <- variable_exposure[wi, ]
      }
  
      #count each rank in each column
      hr <- hh <- md <- vector(length = ncol(matSub))
      for (x in 1:ncol(matSub)) {
        hr[x] <- length(which(matSub[, x] >= 3.5))
        hh[x] <- length(which(matSub[, x] >= 3))
        md[x] <- length(which(matSub[, x] >= 2.5))
      }
  
      #apply logic rule
      expV <- rep(1, times = ncol(matSub))
      expV <- replace(expV, md >= 2, 2)
      expV <- replace(expV, hh >= 2, 3)
      expV <- replace(expV, hr >= 3, 4)
  
      names(expV) <- month.abb
      return(expV)
    } #end if 
  } #end if timeseries
}
