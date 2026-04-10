#' @title Normalize Dynamic Variable Weights
#' @description Normalizes the component model variable importance, as well as the model weights in the final ensemble, to create a weighted average vector of variable importance for the ensemble model. Options allow for static variables, such as time and space, or non-dynamic environmental variables such as bathymetry, to be excluded from the normalization and weighted average.
#'
#' @param vars character vector of variables to generate normalized and averaged ensemble importance for
#' @param ens_weights vector of component model weights in the ensemble
#' @param imp_flist character vector of file paths to variable importance for component models
#'
#' @return A vector representing the weighted average of normalized variable importance, representing variable importance in the final ensemble SDM.

normalize_variable_weights <- function(vars, ens_weights, imp_flist) {
  #set up data frame
  dfI <- as.data.frame(matrix(nrow = length(imp_flist), ncol = length(vars)))
  colnames(dfI) <- vars

  #pull importance vectors and add to data frame
  for (x in 1:length(imp_flist)) {
    load(imp_flist[x]) #imp
    if (inherits(imp, 'data.frame')) {
      imp.vec <- imp$rel.inf ### need to remove spatial variables
      names(imp.vec) <- imp$var
      for (y in 1:length(names(imp.vec))) {
        if (names(imp.vec)[y] %in% vars) {
          i <- which(vars == names(imp.vec)[y])
          dfI[x, i] <- imp.vec[y]
        }
      }
    } else {
      # imp <- range01(imp)
      for (y in 1:length(names(imp))) {
        if (names(imp)[y] %in% vars) {
          i <- which(vars == names(imp)[y])
          dfI[x, i] <- imp[y]
        }
      }
    }
  }

  #we don't have month and year or the static vars so remove those
  dfV <- replace(dfV, is.na(dfV), 0)
  dfV <- t(apply(dfV, 1, FUN = function(x) {
    x / sum(x, na.rm = T)
  }))
  dfV <- replace(dfV, is.na(dfV), 0)
  ws <- apply(dfV, MARGIN = 2, FUN = weighted.mean, w = ens_weights, na.rm = T) #weighted average of weights

  return(ws)
}
