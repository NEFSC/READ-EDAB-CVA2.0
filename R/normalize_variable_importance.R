#' @title Normalize Dynamic Variable Weights
#' @description Normalizes the component model variable importance, as well as the model weights in the final ensemble, to create a weighted average vector of variable importance for the ensemble model. Options allow for static variables, such as time and space, or non-dynamic environmental variables such as bathymetry, to be excluded from the normalization and weighted average.
#'
#' @param vars character vector of variables to generate normalized variable importance for
#' @param ens_weights vector of component model weights in the ensemble
#' @param imp_list a named list of importance outputs from different models. names should match model names function will handle differences between them. 
#'
#' @return A vector representing the weighted average of normalized variable importance, representing variable importance in the final ensemble SDM.
#'
#'@export

normalize_variable_importance <- function(vars, ens_weights, imp_list) {
  #set up data frame
  v <- data.frame(var = vars)
  
  #go through list of outputs and merge with v - since all.x = T for all, this will subset the call to only the desired variables 
  for (x in 1:length(imp_list)) {
    imp <- imp_list[[x]]
    if (inherits(imp, 'data.frame') & grepl('BRT', names(imp_list)[x])) {
      v <- merge(v, imp, by = 'var', all.x = T)
    }
    if (inherits(imp, 'data.frame') & grepl('SDMTMB', names(imp_list)[x])) {
      v <- merge(v, imp[,1:2], by.x = 'var', by.y = 'Variable', all.x = T)
    } 
    if(!inherits(imp, 'data.frame')){
      imp.df <- data.frame(var = names(imp)[!is.na(names(imp))], var.imp = imp[!is.na(names(imp))])
      v <- merge(v, imp.df, by = 'var', all.x = T)
    }
    #print(x)
  }
  
  colnames(v)[-1] <- names(imp_list)
  
  #normalize within models - some are already like this so they won't change 
  v <- replace(v, is.na(v), 0)
  
  dfI <- t(apply(v[,-1], 2, FUN = function(x) {
    x / sum(x)
  }))
  colnames(dfI) <- v$var
  
  #get weighted average of all according to component model weights 
  ws <- apply(dfI, MARGIN = 2, FUN = weighted.mean, w = ens_weights, na.rm = T) #weighted average of weights
  
  #combine 
  dfI <- rbind(dfI, ws)
  rownames(dfI)[nrow(dfI)] <- "ENSEMBLE"
  
  return(dfI)
}
