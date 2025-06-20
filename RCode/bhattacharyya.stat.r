### test Bhattacharyya's coefficient ranges from 0 to 1, with 0 being no overlap in distributions and 1 being perfect overlap. 
#source("http://tguillerme.github.io/R/bhatt.coef.R")
# data must have a presence absence column called "presabs" and column names that match the names in vars
# vars is character vector of the desired variables

bhattacharyya.stat <- function(data, response, vars){

  bh <- list()
  for (i in 1:length(vars)){
    idx.pres <- which(data[,response] == 1)
    data.pres <- data[idx.pres,] %>% dplyr::select(vars[i])
    data.pres <- data.pres[!is.na(data.pres)]
    
    idx.abs <- which(data[,response] == 0)
    data.abs <- data[idx.abs,] %>% dplyr::select(vars[i])
    data.abs <- data.abs[!is.na(data.abs)]
    
    bh[[i]] <- bhatt.coeff(data.pres, data.abs)
    names(bh[[i]]) <- vars[i]
  }

  return(bh)
}


