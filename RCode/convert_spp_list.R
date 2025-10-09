
create_var_list <- function(varLong, varShort, names){
  require(tibble)
  #convert species columns into list to iterate over within targets workflow 
  
  vList <- NULL
  for(x in 1:length(varLong)){
    tb <- tibble::tribble(~Long.Name, ~Short.Name, varLong[x], varShort[x])
    vList <- rbind(vList, tb)
  }
  #remove all spaces 
 # sList$Name <- gsub(' ', '', sList$Name)
  
  names(vList) <- names
  return(vList)
}

