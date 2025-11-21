create_spp_list <- function(species.csv){
  require(tibble)
  #convert species columns into list to iterate over within targets workflow 
  spp <- read.csv(species.csv)
  
  sList <- NULL
    for(x in 1:nrow(spp)){
      tb <- tibble::tribble(~Name, ~Name.Options, spp$Common.Name[x], quote(list(spp$Common.Name[x], spp$SCI_NAME[x], spp$Alternate.Name[x])))
      sList <- rbind(sList, tb)
    }
  #remove all spaces 
  sList$Name <- gsub(' ', '', sList$Name)

  names(sList) <- c('Name', 'Name.Options')
  return(sList)
}

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


combine_vars <- function(vars, names, type = 'norm'){
  varList <- vector(mode = 'list', length = length(vars))
  
  #i <- grepl('norm', names(vars)) #find desired variable type 
  for(x in 1:length(vars)){
    varList[[x]] <- vars[[x]] #add desired variable type to list  
  }
  
  names(varList) <- names
  return(varList)
}

isolate_output <- function(target, name){
  nms <- vector(length = length(target))
  for(x in 1:length(target)){
    nms[x] <- names(target[[x]])
  }
  i <- grep(name, nms)
  out <- target[[i]]
  return(out)
}

save_output <- function(target, path){
  save(target, file = path)
  return(path)
}
