match_guilds <- function(spp_env, spp, spp_col = 'name', spp_guild, feeding_key, feeding_col = 'Feeding.Guild', habitat_key,  habitat_col = 'Habitat.Guild', static_vars = c('x', 'y', 'month', 'year', 'bathy', 'rugosity', 'dist2coast'), pa_col = 'value'){
  #spp_env - a data frame with target species presence/absence matched to all possible environmental covariates
  #spp the name of the target species
  #spp_col - name of the column with the species name to match to the guild csv
  #spp_guild - the name of a csv file containing the species name and associated feeding and habitat guilds
  #feeding/habitat_key - csv files that list the environmental covariates associated with each feeding and habitat guild (two seperate files)
  #feeding/habitat_col - character string (or index?) of the column name associated with feeding/habitat guilds in spp_guild
  #static_vars - a character vector of the variables to retain regardless of feeding guild - should include positional data
  
  guilds <- read.csv(spp_guild) #load in species list 
  
  g <- which(guilds[,spp_col] %in% spp) #find row associated with target species
  
  #isolate feeding guild
  fGuild <- guilds[g, feeding_col]
  
  #isolate habitat guild
  hGuild <- guilds[g, habitat_col]
  
  ##isolate covariates for each guild type based on keys 
  #Feeding guilds: 
  feeding <- read.csv(feeding_key)
  i <- which(colnames(feeding) == fGuild)
  fInd <- feeding[nzchar(feeding[,i]), i] #get covariates 

  
  #Habitat guilds: 
  habitat <- read.csv(habitat_key)
  i <- which(colnames(habitat) == hGuild)
  hInd <- habitat[nzchar(habitat[,i]), i] #get covariates 
  
  ### subset spp_env
  ind <- colnames(spp_env) %in% c(pa_col, static_vars, fInd, hInd)
  se <- spp_env[,ind]
  
  return(se)
  
}

remove_corr <- function(se, pa_col, xy_col, month_col, year_col){
    #but first, remove correlated variables 
    ind <- which(colnames(se) == pa_col | colnames(se) == xy_col[1] | colnames(se) == xy_col[2] | colnames(se) == month_col | colnames(se) == year_col)
    corInd <- findCorrelation(cor(se[,-ind]), names = T) #find correlated variables 
    if(length(corInd) != 0){
      se <- se[,-which(colnames(se) == corInd)]
    }
  return(se)
}

clean_data <- function(se, pa_col){
  paDF <- se[se[,pa_col] != 0,]
  #switch it back to 0/1 now that we only have grid cells that were sampled 
  paDF[,pa_col] <- replace(paDF[,pa_col], paDF[,pa_col] == 1, 0)
  paDF[,pa_col] <- replace(paDF[,pa_col], paDF[,pa_col] == 2, 1)
  
  paDF2 <- paDF[complete.cases(paDF),]
  
  return(paDF)
}