match_guilds <- function(spp_env, spp, spp_col = 'SCI_NAME', spp_guild, feeding_key, feeding_col = 'Feeding.Guild', habitat_key,  habitat_col = 'Habitat.Guild', static_vars = c('x', 'y', 'month_num', 'year', 'bathy', 'rugosity', 'coast_dist'), pa_col = 'value'){
  #spp_env - a data frame with target species presence/absence matched to all possible environmental covariates
  #spp the name of the target species
  #spp_col - name of the column with the species name to match to the guild csv
  #spp_guild - the name of a csv file containing the species name and associated feeding and habitat guilds
  #feeding/habitat_key - csv files that list the environmental covariates associated with each feeding and habitat guild (two seperate files)
  #feeding/habitat_col - character string (or index?) of the column name associated with feeding/habitat guilds in spp_guild
  #static_vars - a character vector of the variables to retain regardless of feeding guild - should include positional data
  
  guilds <- read.csv(spp_guild) #load in species list 
  
  g <- which(guilds[,spp_col] == spp) #find row associated with target species
  
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
  se <- spp_env[,c(pa_col, static_vars, fInd, hInd)]
  
  return(se)
  
}