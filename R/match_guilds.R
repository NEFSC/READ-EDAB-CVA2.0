#' @title Subset Species/Environmental Dataset to Ecologically-Relevant Dynamic Variables
#' @description
#' This function reduces the number of environmental covariates in the species/environmental data frame based on defined feeding and habitat guilds
#'
#' @param spp_env data.frame of fisheries and environmental data
#' @param spp Species name to add to log files and save data to correct directory (see vignette for recommended directory set up)
#' @param spp_col name of column with species name in spp_guild key
#' @param spp_guild file path to csv file containing a list of species and corresponding feeding and habitat guilds. See vignette for recommended set up of this file
#' @param feeding_key,habitat_key file path to csv for feeding and habitat guild keys respectively. See vignette for recommended set up of these files.
#' @param feeding_col,habitat_col columns names indicating the guild column in the feeding and habitat guild csvs
#' @param static_vars a vector with the column names to be retained
#' @param pa_col column name for presence/absence column
#'
#' @return a data frame that contains the static variables listed and only environmental covariates associated with the species' feeding and habitat guilds


match_guilds <- function(spp_env, spp, spp_col = 'name', spp_guild, feeding_key, feeding_col = 'Feeding.Guild', habitat_key,  habitat_col = 'Habitat.Guild', static_vars = c('x', 'y', 'month', 'year', 'bathy', 'rugosity', 'dist2coast'), pa_col = 'value'){

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
