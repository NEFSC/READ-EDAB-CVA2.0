#' @title Make Stock Polygons
#'
#' @description Make polygons for each stock for a species based on provided data.frame. Saves  
#'
#' @param key data.frame listing all of the desired species, stocks, and associated polygons in an ID column 
#' @param species_col,stock_col character string matching columns in \code{key} containing species and stock names 
#' @param id_col character string matching the ID column in \code{key} use to match to polygon names 
#' @param polygons a spatVector containing all the potential polygons to match to 
#' @param poly_id a character string matching the ID column in \code{polygons} to match polygons 
#'
#' @return saves species-specific spatVector to working directory. returns log of species and stocks found in key
#'
#'@export
#'

make_stock_polygons <- function(key, species_col, stock_col, id_col, polygons, poly_id){
  names <- unique(key[,species_col])
  
  log <- NULL
  for(x in names){
    spp <- key[key[,species_col] == x,] #subset key to desired species 
    
    stocks <- unique(spp[,stock_col]) #find unique stocks 
    stock.polys <- NULL 
    for(s in stocks){
      spp.stock <- spp[spp[,stock_col] == s,] #subset to stocks
      
      spp.stock.shp <- polygons[polygons[[poly_id]][,1] %in% spp.stock[,id_col],] #subset geometry 
      
      # Project to EPSG:5070 before processing
      spp.stock.proj <- terra::project(spp.stock.shp, "EPSG:5070")
      
      # Buffer by 100 meters instead of degreesk
      spp.stock.agg <- spp.stock.proj |> 
        terra::buffer(width = 100) |> 
        terra::aggregate() |> 
        terra::fillHoles() |> 
        terra::buffer(width = -100)
      
      # Project back to original CRS if necessary
      spp.stock.agg <- terra::project(spp.stock.agg, terra::crs(polygons))
      
      # Assign the stock name as a new attribute column
      spp.stock.agg$stock_area <- s
      
      # Combine into a single multi-row SpatVector
      if(is.null(stock.polys)){
        stock.polys <- spp.stock.agg
      } else {
        stock.polys <- c(stock.polys, spp.stock.agg)
      }
    } #end s
    names(stock.polys) <- stocks
    stock.polys <- terra::vect(stock.polys)
    
    terra::writeVector(stock.polys, filename = paste0(gsub(' ', '', x), '.shp'), overwrite = T)
    
    log <- rbind(log, c(x, paste(stocks, ' ', collapse = ' ')))
    
  }#end x
  
  
  return(log)
  
} #end function