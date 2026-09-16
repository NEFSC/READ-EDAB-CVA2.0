#' @title Make Stock Polygons
#'
#' @description Make polygons for each stock for a species based on provided data.frame. Saves  
#'
#' @param key data.frame listing all of the desired species, stocks, and associated polygons in an ID column 
#' @param species_col,stock_col character string matching columns in \code{key} containing species and stock names 
#' @param id_col character string matching the ID column in \code{key} use to match to polygon names 
#' @param polygons a spatVector containing all the potential polygons to match to 
#' @param poly_id a character string matching the ID column in \code{polygons} to match polygons 
#' @param plot TRUE/FALSE turning on/off plotting code
#'
#' @return saves species-specific spatVector to working directory. returns log of species and stocks found in key
#'
#'@export
#'

make_stock_polygons <- function(key, species_col, stock_col, id_col, polygons, poly_id, plot = T, bathymetry, coastline){
  names <- unique(key[,species_col])
  
  log <- NULL
  for(x in names){
    spp <- key[key[,species_col] == x,] #subset key to desired species 
    
    stocks <- unique(spp[,stock_col]) #find unique stocks 
    stock.polys <- NULL 
    for(s in stocks){
      spp.stock <- spp[spp[,stock_col] == s,] #subset to stocks
      
      spp.stock.shp <- polygons[polygons[[poly_id]][,1] %in% spp.stock[,id_col],] #subset geometry 
      
      # --- NEW SAFETY CATCH ---
      if (nrow(spp.stock.shp) == 0) {
        warning(paste("No matching polygons found for stock:", s, "- skipping."))
        next 
      }
      # ------------------------
      
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
    
    terra::writeVector(stock.polys, filename = paste0(gsub(' ', '', x), '.shp'), overwrite = TRUE)
    
    log <- rbind(log, c(x, paste(stocks, ' ', collapse = ' ')))
    
    if(plot){
      grDevices::pdf(
        paste0(
          file.path(getwd(), 'maps',
                    paste0('stock_map_',gsub(' ', '', x),'.pdf'))
        ),
        width = 8,
        height = 8
      )
      
      #add bathy contours, coastline, and stocks if necessary
      terra::plot(stock.polys, lwd = 2,
                  col = cmocean::cmocean('rain')(length(stock.polys)),
                  pax = list(cex.axis = 1.5), cex.lab = 1.25,
                  xlab = expression('Longitude (' * degree * ')'),
                  ylab = expression('Latitude (' * degree * ')'),
                  mar = c(3,3,1.5,0.5),
                  ylim = c(35, 45),
                  xlim = c(-78, -65),
                  legend = F)
      terra::contour(bathymetry, filled = F, levels = c(-1000, -100, -50), add = T)
      terra::plot(coastline['id'], col = 'grey', add = T)
      
      #legend
      legend(x = -68, 
             y = 38, 
             legend = stock.polys$stock_area, 
             fill = cmocean::cmocean('rain')(length(stock.polys)), 
             title = "Stocks",
             cex = 1.5,
             bty = 'o',
             bg = 'white')
      
      
      grDevices::dev.off()
    }
    
  }#end x
  
  
  return(log)
  
} #end function