###make sdmtmb abundances 

###sdmspecific functions for container becuase this works better for some reason - only needs to be run once per time period because all models will use the same environmental data
makePredDF <- function(rasts, staticData, bathyR, bathy_max, mask){

      #make static vars (month/year) into rasters
  r <- subset(rasts[[1]][[1]], 1)
  rlon<-rlat<-r #copy r to rlon and rlat rasters [1]][1]which will contain the longitude and latitude
  xy<-xyFromCell(r,1:length(r)) #matrix of longitudes (x) and latitudes(y)
  rlon[]<-xy[,1] #raster of longitudes
  rlat[]<-xy[,2] #raster of latitides
  rMonth <- rYear <- r
  
    load(staticData) #staticVars object containing a list of rasters with the same extent as the environmental variables
    allDF <- NULL
        for(x in 1:nlayers(rasts[[1]][[1]])){ #all the rasters in rasts have the same number of layers so it doesn't matter which one we call
          
          #may need to change if names aren't always going to be month.year 
          mm.year <- strsplit(names(rasts[[1]][[1]])[x], split = '[.]') #all the rasterbricks in normVars also have the same names so again, doesn't matter which one we call
          
          mm <- as.numeric(gsub('X', '', mm.year[[1]][1])) 
          rMonth[] <- mm
          yr <- as.numeric(mm.year[[1]][2])
          rYear[] <- yr
          
          nStack <- vector(mode = 'list', length = length(rasts))
          for(n in 1:length(rasts)){
            nStack[[n]] <- subset(rasts[[n]][[1]], x)
          }
          nStack <- c(nStack, staticVars)
          nStack <- stack(nStack)
          names(nStack) <- c(names(rasts), names(staticVars))
          crs(nStack) <- crs(rasts[[1]][[1]])
          extent(nStack) <- extent(rasts[[1]][[1]])
          
          sr <- stack(rlon, rlat, rMonth, rYear, nStack)
          names(sr)[1:4] <- c("x", "y", "month", "year")
          
          if(mask){ #if mask == T
            #mask off waters deeper than 1000 m
            #i <- which(names(sr) == bathy_nm)
            sr <- replace(sr, abs(bathyR) > bathy_max, NA) #replace values with an absolute value greater than bathy_max with NA
          }
          
          ##convert rasterStack to dataframe to play well with model 
          srDF <- as.data.frame(rasterToPoints(sr)[,-c(1:2)])
          srDF <- srDF[complete.cases(srDF),]
            allDF <- rbind(allDF, srDF)
        }
          save(allDF, file = paste0('./Data/prediction_dataframe_', min(allDF$year, na.rm = T), '_', max(allDF$year, na.rm = T), '.RData'))
          print('allDF created and saved in Data directory')
            return(allDF)
    } #end function


predictSDM <- function(mod, df, staticData){
        load(staticData) #staticVars object containing a list of rasters with the same extent as the environmental variables
          Sys.time()
            require(sdmTMB)
          
          hsm <- vector(mode = 'list', length = length(unique(df$my)))
          for(x in 1:length(unique(df$my))){
            sub <- df[df$my == unique(df$my)[x],]
            coordinates(sub) <- ~x + y
            proj4string(sub) <- CRS("+proj=longlat +datum=WGS84 +no_defs ")
            hsm[[x]] <- raster::rasterize(x = sub, y = staticVars[[1]], field = sub$est)
            print(x)
          }
          
          return(hsm)
    }







