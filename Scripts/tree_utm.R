library(sf)

# Load location data for stands
stand_locations <- read.csv("Data/plot_locations.csv", stringsAsFactors = F)

test <- tree_utm(mapping, stand_locations, 4326, 32610)

tree_utm <- function(tree_x_y, stand_locs, original_crs, utm_crs){
  
  # Create spatial object of stand locations
  stand_locs_sf <- sf::st_as_sf(stand_locs, coords = c("longitude", "latitude"),
                                crs = original_crs)
  
  # Transform coordinate reference system to UTM
  stand_utm <- sf::st_transform(stand_locs_sf, crs = utm_crs)
  
  # Split UTM coordinates into x and y columns 
  stand_utm$x <- as.vector(st_coordinates(stand_utm)[,1])
  stand_utm$y <- as.vector(st_coordinates(stand_utm)[,2])
  
  # Create vector of stands
  stands <- unique(mapping$stand_id)
  
  # Loop through stands
  for(i in 1:length(stands)){
    
    # Isolate focal stand
    focal_stand <- tree_x_y %>%
      filter(stand_id == stands[i])
    
    # Extract y azimuth for focal stand and convert to radians
    y_azim <- stand_locs[stand_locs$stand_id == stands[i], "y_azim"][1]
    y_azim <- y_azim * (pi / 180)
    
    # Calculate relative azimuths (in radians) from stand origin to each tree
    focal_stand$rel_az <- y_azim + atan(focal_stand$x_coord /
                                          focal_stand$y_coord)
    
    # Calculate distances from origin to each tree
    focal_stand$dists <- sqrt((focal_stand$x_coord ^ 2) +
                                (focal_stand$y_coord ^ 2))
    
    # Extract X and Y UTM coordinates of stand origin
    origin_x <- stand_utm[stand_utm$stand_id == stands[i] & 
                            stand_utm$corner_id == "origin", "x"][[1]]
    origin_y <- stand_utm[stand_utm$stand_id == stands[i] & 
                            stand_utm$corner_id == "origin", "y"][[1]]
    
    # Calculate UTM X and Y
    focal_stand$x_utm <- origin_x + (sin(focal_stand$rel_az) *
                                       focal_stand$dists)
    focal_stand$y_utm <- origin_y + (cos(focal_stand$rel_az) *
                                       focal_stand$dists)
    
    # Add to combined output
    if(i == 1){
      all_stands <- focal_stand
    } else {
      all_stands <- bind_rows(all_stands, focal_stand)
    }
  }
  
  # Transform output to spatial object
  all_stands_utm <- sf::st_as_sf(all_stands, coords = c("x_utm", "y_utm"),
                                 crs = utm_crs)
}
