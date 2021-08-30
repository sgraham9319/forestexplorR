#' Get UTM coordinates for individual trees
#' 
#' Calculates UTM coordinates for all trees in a user-provided mapping dataset
#' and returns a spatial object that can be plotted and joined with other
#' spatial datasets such as fine-scale topography data.
#' 
#' Location data for stands (argument \code{stand_locs}) must be provided in
#' either decimal degrees or UTM format. If your data are in another format 
#' such as degrees/minutes/seconds, please convert them to one of the accepted
#' formats. The \code{measurements::conv_unit} function is particularly
#' useful for achieving this.
#' 
#' @param tree_x_y Data frame containing mapping information for trees. Must
#' contain the columns \code{stand_id} (name of stand the tree is located in),
#' \code{x_coord} (local x coordinate), and \code{y_coord} (local y coordinate).
#' Any additional columns will appear unchanged in the output.
#' @param stand_locs Data frame containing locations of stand origins (i.e. 
#' local coordinates (0, 0)) and their y-azimuths (i.e. bearing from (0, 0) to
#' (0, maximum y)). Columns must be named \code{stand_id}, \code{y_azim}, 
#' \code{latitude}, \code{longitude} - see dataset \code{stand_locations} for
#' an example. Latitudes and longitudes must be in decimal degrees or UTM units.
#' @param original_crs Coordinate reference system of the location data in 
#' \code{stand_locs} as an ESPG code.
#' @param utm_crs ESPG code of the UTM zone in which the trees are located. If
#' the location data are in UTM, this will be the same as \code{original_crs}.
#' @return The input data frame \code{tree_x_y} with one additional column 
#' containing spatial information.
#' @examples
#' mapping_utm <- tree_utm(mapping, stand_locations, 4326, 32610)

tree_utm <- function(tree_x_y, stand_locs, original_crs, utm_crs){
  
  # Create spatial object of stand locations
  stand_locs_sf <- sf::st_as_sf(stand_locs, coords = c("longitude", "latitude"),
                                crs = original_crs)
  
  # Transform coordinate reference system to UTM if needed
  if(utm_crs == original_crs){
    stand_utm <- stand_locs_sf
  } else {
    stand_utm <- sf::st_transform(stand_locs_sf, crs = utm_crs)
  }
  
  # Split UTM coordinates into x and y columns 
  stand_utm$x <- as.vector(sf::st_coordinates(stand_utm)[,1])
  stand_utm$y <- as.vector(sf::st_coordinates(stand_utm)[,2])
  
  # Create vector of stands
  stands <- unique(mapping$stand_id)
  
  # Loop through stands
  for(i in 1:length(stands)){
    
    # Isolate focal stand
    focal_stand <- tree_x_y %>%
      filter(stand_id == stands[i])
    
    # Extract y azimuth for focal stand and convert to radians
    y_azim <- stand_locs[stand_locs$stand_id == stands[i], "y_azim"]
    y_azim <- y_azim * (pi / 180)
    
    # Calculate relative azimuths (in radians) from stand origin to each tree
    focal_stand$rel_az <- y_azim + atan(focal_stand$x_coord /
                                          focal_stand$y_coord)
    
    # Calculate distances from origin to each tree
    focal_stand$dists <- sqrt((focal_stand$x_coord ^ 2) +
                                (focal_stand$y_coord ^ 2))
    
    # Extract X and Y UTM coordinates of stand origin
    origin_x <- stand_utm[stand_utm$stand_id == stands[i], "x"][[1]]
    origin_y <- stand_utm[stand_utm$stand_id == stands[i], "y"][[1]]
    
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
  
  # Remove unneeded columns
  all_stands <- all_stands %>%
    select(-c(rel_az, dists))
  
  # Transform output to spatial object
  all_stands_utm <- sf::st_as_sf(all_stands, coords = c("x_utm", "y_utm"),
                                 crs = utm_crs)
}
