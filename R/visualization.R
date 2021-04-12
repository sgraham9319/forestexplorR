#=========================================================
# Convert stand corner locations to decimal degrees format
#=========================================================

stand_dms_to_dd <- function(data){
  
  data %>%
    
    # Convert to long format
    pivot_longer(cols = 3:10, names_to = "location_id",
                 values_to = "deg_min_sec") %>%
    
    # Remove cardinal directions from location IDs
    mutate(corner_id = ifelse(
      stringr::str_detect(location_id, "_lat") == TRUE,
      str_sub(location_id, end = -5),
      str_sub(location_id, end = -6)),
      
      # Add column for lat or long
      coord = ifelse(
        stringr::str_detect(deg_min_sec, "W") == TRUE,
        "lon_dd",
        "lat_dd"),
      
      # Make deg_min_sec values negative for westings and southings
      deg_min_sec = case_when(
        str_detect(deg_min_sec, "W") == TRUE ~
          paste("-", deg_min_sec, sep = ""),
        str_detect(deg_min_sec, "S") == TRUE ~
          paste("-", deg_min_sec, sep = ""),
        TRUE ~ deg_min_sec),
      
      # Remove unwanted symbols from deg_min_sec
      deg_min_sec = str_replace(deg_min_sec, '\"', ""),
      deg_min_sec = str_replace(deg_min_sec, "\'", " "),
      deg_min_sec = str_replace(deg_min_sec, "\\.", " "),
      deg_min_sec = str_sub(deg_min_sec, end = -2),
      
      # Calculate decimal degrees
      dec_deg = as.numeric(conv_unit(deg_min_sec, "deg_min_sec", "dec_deg"))
    ) %>%
    
    # Select required columns
    select(stand_id, y_azim, corner_id, coord, dec_deg) %>%
    
    # Separate lats and longs
    pivot_wider(names_from = coord, values_from = dec_deg)
}

#===========================
# Create polygons for stands
#===========================

polygons <- function(data){
  
  data %>%
    
    # Group by stand to create multipoint objects
    group_by(stand_id) %>%
    
    # Need to include summarize for group_by to take effect
    summarise(y_azim = y_azim[1]) %>%
    
    # Use convex hull method to form polygon around each set of four points
    st_convex_hull() %>%
    
    # Return only stand IDs and geometry
    select(stand_id)
}

#==========================================================
# Plot one stand with trees colored by a specified variable
#==========================================================

utm_mapping <- function(tree_x_y, stand, color_var){
  
  # Isolate focal stand
  focal_stand <- tree_x_y[tree_x_y$stand_id == stand, ]
  
  # Calculate degrees to radians conversion factor
  cf <- pi / 180
  
  # Extract y azimuth for focal stand and convert to radians
  y_azim <- stand_locs_raw[stand_locs_raw$stand_id == stand, "y_azim"]
  y_azim <- y_azim * cf
  
  # Calculate relative azimuths (in radians) from stand origin to each tree
  focal_stand$rel_az <- y_azim + atan(focal_stand$x_coord / focal_stand$y_coord)
  
  # Calculate distances from origin to each tree
  focal_stand$dists <- sqrt((focal_stand$x_coord ^ 2) + (focal_stand$y_coord ^ 2))
  
  # Extract X and Y UTM coordinates of stand origin
  origin_x <- stand_utm[stand_utm$stand_id == stand & 
                          stand_utm$corner_id == "origin", "x"][[1]]
  origin_y <- stand_utm[stand_utm$stand_id == stand & 
                          stand_utm$corner_id == "origin", "y"][[1]]
  
  # Calculate UTM X and Y
  focal_stand$x_utm <- origin_x + (sin(focal_stand$rel_az) * focal_stand$dists)
  focal_stand$y_utm <- origin_y + (cos(focal_stand$rel_az) * focal_stand$dists)
  
  # Transform output to spatial object
  focal_stand_utm <- st_as_sf(focal_stand, coords = c("x_utm", "y_utm"), crs = 32610)
  
  # Plot result
  mapview(focal_stand_utm, zcol = color_var) + 
    mapview(stand_polygons$geometry[stand_polygons$stand_id == stand])
  
}
