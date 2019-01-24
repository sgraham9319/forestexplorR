
#==================================================
# Reorient plots with origin as western-most corner
#==================================================

# The reorient function takes a data frame of the within-plot X, Y locations
# of trees and a data frame containing the y-azimuths of the stands. It returns
# the mapping data frame with new columns for updated X, Y locations of the 
# trees in plots that were reoriented. Reorientation is required for calculation
# of UTM coordinates of individual trees

reorient <- function(tree_x_y, stand_azims){
  
  # Create output table
  output <- tree_x_y %>% filter(TreeID == TreeID[1]) %>% 
    mutate(X_new = NA, Y_new = NA, y_azim = NA)
  
  # Create vector of the stands
  stand_ids <- unique(tree_x_y$StandID)
  
  # Open loop for stands
  for(stand in 1:length(stand_ids)){
    
    # Define action if stand has no y azimuth data
    if(is.na(stand_azims[stand_azims$standid == stand_ids[stand], "y_azim"])) {
      stand_transf <- tree_x_y %>% filter(StandID == stand_ids[stand]) %>%
        mutate(X_new = Xcoord,
               Y_new = Ycoord,
               y_azim = stand_azims[stand_azims$standid == stand_ids[stand], "y_azim"])
      warning(c("Mapping data for ", stand_ids[stand], " was not transformed because no y azimuth information for the stand was provided"))
    } else {
      
      # Determine whether stand needs reorienting (y azimuth greater than 90 degrees)
      if(stand_azims[stand_azims$standid == stand_ids[stand], "y_azim"] > 90) {
        
        # Define action for stands with y azimuths between 90 and 180 degrees 
        if(stand_azims[stand_azims$standid == stand_ids[stand], "y_azim"] <= 180) {
          stand_transf <- tree_x_y %>% filter(StandID == stand_ids[stand]) %>%
            mutate(X_new = Ycoord,
                   Y_new = 100 - Xcoord,
                   y_azim = stand_azims[stand_azims$standid == stand_ids[stand], "y_azim"] - 90)
        } else { 
          
          # Define action for stands with y azimuths between 180 and 270 degrees 
          if(stand_azims[stand_azims$standid == stand_ids[stand], "y_azim"] <= 270) {
            stand_transf <- tree_x_y %>% filter(StandID == stand_ids[stand]) %>%
              mutate(X_new = 100 - Xcoord,
                     Y_new = 100 - Ycoord,
                     y_azim = stand_azims[stand_azims$standid == stand_ids[stand], "y_azim"] - 180)
          } else {
            
            # Define action for stands with y azimuths above 270 degrees
            stand_transf <- tree_x_y %>% filter(StandID == stand_ids[stand]) %>%
              mutate(X_new = 100 - Ycoord,
                     Y_new = Xcoord,
                     y_azim = stand_azims[stand_azims$standid == stand_ids[stand], "y_azim"] - 270)
          }
        }
      } else {
        
        # Define action to be taken if stand does not need reorienting (y azimuth < 90 degrees)
        stand_transf <- tree_x_y %>% filter(StandID == stand_ids[stand]) %>%
          mutate(X_new = Xcoord,
                 Y_new = Ycoord,
                 y_azim = stand_azims[stand_azims$standid == stand_ids[stand], "y_azim"])
      }
    }
    
    # Append reoriented stand to output data frame
    output <- rbind(output, stand_transf)
  }
  
  # Remove first row from output (doesn't contain any information)
  final_result <- output[-1, ]
  
  # Return output data frame
  final_result
}

#=========================================================
# Convert stand corner locations to decimal degrees format
#=========================================================

stand_dms_to_dd <- function(data){
  
  data <- stand_locs_raw %>%
    
    # Convert to long format
    gather("location_id", "deg_min_sec", 3:10) %>%
    
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
      deg_min_sec = ifelse(
        str_detect(deg_min_sec, "W") == TRUE,
        paste("-", deg_min_sec, sep = ""),
        deg_min_sec),
      
      deg_min_sec = ifelse(
        str_detect(deg_min_sec, "S") == TRUE,
        paste("-", deg_min_sec, sep = ""),
        deg_min_sec),
      
      # Remove unwanted symbols from deg_min_sec
      deg_min_sec = str_replace(deg_min_sec, '\"', ""),
      deg_min_sec = str_replace(deg_min_sec, "\'", " "),
      deg_min_sec = str_replace(deg_min_sec, "\\.", " "),
      deg_min_sec = str_sub(deg_min_sec, end = -2),
      
      # Calculate decimal degrees
      dec_deg = as.numeric(conv_unit(deg_min_sec, "deg_min_sec", "dec_deg"))
    ) %>%
    
    # Select required columns
    select(standid, y_azim, corner_id, coord, dec_deg) %>%
    
    # Separate lats and longs
    spread(coord, dec_deg)
}

#===========================
# Create polygons for stands
#===========================

polygons <- function(data){
  
  data %>%
    
    # Group by stand to create multipoint objects
    group_by(standid) %>%
    
    # Need to include summarize for group_by to take effect
    summarise(y_azim = y_azim[1]) %>%
    
    # Use convex hull method to form polygon around each set of four points
    st_convex_hull() %>%
    
    # Return only stand IDs and geometry
    select(standid)
}

#==========================================================
# Plot one stand with trees colored by a specified variable
#==========================================================

utm_mapping <- function(tree_x_y, stand, color_var){
  
  # Isolate focal stand
  focal_stand <- mapping[mapping$StandID == stand, ]
  
  # Calculate degrees to radians conversion factor
  cf <- pi / 180
  
  # Extract y azimuth for focal stand and convert to radians
  y_azim <- stand_locs_raw[stand_locs_raw$standid == stand, "y_azim"]
  y_azim <- y_azim * cf
  
  # Calculate relative azimuths (in radians) from stand origin to each tree
  focal_stand$rel_az <- y_azim + atan(focal_stand$Xcoord / focal_stand$Ycoord)
  
  # Calculate distances from origin to each tree
  focal_stand$dists <- sqrt((focal_stand$Xcoord ^ 2) + (focal_stand$Ycoord ^ 2))
  
  # Extract rows with relative azimuths <= 90 degrees
  sub_90 <- focal_stand[(focal_stand$rel_az / cf) <= 90, ]
  super_90 <- focal_stand[(focal_stand$rel_az / cf) > 90, ]
  
  # Extract X and Y UTM coordinates of stand origin
  origin_x <- stand_utm[stand_utm$stand == stand & 
                          stand_utm$corner_id == "origin", "x"][[1]]
  origin_y <- stand_utm[stand_utm$stand == stand & 
                          stand_utm$corner_id == "origin", "y"][[1]]
  
  # Calculate UTM X and Y for sub 90 group
  sub_90$x_UTM <- origin_x + (sin(sub_90$rel_az) * sub_90$dists)
  sub_90$y_UTM <- origin_y + (cos(sub_90$rel_az) * sub_90$dists)
  
  # Calculate UTM X and Y for super 90 group
  super_90$x_UTM <- origin_x + (sin(super_90$rel_az) * super_90$dists)
  super_90$y_UTM <- origin_y + (cos(super_90$rel_az) * super_90$dists)
  
  # Recombine two groups
  focal_stand_utm <- rbind(sub_90, super_90)
  
  # Transform output to spatial object
  focal_stand_utm <- st_as_sf(focal_stand_utm, coords = c("x_UTM", "y_UTM"), crs = 32610)
  
  # Plot result
  mapview(focal_stand_utm, zcol = color_var) + 
    mapview(stand_polygons[["geometry"]][which(stand_polygons[["standid"]] == stand)])
  
}
