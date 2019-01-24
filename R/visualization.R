
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
  stand_ids <- unique(mapping$StandID)
  
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