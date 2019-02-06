# Load package
devtools::load_all()

# Load new mapping data
new_map <- read.csv("../Data/new_mapping_2017.csv", stringsAsFactors = F)

# Load old mapping data
old_map <- read.csv("../Data/Mapping_2013.csv", stringsAsFactors = F)

# Load location data for stands
stand_locs <- read.csv("../Data/Stand_locations.csv", stringsAsFactors = F)

# Calculate degrees to radians conversion factor
cf <- pi / 180

# Join old_map to new_map using stand id and tag number columns
combined_map <- new_map %>%
  left_join(old_map, by = c("stand_id", "ref_tree" = "tag")) %>%
  
  # Reduce dataset to desired columns
  select(tree_id.x, stand_id, plot.x, tag, year.x, tree_status.x, 
         dbh.x, ref_tree, azim, dist,
         tree_id.y, dbh.y, x_coord, y_coord) %>%

  # Rename columns
  rename(tree_id = tree_id.x,
         plot = plot.x,
         year = year.x,
         tree_status = tree_status.x,
         dbh = dbh.x,
         ref_id = tree_id.y,
         ref_dbh = dbh.y,
         ref_x = x_coord,
         ref_y = y_coord) %>%
  
  # Add stand location data
  left_join(stand_locs[, c("standid", "y_azim")],
            by = c("stand_id" = "standid")) %>%
  
  # Calculate relative azimuths
  mutate(rel_azim = azim - y_azim,
         
         # Correct distance from reference tree using dbh measurements
         tot_dist = dist + ((dbh + ref_dbh) / 200),
         
         # Calculate X, Y coordinates of newly mapped trees
         x_coord = ref_x + (tot_dist * sin(cf * rel_azim)),
         y_coord = ref_y + (tot_dist * cos(cf * rel_azim)))

# Identify trees that are unmappable because their reference tree tag number
# matches multiple trees in their stand
unmappable_trees <- names(which(table(combined_map$tree_id) > 1))

# Remove these trees from dataset
combined_map <- combined_map[-which(combined_map$tree_id %in% unmappable_trees), ]

# Create data frame of trees in new mapping data that could not be mapped
failed_update_ids <- combined_map[which(is.na(combined_map$x_coord)), "tree_id"]
failed_update_data <- new_map[new_map$tree_id %in% failed_update_ids |
                                new_map$tree_id %in% unmappable_trees, ]

# Add column indicating trees for which reference tree did not exist versus
# trees for which the reference tag was associated with multiple tree IDs in
# the stand
failed_update_data$ref_tree_problem <- "Does not exist"
failed_update_data$ref_tree_problem[which(failed_update_data$tree_id %in% 
                                            unmappable_trees)] <- "Ref tag duplicated"

# Write failed mapping update data to file
write.csv(failed_update_data, "../Data/Unusable_new_mapping.csv", row.names = F)

# Remove columns in combined_map that are no longer needed
combined_map <- combined_map %>%
  select(-ref_tree, -azim, -dist, -ref_id, -ref_dbh, -ref_x, -ref_y,
         -y_azim, -rel_azim, -tot_dist)

# Remove columns that appear in old_map but not combined_map
old_map <- old_map[, which(names(old_map) %in% names(combined_map))]

# Remove from old_map the trees that need their mapping updated
old_map <- old_map[-which(old_map$tree_id %in% combined_map$tree_id), ]

# Combine old_map and combined_map to get final new mapping
updated_map <- rbind(old_map, combined_map)

# Remove failed update data from the updated mapping data frame
updated_map <- updated_map[-which(updated_map$tree_id %in% failed_update_ids), ]

# Write updated mapping to file
write.csv(updated_map, "../Data/Mapping_2017.csv", row.names = F)
