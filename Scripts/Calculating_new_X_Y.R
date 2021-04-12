# Load package
devtools::load_all()

# Load new mapping data
new_map <- read.csv("Data/new_mapping_2017.csv", stringsAsFactors = F)

# Load old mapping data
old_map <- read.csv("Data/Mapping_2013.csv", stringsAsFactors = F)

# Load location data for stands
stand_locs <- read.csv("Data/Stand_locations.csv", stringsAsFactors = F)

# Check for repeated tree ids in mapping data
length(unique(new_map$tree_id)) == nrow(new_map)
length(unique(old_map$tree_id)) == nrow(old_map) # repeats in old data

# Identify repeats in old mapping
repeats <- names(which(table(old_map$tree_id) > 1))

# View repeats
View(old_map[old_map$tree_id %in% repeats, ])

# Both repeats appear to be cases of double mapping with very similar x and y
# coordinates for the two mapping events. Average the values
average_trees <- old_map[old_map$tree_id %in% repeats, ]
average_trees <- average_trees %>%
  group_by(tree_id) %>%
  summarize(psp_study_id = psp_study_id[1], stand_id = stand_id[1], plot = plot[1],
            tag = tag[1], species = species[1], year = year[1], tree_status = tree_status[1], 
            live_dead = live_dead[1], mort_year = mort_year[1], dbh = dbh[1],
            x_coord = round(mean(x_coord), 2) , y_coord = round(mean(y_coord), 2))

# Replace repeats with averages in old mapping
old_map <- old_map[old_map$tree_id %in% repeats == F, ]
old_map <- rbind(old_map, average_trees)

# Recheck for repeated mapping data
length(unique(old_map$tree_id)) == nrow(old_map)

# Calculate degrees to radians conversion factor
cf <- pi / 180

# Join old_map to new_map using stand id and tag number columns
combined_map <- new_map %>%
  left_join(old_map, by = c("stand_id", "ref_tree" = "tag")) %>%
  
  # Reduce dataset to desired columns
  select(tree_id.x, stand_id, plot.x, species.x, tag, year.x, tree_status.x, 
         dbh.x, ref_tree, azim, dist,
         tree_id.y, dbh.y, x_coord, y_coord) %>%

  # Rename columns
  rename(tree_id = tree_id.x,
         plot = plot.x,
         species = species.x,
         year = year.x,
         tree_status = tree_status.x,
         dbh = dbh.x,
         ref_id = tree_id.y,
         ref_dbh = dbh.y,
         ref_x = x_coord,
         ref_y = y_coord) %>%
  
  # Add stand location data
  left_join(stand_locs[, c("stand_id", "y_azim")],
            by = "stand_id") %>%
  
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

# View these cases
View(combined_map[combined_map$tree_id %in% unmappable_trees, ])

# This single case matches two totally different trees, so must remove it
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
write.csv(failed_update_data, "Data/Unusable_new_mapping.csv", row.names = F)

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

# Check for any tree IDs that appear in dataset multiple times
length(unique(updated_map$tree_id)) == nrow(updated_map)
# No more repeats present

# Write updated mapping to file
write.csv(updated_map, "Data/Mapping_2017.csv", row.names = F)
