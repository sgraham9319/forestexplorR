##############################
# Data cleaning and formatting
##############################

# Author: Stuart Graham
# Created: 4/8/2020
# Last edited: 4/8/2020

# Load required packages
library(dplyr)

# Load raw mapping data
mapping <- read.csv("Data/Mapping_2017.csv", stringsAsFactors = F)

# Load raw tree measurement data
tree <- read.csv("Data/Tree_growth_2017.csv", stringsAsFactors = F)

#==============
# Data cleaning
#==============

# Trees with no mapping data are not useful but trees for which we have size data but
# no growth data (because measured only once) can still be used as competitors (but
# not focals). Therefore, we will clean the mapping data first, then remove trees
# with no mapping from the tree data before cleaning the tree data

# Q1. Are there any tree ID numbers that refer to multiple trees?
length(unique(mapping$tree_id)) == nrow(mapping)
# All trees in mapping data have unique ID numbers

# Q2. Do all trees have a stand_id, species, dbh, and x and y coordinates?
sum(!is.na(mapping$stand_id)) == nrow(mapping)
sum(!is.na(mapping$species)) == nrow(mapping)
sum(!is.na(mapping$dbh)) == nrow(mapping)
sum(!is.na(mapping$x_coord)) == nrow(mapping)
sum(!is.na(mapping$y_coord)) == nrow(mapping)
# All trees in mapping data have a stand_id, species, dbh, and x and y coordinates

# Q3. Are any x or y coordinates illogical?
max(mapping$x_coord)
min(mapping$x_coord)
max(mapping$y_coord)
min(mapping$y_coord)
# All x and y coordinates seem reasonable (OK if slightly beyond stand boundary)

# Q4. Are there any tree ID numbers in tree data that refer to multiple trees?
tree_ids_summ <- tree %>%
  group_by(tree_id) %>%
  summarize(num_stands = length(unique(stand_id)), num_species = length(unique(species)))
max(tree_ids_summ$num_species)
max(tree_ids_summ$num_stands)
# All tree ids in tree data refer to a single species in a single stand so we
# can conclude that all refer to distinct trees

# Q5. How many trees in tree data do not have mapping data?
tree_data_ids <- unique(tree$tree_id)
length(tree_data_ids) - sum(tree_data_ids %in% mapping$tree_id)
# 86 tree_ids in tree data (<1% of trees) do not have mapping data.
# This will add noise to the data unfortunately, but they must be
# dropped

# Drop trees without mapping data from tree data
tree <- tree[tree$tree_id %in% mapping$tree_id, ]

# Q6. Are there trees with missing dbh data?
sum(is.na(tree$dbh))
# Yes, 2053 occurences of missing dbh data. Why is data missing?
no_dbh <- tree[is.na(tree$dbh), ]
table(no_dbh$tree_status)
# 2033 were dead (code 6) and 20 were not found (code 9)

# Q7. Do we have any dbh data for trees reported missing?
missing_ids <- no_dbh[no_dbh$tree_status == 9, "tree_id"]
missing <- tree %>%
  filter(tree_id %in% missing_ids) %>%
  group_by(tree_id) %>%
  summarize(dbh_vals = sum(!is.na(dbh)))
min(missing$dbh_vals)
# We have at least 1 dbh value for each tree reported missing so all 
# can be used as competitors still. Good to know we aren't losing 
# competitors from our data by excluding these. Go ahead and remove
# all cases of NA dbh values
tree <- tree[!is.na(tree$dbh), ]

# Q8. Do all trees in tree data have species information?
sum(!is.na(tree$species)) == nrow(tree)
# Yes, all trees have species information

# Q9. Are there any trees that were never observed to have dbh of 15cm or greater?
small_trees_data <- tree %>%
  group_by(tree_id) %>% 
  summarize(max_size = max(dbh), stand = stand_id[1]) %>%
  filter(max_size < 15)
# Yes, there are many (2577) because of the detail plots where 5 cm was the
# minimum dbh.

# Q10. Are there detail plots in every stand?
table(small_trees_data$stand)
# Yes, there is a good number of small trees in each stand

#====================
# Saving cleaned data
#====================

# Before saving to file, remove trees from mapping data that don't occur
# in tree data
mapping <- mapping[mapping$tree_id %in% tree$tree_id, ]

# Identify mapped trees that were never measured as > 15cm dbh
mapping$size_cat <- "regular"
mapping$size_cat[mapping$tree_id %in% small_trees_data$tree_id] <- "small"

# Save cleaned mapping data
write.csv(mapping, "Data/Cleaned_mapping_2017.csv", row.names = F)

# Save cleaned tree data
write.csv(tree, "Data/Cleaned_tree_growth_2017.csv", row.names = F)
