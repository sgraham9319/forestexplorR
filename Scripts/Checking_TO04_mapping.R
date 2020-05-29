

# Load TO04 data
tO04_growth <- read.csv("../Data/Tree_growth_2017_TO04.csv", stringsAsFactors = F)

# Load mapping 2013 data
old_map <- read.csv("../Data/Mapping_2013.csv", stringsAsFactors = F)

# Get unique tree ids from TO04
TO04_ids <- unique(tO04_growth$treeid)

# Subset mapping data to trees matching ids from TO04
TO04_old_map <- old_map[old_map$tree_id %in% TO04_ids, ]

# Check these ids are all actually in TO04
sum(TO04_old_map$stand_id == "TO04") == nrow(TO04_old_map)

# Find tree ids from TO04 not in old mapping data
TO04_miss_map <- TO04_ids[TO04_ids %in% TO04_old_map$tree_id == F]

# Add new mapping column to TO04 growth data
tO04_growth$new_mapping <- ""

# Indicate which trees are missing mapping data  
tO04_growth$new_mapping[which(tO04_growth$treeid %in% TO04_miss_map)] <- "Y"

# Load growth data missing TO04
growth_rest <- read.csv("../Data/Tree_growth_2017_no_TO04.csv", stringsAsFactors = F)

# Compare column names between two growth datasets
names(tO04_growth)
names(growth_rest)
# Slight differences in names but columns match

# Change names of TO04 dataset to match
names(tO04_growth) <- names(growth_rest)

# Bind the two datasets together
all_growth <- rbind(growth_rest, tO04_growth)

# Write out complete growth dataset
write.csv(all_growth, "../Data/Tree_growth_2017.csv", row.names = F)

# Subset TO04 growth data to only those missing from mapping data
TO04_new_map <- tO04_growth[tO04_growth$new_mapping == "Y", ]

# Load new 2017 mapping data for all stands but TO04
new_map_2017 <- read.csv("../Data/new_mapping_2017_no_TO04.csv", stringsAsFactors = F)

# Compare names of two new mapping datasets
names(TO04_new_map)
names(new_map_2017)

# Update columns of TO04 mapping to match new mapping
names(TO04_new_map)[22] <- "sample_date"
TO04_new_map <- TO04_new_map[, which(names(TO04_new_map) %in% names(new_map_2017))]
TO04_new_map$ref_tree <- NA
TO04_new_map$azim <- NA
TO04_new_map$dist <- NA

# Manually add reference tree, azim and dist information to TO04 mapping by reference
# to check_notes column - only do one row for tree ids appearing multiple times
TO04_new_map[1, 14:16] <- c(4392, 18, 3.6)
TO04_new_map[2, 14:16] <- c(4390, 175, 4.1)
TO04_new_map[6, 14:16] <- c(4442, 64, 3.15)
TO04_new_map[7, 14:16] <- c(2524, 110, 2.8)
TO04_new_map[8, 14:16] <- c(335, 13, 3.09)
TO04_new_map[9, 14:16] <- c(435, 210, 5.07)
TO04_new_map[10, 14:16] <- c(4450, 174, 4.98)
TO04_new_map[11, 14:16] <- c(438, 225, 4.22)
TO04_new_map[12, 14:16] <- c(4450, 123, 8.04)
TO04_new_map[14, 14:16] <- c(3647, 275, 4.38)
TO04_new_map[15, 14:16] <- c(2835, 13, 2.51)
TO04_new_map[16, 14:16] <- c(3645, 255, 5.00)
TO04_new_map[17, 14:16] <- c(15, 124, 0.33)
TO04_new_map[18, 14:16] <- c(3645, 216, 4.34)
TO04_new_map[22, 14:16] <- c(16, 175, 2.37)
TO04_new_map[25, 14:16] <- c(3640, 141, 1.51)
TO04_new_map[28, 14:16] <- c(4663, 316, 5.2)
TO04_new_map[29, 14:16] <- c(4663, 294, 2.6)
TO04_new_map[30, 14:16] <- c(4685, 33, 2.9)
TO04_new_map[31, 14:16] <- c(3637, 344, 1.9)
TO04_new_map[32, 14:16] <- c(2546, 12, 2.3)
TO04_new_map[33, 14:16] <- c(4692, 63, 3.5)
TO04_new_map[34, 14:16] <- c(4692, 75, 4.8)
TO04_new_map[35, 14:16] <- c(4515, 145, 5.3)

# Subset TO04 to trees with usuable mapping data
new_TO04 <- TO04_new_map[!is.na(TO04_new_map$azim), ]

# Combine with other mapping
all_new_map <- rbind(new_map_2017, new_TO04)

# Write output
write.csv(all_new_map, "../Data/new_mapping_2017.csv", row.names = F)
