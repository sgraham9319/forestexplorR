
# Check mapping data first

# Create some messy data
map_old <- read.csv("Data/Mapping_2017.csv", stringsAsFactors = F)

# Repeat some rows
map_add <- map_old[sample(1:nrow(map_old), 10), ]

# Combine repeats and original
map_old <- bind_rows(map_old, map_add)

# Give some trees NAs for x, y, or x and y coordinates
map_old[sample(1:nrow(map_old), 5), c("x_coord", "y_coord")] <- NA
map_old[sample(1:nrow(map_old), 5), "x_coord"] <- NA
map_old[sample(1:nrow(map_old), 5), "y_coord"] <- NA

# Give some trees x or y coordinates beyond expected range
# No need to do this because they already exist!

# Give some trees an NA for stand_id
map_old[sample(1:nrow(map_old), 10), "stand_id"] <- NA

# Give some trees an NA for species
map_old[sample(1:nrow(map_old), 10), "species"] <- NA

# Write messy mapping to file
write.csv(map_old, "Data/messy_mapping.csv", row.names = F)



# Check tree data

#----------------------------
# Create some messy tree data
#----------------------------

# Load cleaned tree data
messy_tree <- tree

# Subset to two stands to reduce dataset size
messy_tree <- messy_tree %>%
  filter(stand_id %in% c("AB08", "PP17"))

# Duplicate some trees
set.seed(100)
dups <- messy_tree[sample(1:nrow(messy_tree), 10), ]

# Make changes so all have errors in stand, species, or both
dups$stand_id <- "AB08"
dups$species <- "PSME"

# Rejoin with other tree data
messy_tree <- bind_rows(messy_tree, dups)

# Create some trees without mapping data
set.seed(500)
no_map <- messy_tree[sample(1:nrow(messy_tree), 10), ]
no_map$tree_id <- paste("TEST000", 0:9, sep = "")

# Rejoin with other tree data
messy_tree <- bind_rows(messy_tree, no_map)

# Give some tree records an NA for dbh
messy_tree[sample(1:nrow(messy_tree), 10), "dbh"] <- NA

# Give some tree records an NA for stand_id
messy_tree[sample(1:nrow(messy_tree), 10), "stand_id"] <- NA

# Give some tree records an NA for species
messy_tree[sample(1:nrow(messy_tree), 10), "species"] <- NA

# Give some tree records an NA for year
messy_tree[sample(1:nrow(messy_tree), 10), "year"] <- NA

usethis::use_data(messy_tree)





# Are required columns present? tree_id, stand_id, species, year, dbh

# Tree ids that refer to multiple trees? i.e. duplicate ids

# Trees without mapping data

# Missing dbh data

# Missing stand_id

# Missing species

# Missing year

# Negative growth? No, check in growth calculation function only (add this to
# documention of tree checking function)


# Find trees that have any missing data and highlight them for investigation
# No need to make the call on whether the missing data makes those trees 
# totally unusable


