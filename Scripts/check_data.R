
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


# User defined parameters
map_data <- map_old
#map_data <- map_old %>% select(-species)
max_x <- 100
max_y <- 100
test <- mapping_check(map_old, 100, 100)

