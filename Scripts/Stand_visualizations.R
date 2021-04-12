# Load required packages
library(measurements)
library(sf)
library(mapview)
library(plotly)
library(tidyverse)

# Load package
devtools::load_all()

#=====================================
# Loading and formatting location data
#=====================================

# Load 2013 within-stand location data for individual trees
mapping_raw <- read.csv("Data/Mapping_2017.csv", stringsAsFactors = F)

# Load location data for stands
stand_locs_raw <- read.csv("Data/Stand_locations.csv", stringsAsFactors = F)

# Remap individual trees with the western-most corner of each stand as
# the origin
mapping <- mapping_raw

# Convert stand corner locations to decimal degrees
stand_locs <- stand_dms_to_dd(stand_locs_raw)

# Create spatial object of stand locations
stand_locs_sf <- st_as_sf(stand_locs, coords = c("lon_dd", "lat_dd"), crs = 4326)

# Check points appear in expected places
mapview(stand_locs_sf, zcol = "stand_id")

# Create polygons for each stand using convex hull method
stand_polygons <- polygons(stand_locs_sf) # message is not an error
  
# Plot polygons
mapview(stand_polygons, zcol = "stand_id")

# Transform coordinate reference system to UTM
stand_utm <- st_transform(stand_locs_sf, crs = 32610)

# Split UTM coordinates into x and y columns 
stand_utm$x <- as.vector(st_coordinates(stand_utm)[,1])
stand_utm$y <- as.vector(st_coordinates(stand_utm)[,2])

#===============================================
# Obtaining UTM coordinates for individual trees
#===============================================

# Load growth data
growth_data <- read.csv("Data/Tree_growth_2017.csv", stringsAsFactors = F)

# Remove rows without dbh measurement
growth_data <- growth_data[!is.na(growth_data$dbh), ]

# Find trees that were never measured as 15 cm dbh or larger
small_trees_data <- growth_data %>% group_by(tree_id) %>% 
  summarize(max_size = max(dbh)) %>%
  filter(max_size < 15)

# Exclude these trees
small_trees <- unique(small_trees_data$tree_id)
growth_data <- growth_data[-which(growth_data$tree_id %in% small_trees), ]

# Exclude small trees from mapping data
mapping <- mapping[-which(mapping$tree_id %in% small_trees), ]

# Calculate annual growth over entire measurement period
overall_growth <- growth_summary(growth_data)

# Add growth data to mapping
mapping$ann_growth <- overall_growth$annual_growth[match(mapping$tree_id, overall_growth$tree_id)]
mapping$size_corr_growth <- overall_growth$size_corr_growth[match(mapping$tree_id, overall_growth$tree_id)]

# Check distributions of growth data
hist(mapping$ann_growth)
hist(mapping$size_corr_growth)

# Plot a stand with color representing sqrt(annual growth / initial size)
utm_mapping(tree_x_y = mapping, stand = "PP17", color_var = "species")

#=====================================
# Calculating density around each tree
#=====================================

# Calculate density at local coordinates 50, 50 in stand AB08 given
# neighborhood radius of 10 m
nbhd_density(mapping_data = mapping_raw, stand = "AB08", x = 50, y = 50, 
             nbhd_radius = 10)

# Calculate density at every intersection of a 5 x 5 grid across AB08 given
# neighborhood radius of 10 m
des_x <- rep(seq(0, 100, 5), each = 21)
des_y <- rep(seq(0, 100, 5), times = 21)
result <- nbhd_density(mapping_data = mapping_raw, stand = "AB08", x = des_x, 
                       y = des_y, nbhd_radius = 10)
result <- density_specific(mapping_raw, "AB08", 10, "grid")

# Create a contour plot of density for AB08
result <- density_specific(mapping_raw, "AB08", 10, "grid")
plot_ly(
  x = result$x_coord, 
  y = result$y_coord, 
  z = result$all_density, 
  type = "contour" 
)

result1 <- result
result1[, 3:ncol(result1)] <- (result1[, 3:ncol(result1)] / 10000) * (10000 / circ_area(10))
plot_ly(
  x = a[, "x_coord"], 
  y = a[, "y_coord"], 
  z = a[, "all_density"], 
  type = "contour" 
)

# Calculate neighborhood density for each tree in stand AB08
ab08 <- mapping[mapping$stand_id == "AB08", ]

ab08_density <- nbhd_density(mapping_data = ab08, stand = "AB08", x = ab08$x_coord, 
                       y = ab08$y_coord, nbhd_radius = 10)
ab08_density <- cbind(ab08$tree_id, ab08_density)
colnames(ab08_density)[1] <- "tree_id"

# Combine density with annual growth
ab08_density$ann_growth <- 
  overall_growth$size_corr_growth[match(ab08_density$tree_id,
                                            overall_growth$tree_id)]

# Plot annual growth against density
plot(ab08_density$ann_growth ~ ab08_density$all_density)
tshe_ids <- mapping[mapping$stand_id == "AB08" & mapping$species == "TSHE",
        "tree_id"]
length(tshe_ids)
ab08_tshe_density <- ab08_density[ab08_density$tree_id %in% tshe_ids, ]
plot(ab08_tshe_density$ann_growth ~ ab08_tshe_density$all_density)
summary(lm(ann_growth ~ all_density, data = ab08_tshe_density))

# Repeat for AV06
av06 <- mapping[mapping$stand_id == "AV06", ]
av06_density <- nbhd_density(mapping_data = av06, stand = "AV06", x = av06$x_coord, 
                             y = av06$y_coord, nbhd_radius = 10)
av06_density <- cbind(av06$tree_id, av06_density)
colnames(av06_density)[1] <- "tree_id"
av06_density$ann_growth <- 
  overall_growth$size_corr_growth[match(av06_density$tree_id,
                                            overall_growth$tree_id)]
plot(av06_density$ann_growth ~ av06_density$all_density)
tshe_ids <- mapping[mapping$stand_id == "AV06" & mapping$species == "TSHE",
                    "tree_id"]
length(tshe_ids)
av06_tshe_density <- av06_density[av06_density$tree_id %in% tshe_ids, ]
plot(av06_tshe_density$ann_growth ~ av06_tshe_density$all_density)
summary(lm(ann_growth ~ all_density, data = av06_tshe_density))
