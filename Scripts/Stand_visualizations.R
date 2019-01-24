# Load required packages
library(measurements)
library(sf)
library(mapview)

# Load package
devtools::load_all()

#=====================================
# Loading and formatting location data
#=====================================

# Load 2013 within-stand location data for individual trees
mapping_raw <- read.csv("../Data/Mapping_2013.csv", stringsAsFactors = F)

# Load location data for stands
stand_locs_raw <- read.csv("../Data/Stand_locations.csv", stringsAsFactors = F)

# Remap individual trees with the western-most corner of each stand as
# the origin
mapping <- reorient(tree_x_y = mapping_raw, stand_azims = stand_locs_raw)

# Convert stand corner locations to decimal degrees
stand_locs <- stand_dms_to_dd(stand_locs_raw)

# Create spatial object of stand locations
stand_locs_sf <- st_as_sf(stand_locs, coords = c("lon_dd", "lat_dd"), crs = 4326)

# Check points appear in expected places
mapview(stand_locs_sf, zcol = "standid")

# Create polygons for each stand using convex hull method
stand_polygons <- polygons(stand_locs_sf)
  
# Plot polygons
mapview(stand_polygons, zcol = "standid")

# Transform coordinate reference system to UTM
stand_utm <- st_transform(stand_locs_sf, crs = 32610)

# Split UTM coordinates into x and y columns 
stand_utm$x <- as.vector(st_coordinates(stand_utm)[,1])
stand_utm$y <- as.vector(st_coordinates(stand_utm)[,2])

#===============================================
# Obtaining UTM coordinates for individual trees
#===============================================

# Load growth data
growth_data <- read.csv("../Data/Tree_growth_2017.csv")

# Remove rows without dbh measurement
growth_data <- growth_data[!is.na(growth_data$dbh), ]

# Find trees that were never measured as 15 cm dbh or larger
small_trees <- growth_data %>% group_by(treeid) %>% 
  summarize(max_size = max(dbh)) %>%
  filter(max_size < 15)

length(unique(growth_data$treeid))


# Exclude rows with trees measured at < 15 cm dbh
growth_data <- growth_data[growth_data$dbh >= 15, ]

# Calculate annual growth over entire measurement period
overall_growth <- overall_annual_growth(growth_data)

# Add growth data to mapping
mapping$ann_growth <- overall_growth$annual_growth[match(mapping$TreeID, overall_growth$treeid)]
mapping$sqrt_ann_growth <- overall_growth$sqrt_annual_growth[match(mapping$TreeID, overall_growth$treeid)]
mapping$size_corr_growth <- overall_growth$size_adj_sqrt_growth[match(mapping$TreeID, overall_growth$treeid)]

hist(mapping$ann_growth)
hist(mapping$sqrt_ann_growth)
hist(mapping$size_corr_growth)

# Plot stand with color representing size corrected growth 
utm_mapping(tree_x_y = mapping, stand = "AV06", color_var = "size_corr_growth")


# Calculate annual growth over a specific time period

# Calculate annual growth between 2005 and 2015
specific_annual_growth <- defined_period_annual_growth(data = growth_data,
                                                       begin = 2005, end = 2015)

# Add growth data to mapping
mapping$sqrt_ann_growth <- specific_annual_growth$sqrt_annual_growth[match(mapping$TreeID, specific_annual_growth$treeid)]

# Plot one stand
utm_mapping(tree_x_y = mapping, stand = "AX15", color_var = "sqrt_ann_growth")
