# Load required packages
library(measurements)
library(sf)
library(mapview)
library(plotly)

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
growth_data <- read.csv("../Data/Tree_growth_2017.csv", stringsAsFactors = F)

# Remove rows without dbh measurement
growth_data <- growth_data[!is.na(growth_data$dbh), ]

# Find trees that were never measured as 15 cm dbh or larger
small_trees_data <- growth_data %>% group_by(treeid) %>% 
  summarize(max_size = max(dbh)) %>%
  filter(max_size < 15)

# Exclude these trees
small_trees <- unique(small_trees_data$treeid)
growth_data <- growth_data[-which(growth_data$treeid %in% small_trees), ]

# Exclude small trees from mapping data
mapping <- mapping[-which(mapping$TreeID %in% small_trees), ]

# Calculate annual growth over entire measurement period
overall_growth <- overall_annual_growth(growth_data)

# Add growth data to mapping
mapping$ann_growth <- overall_growth$annual_growth[match(mapping$TreeID, overall_growth$treeid)]
mapping$sqrt_ann_growth <- overall_growth$sqrt_annual_growth[match(mapping$TreeID, overall_growth$treeid)]
mapping$size_corr_growth <- overall_growth$size_adj_sqrt_growth[match(mapping$TreeID, overall_growth$treeid)]

# Check distributions of growth data
hist(mapping$ann_growth)
hist(mapping$sqrt_ann_growth)
hist(mapping$size_corr_growth)

# Plot a stand with color representing sqrt(annual growth / initial size)
utm_mapping(tree_x_y = mapping, stand = "AB08", color_var = "size_corr_growth")

utm_mapping1(tree_x_y = mapping, stand = "AB08", color_var = "size_corr_growth")


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

# Create a contour plot of density for AB08
plot_ly(
  x = result$x_coord, 
  y = result$y_coord, 
  z = result$total_abh, 
  type = "contour" 
)





