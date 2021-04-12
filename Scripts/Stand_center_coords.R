# Load TreeNeighborhood
devtools::load_all()

# Load required packages
library(tidyr)
library(measurements)
library(sf)
library(mapview)

# Load stand location data
stand_locs_raw <- read.csv("Data/Stand_locations.csv", stringsAsFactors = F)

# Convert stand corner locations to decimal degrees
stand_locs <- stand_dms_to_dd(stand_locs_raw)

# Create spatial object of stand locations
stand_locs_sf <- st_as_sf(stand_locs, coords = c("lon_dd", "lat_dd"), crs = 4326)

# Transform coordinate reference system to UTM
stand_utm <- st_transform(stand_locs_sf, crs = 32610)

# Split UTM coordinates into x and y columns 
stand_utm$x <- as.vector(st_coordinates(stand_utm)[,1])
stand_utm$y <- as.vector(st_coordinates(stand_utm)[,2])

# Create data frame for stand center UTM locations
stand_centers <- as.data.frame(matrix(NA, nrow = length(unique(stand_locs_raw$stand_id)), ncol = 4))
names(stand_centers) <- c("stand_id", "y_azim", "origin_x", "origin_y")

# Calculate degrees to radians conversion factor
cf <- pi / 180

# Identify stands
stands <- unique(stand_locs_raw$stand_id)

# Add stand id, y azimuth, and origin coordinates of all stands
for(stand in 1:length(stands)){
  stand_centers[stand, "stand_id"] <- stands[stand]
  stand_centers[stand, "y_azim"] <- stand_locs_raw[stand_locs_raw$stand_id == stands[stand], "y_azim"] * cf
  stand_centers[stand, "origin_x"] <- stand_utm[stand_utm$stand == stands[stand] & 
                                                  stand_utm$corner_id == "origin", "x"][[1]]
  stand_centers[stand, "origin_y"] <- stand_utm[stand_utm$stand == stands[stand] & 
                                                  stand_utm$corner_id == "origin", "y"][[1]]
}

# Calculate relative azimuth (in radians) from origin to center for each stand
stand_centers$rel_az <- stand_centers$y_azim + atan(50 / 50)

# Calculate distance from origin to center (same for all stands)
dist <- sqrt((50 ^ 2) + (50 ^ 2))

# Calculate UTM coordinates of stand centers
stand_centers$center_x <- stand_centers$origin_x + (sin(stand_centers$rel_az) * dist)
stand_centers$center_y <- stand_centers$origin_y + (cos(stand_centers$rel_az) * dist)

# Transform output to spatial object
center_utm <- st_as_sf(stand_centers, coords = c("center_x", "center_y"), crs = 32610)

# Create polygons for each stand using convex hull method
stand_polygons <- polygons(stand_locs_sf)

# Plot result
mapview(center_utm) + 
  mapview(stand_polygons[["geometry"]])

# Save as shape file
st_write(center_utm, "Data/stand_centers.shp")
