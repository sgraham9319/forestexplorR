
# Load required packages
library(measurements)
library(sp)
library(sf)
library(mapview)
library(rgdal)

# Load stand location data
stand <- read.csv("../Data/Stand_locations.csv")

# Change degrees/minutes/seconds to decimal degrees
stand$northing <- as.character(stand$northing)
stand$westing <- as.character(stand$westing)
stand$lat <- as.numeric(conv_unit(stand$northing, "deg_min_sec", "dec_deg"))
stand$lon <- as.numeric(conv_unit(stand$westing, "deg_min_sec", "dec_deg"))
stand$lon <- -stand$lon # longitudes should be negative

# Transform to spatial object
coordinates(stand) <- ~ lon + lat

# Plot points
plot(stand)

# Check coordinate reference system (CRS)
proj4string(stand) # right now is NA, need to set

# Assign WGS84 EPSG code to data
proj4string(stand) <- CRS("+init=epsg:4326")

# Transform to UTM Zone 10
stand_UTM <- spTransform(stand, CRS("+init=epsg:32610"))

# Double check in mapview that coordinates are correct
mapview(stand_UTM)

# Now convert to sf object (easier to work with, issues before with starting as sf which I didn't have time to trouble shoot so kept as sp)
new_stand_UTM <- st_as_sf(stand_UTM)

# Get x and y columns for each point, so you can add meters 
new_stand_UTM$x <- as.vector(st_coordinates(new_stand_UTM)[,1])
new_stand_UTM$y <- as.vector(st_coordinates(new_stand_UTM)[,2])

# Plot stands
mapview(new_stand_UTM)

#===============================================
# Obtaining UTM coordinates for individual trees
#===============================================

# Load 2013 mapping data
mapping <- read.csv("../Data/Mapping_2013.csv")

# Isolate AG05 stand
ag05 <- mapping[mapping$StandID == "AG05", ]

# Reload y azimuth data
stand_azimuths <- read.csv("../Data/Stand_locations.csv")

# Calculate degrees to radians conversion factor
cf <- pi / 180

# Extract y azimuth for AG05 and convert to radians
y_azim <- stand_azimuths[stand_azimuths$standid == "AG05", "y_azim"]
y_azim <- y_azim * cf

# Calculate relative azimuths (in radians) from stand origin to each tree
ag05$rel_az <- y_azim + atan(ag05$Xcoord / ag05$Ycoord)

# Calculate distances from origin to each tree
ag05$dists <- sqrt((ag05$Xcoord ^ 2) + (ag05$Ycoord ^ 2))

# Extract rows with relative azimuths <= 90 degrees
sub_90 <- ag05[(ag05$rel_az / cf) <= 90, ]
super_90 <- ag05[(ag05$rel_az / cf) > 90, ]

# Extract X and Y UTM coordinates of stand origin
origin_x <- new_stand_UTM[new_stand_UTM$standid == "AG05", "x"][[1]]
origin_y <- new_stand_UTM[new_stand_UTM$standid == "AG05", "y"][[1]]

# Calculate UTM X and Y for sub 90 group
sub_90$x_UTM <- origin_x + (sin(sub_90$rel_az) * sub_90$dists)
sub_90$y_UTM <- origin_y + (cos(sub_90$rel_az) * sub_90$dists)

# Calculate UTM X and Y for sub 90 group
super_90$x_UTM <- origin_x + (sin(super_90$rel_az) * super_90$dists)
super_90$y_UTM <- origin_y - (cos(super_90$rel_az) * super_90$dists)

# Recombine two groups
ag05_UTM <- rbind(sub_90, super_90)

# Transform to spatial object
coordinates(ag05_UTM) <- ~ x_UTM + y_UTM

# Assign UTM Zone 10 to data
proj4string(ag05_UTM) <- CRS("+init=epsg:32610")

# Plot trees
mapview(ag05_UTM)
