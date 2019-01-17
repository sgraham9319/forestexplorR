
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
