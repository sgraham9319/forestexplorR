# Load required packages
library(measurements)
library(sf)
library(mapview)
library(tidyverse)

# Load stand location data
stand_locs_raw <- read.csv("../Data/all_stand_location_data.csv", stringsAsFactors = F)

# Calculate decimal degrees and format for analysis
stand_locs <- stand_locs_raw %>%
  
  # Convert to long format
  gather("location_id", "deg_min_sec", 3:10) %>%
  
  # Remove cardinal directions from location IDs
  mutate(corner_id = str_sub(location_id, end = -3),
         
         # Add column for lat or long
         coord = ifelse(
           stringr::str_detect(deg_min_sec, "W") == TRUE,
           "lon",
           "lat"),
         
         # Make deg_min_sec values negative for westings
         deg_min_sec = ifelse(
           str_detect(deg_min_sec, "W") == TRUE,
           paste("-", deg_min_sec, sep = ""),
           deg_min_sec),
         
         # Remove unwanted symbols from deg_min_sec
         deg_min_sec = str_replace(deg_min_sec, '\"', ""),
         deg_min_sec = str_replace(deg_min_sec, "\'", " "),
         deg_min_sec = str_replace(deg_min_sec, "\\.", " "),
         deg_min_sec = str_sub(deg_min_sec, end = -2),
         
         # Calculate decimal degrees
         dec_deg = as.numeric(conv_unit(deg_min_sec, "deg_min_sec", "dec_deg"))
  ) %>%
  
  # Select required columns
  select(stand, y_azim, corner_id, coord, dec_deg) %>%
  
  # Separate lats and longs
  spread(coord, dec_deg)

# Create spatial object
stand_locs_sf <- st_as_sf(stand_locs, coords = c("lon", "lat"), crs = 4326)

# Check points appear in expected places
mapview(stand_locs_sf, zcol="stand")

# Create polygons for each stand using convex hull method
stand_polygons <- stand_locs_sf %>%
  
  # Group by stand to create multipoint objects
  group_by(stand) %>%
  
  # Need to include summarize for group_by to take effect
  summarise(y_azim = y_azim[1]) %>%
  
  # Use convex hull method to form polygon around each set of four points
  st_convex_hull() %>%
  
  # Return only stand IDs and geometry
  select(stand) #get only stand + geometry

# Plot polygons
mapview(stand_polygons, zcol="stand")

# Transform coordinate reference system to UTM
stand_utm <- st_transform(stand_locs_sf, crs = 32610)

# Split UTM coordinates into x and y columns 
stand_utm$x <- as.vector(st_coordinates(stand_utm)[,1])
stand_utm$y <- as.vector(st_coordinates(stand_utm)[,2])

#==============================================
# Reorienting plots for y azimuths < 90 degrees
#==============================================

# Load 2013 mapping data
mapping <- read.csv("../Data/Mapping_2013.csv", stringsAsFactors = F)

# Reload y azimuth data
stand_azimuths <- read.csv("../Data/Stand_locations.csv")

stand <- unique(mapping$StandID)[1]

reorient <- function(tree_x_y, stand_azims, stand){
  
  # Define action for stands with no y azimuth data
  if(is.na(stand_azims[stand_azims$standid == stand, "y_azim"])) {
    stand_transf <- tree_x_y %>% filter(StandID == stand) %>%
      mutate(X_new = Xcoord,
             Y_new = Ycoord,
             y_azim = stand_azims[stand_azims$standid == stand, "y_azim"])
    warning("Mapping data for this stand was not transformed because no y azimuth information for the stand was provided")
  } else {
    
    # Determine whether stand needs reorienting (y azimuth greater than 90 degrees)
    if(stand_azims[stand_azims$standid == stand, "y_azim"] > 90) {
      
      # Define action for stands with y azimuths between 90 and 180 degrees 
      if(stand_azims[stand_azims$standid == stand, "y_azim"] <= 180) {
        stand_transf <- tree_x_y %>% filter(StandID == stand) %>%
          mutate(X_new = Ycoord,
                 Y_new = 100 - Xcoord,
                 y_azim = stand_azims[stand_azims$standid == stand, "y_azim"] - 90)
      } else { 
        
        # Define action for stands with y azimuths between 180 and 270 degrees 
        if(stand_azims[stand_azims$standid == stand, "y_azim"] <= 270) {
          stand_transf <- tree_x_y %>% filter(StandID == stand) %>%
            mutate(X_new = 100 - Xcoord,
                   Y_new = 100 - Ycoord,
                   y_azim = stand_azims[stand_azims$standid == stand, "y_azim"] - 180)
        } else {
          
          # Define action for stands with y azimuths above 270 degrees
          stand_transf <- tree_x_y %>% filter(StandID == stand) %>%
            mutate(X_new = 100 - Ycoord,
                   Y_new = Xcoord,
                   y_azim = stand_azims[stand_azims$standid == stand, "y_azim"] - 270)
        }
      }
    } else {
      
      # Define action to be taken if stand does not need reorienting
      stand_transf <- tree_x_y %>% filter(StandID == stand) %>%
        mutate(X_new = Xcoord,
               Y_new = Ycoord,
               y_azim = stand_azims[stand_azims$standid == stand, "y_azim"])
    }
    
  }
  # Return new data frame
  stand_transf
}

test <- reorient(tree_x_y = mapping, stand_azims = stand_azimuths, stand = "AV02")
      


#===============================================
# Obtaining UTM coordinates for individual trees
#===============================================

# Isolate AG05 stand
ag05 <- mapping[mapping$StandID == "AG05", ]

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
origin_x <- stand_utm[stand_utm$stand == "AG05" & 
                        stand_utm$corner_id == "origin", "x"][[1]]
origin_y <- stand_utm[stand_utm$stand == "AG05" & 
                        stand_utm$corner_id == "origin", "y"][[1]]

# Calculate UTM X and Y for sub 90 group
sub_90$x_UTM <- origin_x + (sin(sub_90$rel_az) * sub_90$dists)
sub_90$y_UTM <- origin_y + (cos(sub_90$rel_az) * sub_90$dists)

# Calculate UTM X and Y for super 90 group
super_90$x_UTM <- origin_x + (sin(super_90$rel_az) * super_90$dists)
super_90$y_UTM <- origin_y + (cos(super_90$rel_az) * super_90$dists)

# Recombine two groups
ag05_utm <- rbind(sub_90, super_90)

# Transform to spatial object
ag05_utm <- st_as_sf(ag05_utm, coords = c("x_UTM", "y_UTM"), crs = 32610)

# Plot trees
mapview(ag05_utm, zcol = "Species") + 
  mapview(stand_polygons[stand_polygons$stand == "AG05",], legend = F)
