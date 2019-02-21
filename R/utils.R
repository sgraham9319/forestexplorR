

if(getRversion() >= "2.15.1")  
  utils::globalVariables(c(".", "year", "tree_id", "begin_dif", "end_dif", 
                           "stand_id", "species", "dbh", "annual_growth", 
                           "dbh_diff", "year_diff", "begin_size", "dist", 
                           "abh", "y_azim", "st_convex_hull", "stands_locs_raw", 
                           "gather", "location_id", "deg_min_sec", "conv_unit", 
                           "corner_id", "coord", "dec_deg", "spread", "stand_locs_raw", 
                           "stand_utm", "st_as_sf", "mapview", "stand_polygons"))

#======================================================
# Calculate annual growth over entire measurment period
#======================================================

# Function requires a data frame where separate measurements of the same tree
# appear in different rows. The data frame needs to have a column called "treeid"
# that contains unique tree ID values, a column named "year" containing the year
# of the measurement, and a column named "dbh" containing the dbh measurements.

growth_summary <- function(data){
  data %>% 
    group_by(tree_id) %>% 
    arrange(year) %>%
    summarize(
      stand_id = stand_id[1],
      species = species[1],
      first_record = year[1],
      last_record = year[n()],
      begin_size = dbh[1],
      mean_size = mean(dbh),
      midpoint_size = (min(dbh) + max(dbh)) / 2, 
      annual_growth = (dbh[n()] - dbh[1]) / (year[n()] - year[1])
    ) %>% 
    mutate(size_corr_growth = sqrt(annual_growth / begin_size))
}


#==================================================================
# Calculate annual growth between each pair of consecutive censuses
#==================================================================

detailed_growth <- function(data){
  data %>% 
    group_by(tree_id) %>% 
    arrange(tree_id, year) %>%
    mutate(year_diff = c(diff(year), NA),
           dbh_diff = c(diff(dbh), NA),
           annual_growth = dbh_diff / year_diff) %>%
    select(-year_diff, -dbh_diff)
}


#==============================================
# Calculate annual growth over a defined period
#==============================================

# Function requires a data frame where separate measurements of the same tree
# appear in different rows. The data frame needs to have a column called "tree_id"
# that contains unique tree ID values, a column named "year" containing the year
# of the measurement, and a column named "dbh" containing the dbh measurements. The
# function also requires the user to provide the beginning and end date (as a 4 
# digit year) of the period over which annual growth rates are desired

defined_period_annual_growth <- function(data, begin, end){
  
  # Calculate differences between each measurement and begin and end
  data <- data %>%
    mutate(
      begin_dif = abs(begin - year),
      end_dif = abs(end - year))
  
  # Subset measurements earlier than begin
  before <- data %>%
    filter(year <= begin) %>%
    group_by(tree_id) %>%
    filter(begin_dif == min(begin_dif))
  
  # Subset measurements later than end
  after <- data %>%
    filter(year >= end) %>%
    group_by(tree_id) %>%
    filter(end_dif == min(end_dif))
  
  # Combine subsets
  before_after <- rbind(before, after)
  
  # Calculate annual growth
  period_growth = before_after %>% 
    group_by(tree_id) %>% 
    arrange(year) %>%
    summarize(
      stand_id = stand_id[1],
      species = species[1],
      begin_size = dbh[1],
      annual_growth = (dbh[2] - dbh[1]) / (year[2] - year[1])
    ) %>% 
    filter(!is.na(annual_growth)) %>%
    mutate(size_corr_growth = sqrt(annual_growth / begin_size))
  
  # Return new data frame
  period_growth
}

#======================
# Calculate circle area
#======================

circ_area <- function(radius){
  pi * (radius ^ 2)
}

#=================================================
# Summarize density by species from abh data frame
#=================================================

# Takes a data frame where each row represents a different tree in the stand
# and each column represents a different set of coordinates for which a summary
# of neighborhood density is desired. In each column, the rows representing
# trees that are not in the neighborhood of that tree should contain NA, while
# other cells contain the size of the tree in that row. There is also an 
# additional column for the species identity of the trees. The function returns
# a data frame where each row represents a set of coordinates and each column
# represents a different summary of density at those coordinates

density_calc <- function(data){
  
  # Create list of species
  sps_list <- levels(data$species)
  
  # Create output data frame
  output <- as.data.frame(matrix(nrow = ncol(data) - 1,
                                 ncol = length(sps_list) + 1))
  colnames(output)[1] <- "all_density"
  for(i in 1:length(sps_list)){
    colnames(output)[i + 1] <- paste(sps_list[i], "density", sep = "_")
  }
  
  if(nrow(output) > 1){
    
    # Add density includng all species to output
    output[, "all_density"] <- apply(data[, 2:ncol(data)], 2, sum, na.rm = T)
    
    # Loop through species adding species-specific densities to output
    for(sps in 1:length(sps_list)){
      one_sps <- data %>%
        filter(species == sps_list[sps])
      output[, sps + 1] <- apply(one_sps[, 2:ncol(data)], 2, sum, na.rm = T)
    }
  } else {
    
    # Add density includng all species to output
    output[, "all_density"] <- sum(data[, 2], na.rm = T)
    
    # Loop through species adding species-specific densities to output
    for(sps in 1:length(sps_list)){
      one_sps <- data %>%
        filter(species == sps_list[sps])
      output[, sps + 1] <- sum(one_sps[, 2], na.rm = T)
    }
  }
  
  # Return output
  output
}

#===============================
# Calculate neighborhood density
#===============================

# Takes a set of mapping data and returns a species-specific summary of
# neighborhood density for each tree in the mapping dataset. Required inputs
# are a data frame of mapping data, the name of the stand for which 
# neighborhood densities are to be calculated (to iterate over all stands,
# see "density_all_stands" function which calls "density_summary") entered as
# a character string, and the radius of the desired neighborhood size.

density_summary <- function(mapping, stand, radius){
  
  # Convert species column to factor
  mapping$species <- as.factor(mapping$species)
  
  # Extract relevant coordinates
  coords <- mapping %>%
    filter(stand_id == stand) %>%
    mutate(abh = circ_area(dbh / 2) / circ_area(radius)) %>%
    select(tree_id, species, abh, x_coord, y_coord)
  
  # Create matrix of distances between each tree pair
  dist_mat <- as.matrix(dist(coords[, 4:5], diag = T, upper = T))
  
  # Create matrix where each column contains abh measurements of all trees
  nbhds <- matrix(coords$abh, nrow = nrow(coords), ncol = nrow(coords))
  
  # Change abh values of trees not in each neighborhood to NA
  nbhds[which(dist_mat > radius)] <- NA
  
  # Add species information
  all <-  cbind(coords$species, as.data.frame(nbhds))
  colnames(all)[1] <- "species"
  
  # Summarize density data
  output <- density_calc(all)
  
  # Add tree_id column
  output <- cbind(coords$tree_id, output)
  colnames(output)[1] <- "tree_id"
  
  # Return output
  output
}

#======================================================================
# Calculate neighborhood density for all trees in a multi-stand dataset
#======================================================================

density_all_stands <- function(all_stands, radius){
  
  # Identify unique stands in dataset
  stand_ids <- unique(all_stands$stand_id)
  
  # Loop through stands
  for(stand_num in 1:length(stand_ids)){
    
    # If first stand, make the output table
    if(stand_num == 1){
      
      output <- density_summary(all_stands,
                                stand = stand_ids[stand_num],
                                radius = radius)
    } else { 
      
      # If not first stand, append results to those from earlier stands
      new_stand <- density_summary(all_stands,
                                   stand = stand_ids[stand_num],
                                   radius = radius)
      
      output <- rbind(output, new_stand)
    }
  }
  
  # Return output
  output
}

#=================================================================
# Calculate neighborhood density for any user-provided coordinates
#=================================================================

# Returns neighborhood density for either a single set of coordinates or
# a vector of provided coordinates. By setting the focal_coords argument
# equal to "grid" you can also return estimates of density at every
# intersection of a 5 m grid. Required inputs are a mapping data frame,
# a stand name (as character string), neighborhood radius, and the
# coordinates for which density measurements are desired (entered as a 
# data frame with two columns: x_coord, and y_coord)

density_specific <- function(mapping, stand, radius, focal_coords){
  
  # Extract relevant tree data
  coords <- mapping %>%
    filter(stand_id == stand) %>%
    mutate(abh = circ_area(dbh / 2) / circ_area(radius)) %>%
    select(species, abh, x_coord, y_coord)
  
  if(any(focal_coords == "grid")){
    
    x_coord <- rep(seq(0, 100, 5), each = 21)
    y_coord <- rep(seq(0, 100, 5), times = 21)
    input_coords <- cbind(x_coord, y_coord)
  } else {
    
    input_coords <- focal_coords
  }
  # Create matrix of distances between each focal point and each tree
  all_coords <- rbind(coords[, 3:4], input_coords) 
  dist_mat <- as.matrix(dist(all_coords, diag = T, upper = T))
  dist_mat <- dist_mat[1:nrow(coords), (nrow(coords) + 1):ncol(dist_mat)]
  
  # Create matrix where each column contains abh measurements of all trees
  nbhds <- matrix(coords$abh, nrow = nrow(coords), ncol = nrow(input_coords))
  
  # Change abh values of trees not in each neighborhood to NA
  nbhds[which(dist_mat > radius)] <- NA
  
  # Add species information
  all <-  cbind(coords$species, as.data.frame(nbhds))
  colnames(all)[1] <- "species"
  
  # Summarize density data
  output <- density_calc(all)
  
  # Add input coordinates to output
  output <- cbind(input_coords, output)
  
  # Return output
  output
}
