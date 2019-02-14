

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

#=======================
# Calculate tree density
#=======================

# This function calculates tree density as both the number of trees and the
# summed area at breast height of trees within a user-defined distance of a
# provided set of coordinates.

nbhd_density <- function(mapping_data, stand, x, y, nbhd_radius){
  
  # Subset to stand of interest
  focal_stand <- mapping_data %>%
    filter(stand_id == stand)
  
  # Create output vectors
  x_coord <- c()
  y_coord <- c()
  num_trees <- c()
  all_density <- c()
  tshe_density <- c()
  abam_density <- c()
  thpl_density <- c()
  tsme_density <- c()
  cano_density <- c()
  pico_density <- c()
  psme_density <- c()
  
  # Begin looping through input coordinates
  for(coord_num in 1:length(x)){
    
    focal_stand_summary <- focal_stand %>%
      
      # Calculate distance of each tree from input coordinates
      mutate(dist = sqrt((abs(x_coord - x[coord_num]) ^ 2) +
                           (abs(y_coord - y[coord_num]) ^ 2))) %>%
      
      # Subset to trees within input radius of input coordinates
      filter(dist <= nbhd_radius) %>%
      
      # Calculate cross-sectional area at breast height (abh) for each tree
      mutate(abh = pi * ((dbh / 2) ^ 2)) %>%
      
      # Extract number of trees and total abh within neighborhood
      summarize(x_coord = x[coord_num],
                y_coord = y[coord_num],
                num_trees = n(),
                all_density = sum(abh) / (pi * (nbhd_radius ^ 2)),
                tshe_density = sum(abh[which(species == "TSHE")]) / (pi * (nbhd_radius ^ 2)),
                abam_density = sum(abh[which(species == "ABAM")]) / (pi * (nbhd_radius ^ 2)),
                thpl_density = sum(abh[which(species == "THPL")]) / (pi * (nbhd_radius ^ 2)),
                tsme_density = sum(abh[which(species == "TSME")]) / (pi * (nbhd_radius ^ 2)),
                cano_density = sum(abh[which(species == "CANO9")]) / (pi * (nbhd_radius ^ 2)),
                pico_density = sum(abh[which(species == "PICO")]) / (pi * (nbhd_radius ^ 2)),
                psme_density = sum(abh[which(species == "PSME")]) / (pi * (nbhd_radius ^ 2)))
                
    # Append output data to vectors
    x_coord <- c(x_coord, focal_stand_summary$x_coord)
    y_coord <- c(y_coord, focal_stand_summary$y_coord)
    num_trees <- c(num_trees, focal_stand_summary$num_trees)
    all_density <- c(all_density, focal_stand_summary$all_density)
    tshe_density <- c(tshe_density, focal_stand_summary$tshe_density)
    abam_density <- c(abam_density, focal_stand_summary$abam_density)
    thpl_density <- c(thpl_density, focal_stand_summary$thpl_density)
    tsme_density <- c(tsme_density, focal_stand_summary$tsme_density)
    cano_density <- c(cano_density, focal_stand_summary$cano_density)
    pico_density <- c(pico_density, focal_stand_summary$pico_density)
    psme_density <- c(psme_density, focal_stand_summary$psme_density)
    
  }
  
  # Combine vectors into output table
  output <- data.frame(x_coord, y_coord, num_trees, all_density,
                tshe_density, abam_density, thpl_density, tsme_density,
                cano_density, pico_density, psme_density)
  
  # Return output
  output
}

#======================================================================
# Calculate neighborhood density for all trees in a multi-stand dataset
#======================================================================

nbhd_density_all <- function(all_stands){
  
  # Identify unique stands in dataset
  stand_ids <- unique(all_stands$stand_id)
  
  # Loop through stands
  for(stand in 1:length(stand_ids)){
    
    # Isolate focal stand
    single_stand <- all_stands[all_stands$stand_id == stand_ids[stand], ]
    
    # If first stand, make the output table
    if(stand == 1){
      
      output <- nbhd_density(mapping_data = single_stand,
                             stand = stand_ids[stand],
                             x = single_stand[, "x_coord"], 
                             y = single_stand[, "y_coord"],
                             nbhd_radius = 10)
    } else { 
      
      # If not first stand, append results to those from earlier stands
      new_stand <- nbhd_density(mapping_data = single_stand,
                                stand = stand_ids[stand],
                                x = single_stand[, "x_coord"], 
                                y = single_stand[, "y_coord"],
                                nbhd_radius = 10)
      
      output <- rbind(output, new_stand)
    }
  }
  
  # Return output
  output
  
}



#===============================================
# Vectorized calculation of neighborhood density
#===============================================

density_summary <- function(mapping, stand, radius){
  
  # Extract relevant coordinates
  coords <- mapping %>%
    filter(stand_id == stand) %>%
    mutate(abh = (pi * ((dbh / 2) ^ 2)) / (pi * (radius ^ 2))) %>%
    select(tree_id, species, abh, x_coord, y_coord)

  # Create matrix of distances between each tree pair
  dist_mat <- as.matrix(dist(coords[, 4:5], diag = T, upper = T))
  
  # Create matrix where each column contains dbh measurements of all trees
  nbhds <- matrix(coords$abh, nrow = nrow(coords), ncol = nrow(coords))
  
  # Make each column contain only the dbh values of the trees in the
  # neighborhood represented by that column by changing unwanted dbh values
  # to NAs
  nbhds[which(dist_mat > radius)] <- NA
  
  # Add species information
  all <-  cbind(coords$species, as.data.frame(nbhds))
  colnames(all)[1] <- "species"
  
  # Create subsets
  tshe <- all %>%
    filter(species == "TSHE")
  abam <- all %>%
    filter(species == "ABAM")
  thpl <- all %>%
    filter(species == "THPL")
  tsme <- all %>%
    filter(species == "TSME")
  cano <- all %>%
    filter(species == "CANO9")
  pico <- all %>%
    filter(species == "PICO")
  psme <- all %>%
    filter(species == "PSME")
    
  # Calculate densities
  all_density <- apply(all[, 2:ncol(all)], 2, sum, na.rm = T)
  tshe_density <- apply(tshe[, 2:ncol(all)], 2, sum, na.rm = T)
  abam_density <- apply(abam[, 2:ncol(all)], 2, sum, na.rm = T)
  thpl_density <- apply(thpl[, 2:ncol(all)], 2, sum, na.rm = T)
  tsme_density <- apply(tsme[, 2:ncol(all)], 2, sum, na.rm = T)
  cano_density <- apply(cano[, 2:ncol(all)], 2, sum, na.rm = T)
  pico_density <- apply(pico[, 2:ncol(all)], 2, sum, na.rm = T)
  psme_density <- apply(psme[, 2:ncol(all)], 2, sum, na.rm = T)
  
  # Combine different density measures into data frame
  output <- data.frame(coords$tree_id, all_density, tshe_density,
                       abam_density, thpl_density, tsme_density,
                       cano_density, pico_density, psme_density)
  colnames(output)[1] <- "tree_id"
  
  # Return output
  output
}

#======================================================================
# Calculate neighborhood density for all trees in a multi-stand dataset
#======================================================================

density_all_stands <- function(all_stands){
  
  # Identify unique stands in dataset
  stand_ids <- unique(all_stands$stand_id)
  
  # Loop through stands
  for(stand_num in 1:length(stand_ids)){
    
    # If first stand, make the output table
    if(stand_num == 1){
      
      output <- density_summary(all_stands,
                                stand = stand_ids[stand_num],
                                radius = 10)
    } else { 
      
      # If not first stand, append results to those from earlier stands
      new_stand <- density_summary(all_stands,
                                   stand = stand_ids[stand_num],
                                   radius = 10)
      
      output <- rbind(output, new_stand)
    }
  }
  
  # Return output
  output
}