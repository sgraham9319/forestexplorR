#####################################################
# Functions for rearranging data for likelihood model
#####################################################

# Author: Stuart Graham
# Created: 4/3/2020
# Last edited: 4/3/2020

# Functions list:
#     neighborhoods - makes a separate row for every focal-competitor pair
#                     in a given stand with given neighborhood radius
#     neighborhoods_all - implements neighborhoods function over all stands in
#                         a multiple-stand dataset with a provided neighborhood
#                         radius


#==============
# neighborhoods
#==============

# Note that the TreeNeighborhood package already contains a similar function
# but it returns abh (proportion of neighborhood area occupied by the tree)
# rather than dbh. We want to use dbh for the likelihood model to be consistent
# with Fortunel et al. 2016

# Function for creating neighborhoods table with dbh values
neighborhoods <- function(mapping, stand, radius) {
  
  # Convert species column to factor
  mapping$species <- as.factor(mapping$species)
  
  # Subset to focal stand and remove unneeded columns
  one_stand <- mapping %>%
    filter(stand_id == stand) %>%
    select(tree_id, stand_id, species, dbh, x_coord, y_coord, size_cat)
  
  # Create distance matrix - element i,j is the distance (in meters) between trees in
  # rows i and j in one_stand
  dist_mat <- as.matrix(dist(one_stand[, c("x_coord", "y_coord")]))
  # Create species matrix - each column is a replica of one_stand$species
  sps_mat <- matrix(one_stand$species, nrow = nrow(one_stand), ncol = nrow(one_stand))
  # Create size matrix - each column is a replica of one_stand$dbh
  size_mat <- matrix(one_stand$dbh, nrow = nrow(one_stand), ncol = nrow(one_stand))
  # Create size category matrix - each column is a replica of one_stand$size_cat
  cat_mat <- matrix(one_stand$size_cat, nrow = nrow(one_stand), ncol = nrow(one_stand))
  
  # Find distances greater than radius or equal to zero in distance matrix - these
  # elements represent pairs of trees that are not in each other's neighborhood or
  # are the same tree paired with itself
  non_comp <- which(dist_mat > radius | dist_mat == 0)
  
  # Convert values of these non-competing pairs to NA in all matrices
  dist_mat[non_comp] <- NA
  sps_mat[non_comp] <- NA
  size_mat[non_comp] <- NA
  cat_mat[non_comp] <- NA
  
  # Create dist, species, size and size category vectors not including NAs
  prox <- dist_mat[!is.na(dist_mat)]
  sps_comp <- sps_mat[!is.na(sps_mat)]
  size_comp <- size_mat[!is.na(size_mat)]
  size_cat_comp <- cat_mat[!is.na(cat_mat)]
  comp_dat <- data.frame(prox, sps_comp, size_comp, size_cat_comp)
  
  # Get vector of how many rows (competitors) each focal tree should have
  repeats <- unname(apply(dist_mat, 2, non_na_len))
  
  # Repeat rows of one_stand appropriate number of times
  all_rows <- one_stand[rep(1:nrow(one_stand), repeats), ]
  
  # Bind competitor data to one_stand
  output <- cbind(all_rows, comp_dat)
}


#==================
# neighborhoods_all
#==================

# This function takes as input a table of mapping data for trees distributed
# across multiple stands, implements the above neighborhoods function for 
# each stand, and then binds results together in a single output table

# Applying neighborhoods function to all stands
neighborhoods_all <- function(all_stands, radius){
  
  # Identify unique stands in dataset
  stand_ids <- unique(all_stands$stand_id)
  
  # Loop through stands
  for(stand_num in 1:length(stand_ids)){
    
    # If first stand, make the output table
    if(stand_num == 1){
      
      output <- neighborhoods(all_stands,
                              stand = stand_ids[stand_num],
                              radius = radius)
    } else { 
      
      # If not first stand, append results to those from earlier stands
      new_stand <- neighborhoods(all_stands,
                                 stand = stand_ids[stand_num],
                                 radius = radius)
      
      output <- rbind(output, new_stand)
    }
  }
  
  # Return output
  output
}
