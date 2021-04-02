

if(getRversion() >= "2.15.1")  
  utils::globalVariables(c(".", "year", "tree_id", "begin_dif", "end_dif", 
                           "stand_id", "species", "dbh", "annual_growth", 
                           "dbh_diff", "year_diff", "begin_size", "dist", 
                           "abh", "y_azim", "st_convex_hull", "stands_locs_raw", 
                           "gather", "location_id", "deg_min_sec", "conv_unit", 
                           "corner_id", "coord", "dec_deg", "spread", "stand_locs_raw", 
                           "stand_utm", "st_as_sf", "mapview", "stand_polygons"))

#==================================
# Get variable length excluding NAs
#==================================

non_na_len <- function(x){length(na.omit(x))}

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

#================================================
# Calculate densities (m^2/hectare) from abh sums
#================================================

# This function takes abh values and converts them into densities with units
# of m^2/hectare. The formula looks surprisingly simple because it is a 
# simplified version of: abh/10000 * 10000/neighborhood_area. The abh values
# are in cm^2 so must be divided by 10000 to convert to m^2, then to get the 
# density for a hectare we divide a hectare (10000 m^2) by the size of our
# sample area (neighborhood). The 10000s cancel to give this simple formula

density_conv <- function(dens, radius){
  dens / circ_area(radius)
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

density_calc <- function(data, radius){
  
  # Create list of species
  sps_list <- sort(unique(data$species))
  
  # Create output data frame
  output <- as.data.frame(matrix(nrow = ncol(data) - 1,
                                 ncol = length(sps_list) + 1))
  colnames(output)[1] <- "all_density"
  for(i in 1:length(sps_list)){
    colnames(output)[i + 1] <- paste(sps_list[i], "density", sep = "_")
  }
  
  if(nrow(output) > 1){
    
    # Add density including all species to output
    output[, "all_density"] <- apply(data[, 2:ncol(data)], 2, sum, na.rm = T)
    
    # Loop through species adding species-specific densities to output
    for(sps in 1:length(sps_list)){
      one_sps <- data %>%
        filter(species == sps_list[sps])
      output[, sps + 1] <- apply(one_sps[, 2:ncol(data)], 2, sum, na.rm = T)
    }
  } else {
    
    # Add density including all species to output
    output[, "all_density"] <- sum(data[, 2], na.rm = T)
    
    # Loop through species adding species-specific densities to output
    for(sps in 1:length(sps_list)){
      one_sps <- data %>%
        filter(species == sps_list[sps])
      output[, sps + 1] <- sum(one_sps[, 2], na.rm = T)
    }
  }
  
  # Convert densities to m^2/hectare
  output <- output %>%
    mutate_all(density_conv, radius)
  
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
    mutate(abh = circ_area(dbh / 2)) %>%
    select(tree_id, species, abh, x_coord, y_coord)
  
  # Create matrix of distances between each tree pair
  dist_mat <- as.matrix(dist(coords[, 4:5], diag = T, upper = T))
  
  # Create matrix where each column contains abh measurements of all trees
  nbhds <- matrix(coords$abh, nrow = nrow(coords), ncol = nrow(coords))
  
  # Change abh values of trees not in each neighborhood to NA
  nbhds[which(dist_mat > radius | dist_mat == 0)] <- NA
  
  # Add species information
  all <-  cbind(coords$species, as.data.frame(nbhds))
  colnames(all)[1] <- "species"
  
  # Summarize density data
  output <- density_calc(all, radius)
  
  # Add tree_id and x/y coordinate columns
  output <- cbind(coords$tree_id, coords$x_coord, coords$y_coord, output)
  colnames(output)[1:3] <- c("tree_id", "x_coord", "y_coord")
  
  # Remove focals whose neighborhood overlaps with stand boundary
  output <- output %>%
    filter(x_coord >= radius & x_coord <= 100 - radius & y_coord >= radius & y_coord <= 100 - radius)
  
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
      
      output <- bind_rows(output, new_stand)
    }
  }
  
  # Change any NA to 0
  output[, 4:ncol(output)][is.na(output[, 4:ncol(output)])] <- 0
  
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
    mutate(abh = circ_area(dbh / 2)) %>%
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
  output <- density_calc(all, radius)
  
  # Add input coordinates to output
  output <- cbind(input_coords, output)
  
  # Return output
  output
}

#==========================================================================
# Generate matrix from mapping data with one row for every interacting pair
#==========================================================================

graph_matrix <- function(mapping, stand, radius){
  
  # Determine if stand exists
  if(stand %in% mapping$stand_id == F){
    print("Error: stand ID not found in mapping data")
  } else {
    
    # Convert species column to factor
    mapping$species <- as.factor(mapping$species)
    
    # Subset to focal stand
    one_stand <- mapping %>%
      filter(stand_id == stand) %>%
      mutate(abh = circ_area(dbh / 2)) %>%
      select(tree_id, stand_id, species, dbh, abh, size_cat, x_coord, y_coord)
    
    # Create distance, species, and size matrices
    dist_mat <- as.matrix(dist(one_stand[, c("x_coord", "y_coord")], diag = T, upper = T))
    sps_mat <- matrix(one_stand$species, nrow = nrow(one_stand), ncol = nrow(one_stand))
    size_mat_dens <- matrix(one_stand$abh, nrow = nrow(one_stand), ncol = nrow(one_stand))
    size_mat <- matrix(one_stand$dbh, nrow = nrow(one_stand), ncol = nrow(one_stand))
    cat_mat <- matrix(one_stand$size_cat, nrow = nrow(one_stand), ncol = nrow(one_stand))
    
    # Convert values of non-competitors to NA
    non_comp <- which(dist_mat > radius | dist_mat == 0)
    dist_mat[non_comp] <- NA
    sps_mat[non_comp] <- NA
    size_mat_dens[non_comp] <- NA
    size_mat[non_comp] <- NA
    cat_mat[non_comp] <- NA
    
    # Create dist, species and size vectors not including NAs
    prox <- dist_mat[!is.na(dist_mat)]
    sps_comp <- sps_mat[!is.na(sps_mat)]
    size_comp <- size_mat[!is.na(size_mat)]
    size_cat_comp <- cat_mat[!is.na(cat_mat)]
    comp_dat <- data.frame(prox, sps_comp, size_comp, size_cat_comp)
    
    # Get vector of how many rows each focal tree should have
    repeats <- unname(apply(dist_mat, 2, non_na_len))
    
    # Create data frame with repeated rows
    all_rows <- one_stand[rep(1:nrow(one_stand), repeats), ]
    
    # Combine focal and comp data
    all_dat <- cbind(all_rows, comp_dat)
    
    # Add species information
    nbhds <-  cbind(one_stand$species, as.data.frame(size_mat_dens))
    colnames(nbhds)[1] <- "species"
    
    # Summarize density data
    dens <- density_calc(nbhds, radius)
    
    # Add tree_id column
    dens <- cbind(one_stand$tree_id, dens)
    colnames(dens)[1] <- "tree_id"
    dens$tree_id <- as.character(dens$tree_id)
    
    # Add density information to all_dat
    output <- left_join(all_dat, dens)
    
    # Return output
    output
  }
}

#=================================
# Apply graph_matrix to all stands
#=================================

graph_mat_all <- function(all_stands, radius){
  
  # Identify unique stands in dataset
  stand_ids <- unique(all_stands$stand_id)
  
  # Loop through stands
  for(stand_num in 1:length(stand_ids)){
    
    # If first stand, make the output table
    if(stand_num == 1){
      
      output <- graph_matrix(all_stands,
                                stand = stand_ids[stand_num],
                                radius = radius)
    } else { 
      
      # If not first stand, append results to those from earlier stands
      new_stand <- graph_matrix(all_stands,
                                   stand = stand_ids[stand_num],
                                   radius = radius)
      
      output <- bind_rows(output, new_stand)
    }
  }
  
  # Replace NA densities with zero
  dens_cols <- grep("density", names(output))
  output[, dens_cols][is.na(output[, grep("density", names(output))])] <- 0
  
  # Return output
  output
}

#===============================
# Standardize values in a matrix
#===============================

z_trans <- function(x){(x - mean(x)) / sd(x)}

#=======================================
# Calculate coefficient of determination
#=======================================

coef_det <- function(x){
  1 - (sum((x$observations - x$predictions)^2) / 
         sum((x$observations - mean(x$observations))^2))
}

#=========================
# Calculate standard error
#=========================

st_err <- function(x){sd(x, na.rm = T) / sqrt(length(na.omit(x)))}

#============================================================================
# Calculate percentage of values in a vector greater or less than a threshold
#============================================================================

pct_dif <- function(x, direction, value){
  if(direction == "greater"){
    sum(x > value) * 100 / length(x)
  } else if(direction == "less"){
    sum(x < value) * 100 / length(x)
  } else {
    print("direction input not recognized")
  }
}


#============================================
# Fit regularized regression model many times
#============================================

# This function takes as input a design matrix already formatted for the 
# regularized regression model, a data frame of the tree growth and neighborhood
# data used to build the design matrix, and the desired number of runs for the
# model. It fits the model this many times and outputs a matrix containing the
# fitted parameter values, coefficient of determination for the training data, 
# and the cross-validated mean square error for each model run. The function 
# also keeps track of the best model (lowest MSE) and includes this model object
# in the list output

mod_iter <- function(des_mat, growth_dat, runs){
  
  # Create output table
  results <- matrix(NA, nrow = runs, ncol = ncol(des_mat) + 3)
  
  # Iterate model fitting process
  for(i in 1:runs){
    
    # Fit model
    mod <- cv.glmnet(des_mat, y = growth_dat$size_corr_growth, family = "gaussian")
    
    # Store model coefficients
    results[i, 1:(ncol(des_mat) + 1)] <- as.vector(as.matrix(coef(mod)))
    
    # Make data frame of growth predictions and observations
    int_pred <- predict(mod, newx = des_mat, s = "lambda.1se")
    int_obs <- growth_dat %>% select(tree_id, size_corr_growth)
    comparison <- cbind(int_obs, int_pred)
    colnames(comparison)[2:3] <- c("obs", "pred")
    comparison <- comparison %>%
      group_by(tree_id) %>%
      summarize(observations = obs[1],
                predictions = mean(pred))
    
    # Store coefficient of determination
    results[i, ncol(des_mat) + 2] <- coef_det(comparison)
    
    # Store mean square error of 1se model
    results[i, ncol(des_mat) + 3] <- mod$cvm[mod$lambda == mod$lambda.1se]
    
    # Remove added intercept column and name remaining columns
    colnames(results) <- c(rownames(coef(mod)), "R_square", "MSE")
    
    # Store/update best model
    if(i == 1){
      best_mod <- mod
    } else if(mod$cvm[mod$lambda == mod$lambda.1se] < 
              best_mod$cvm[best_mod$lambda == best_mod$lambda.1se]){
      best_mod <- mod
    }
    
  }
  
  # Convert results to data frame and output results table and best mod
  results <- as.data.frame(results)
  output <- list(results, best_mod)
  
}