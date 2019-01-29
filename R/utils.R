
#======================================================
# Calculate annual growth over entire measurment period
#======================================================

# Function requires a data frame where separate measurements of the same tree
# appear in different rows. The data frame needs to have a column called "treeid"
# that contains unique tree ID values, a column named "year" containing the year
# of the measurement, and a column named "dbh" containing the dbh measurements.

overall_annual_growth <- function(data){
  data %>% 
    group_by(treeid) %>% 
    arrange(year) %>%
    filter(
      year == max(year) | 
        year == min(year)
    ) %>%
    summarize(
      annual_growth = (dbh[2] - dbh[1]) / (year[2] - year[1]),
      begin_size = dbh[1]
    ) %>% 
    mutate(sqrt_annual_growth = sqrt(annual_growth),
           size_adj_sqrt_growth = sqrt(annual_growth / begin_size))
}


#==============================================
# Calculate annual growth over a defined period
#==============================================

# Function requires a data frame where separate measurements of the same tree
# appear in different rows. The data frame needs to have a column called "treeid"
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
    group_by(treeid) %>%
    filter(begin_dif == min(begin_dif))
  
  # Subset measurements later than end
  after <- data %>%
    filter(year >= end) %>%
    group_by(treeid) %>%
    filter(end_dif == min(end_dif))
  
  # Combine subsets
  before_after <- rbind(before, after)
  
  # Calculate annual growth
  period_growth = before_after %>% 
    group_by(treeid) %>% 
    arrange(year) %>%
    summarize(
      annual_growth = (dbh[2] - dbh[1]) / (year[2] - year[1])
    ) %>% 
    filter(!is.na(annual_growth)) %>%
    mutate(sqrt_annual_growth = sqrt(annual_growth))
  
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
    filter(StandID == stand)
  
  # Create output table
  x_coord <- NA
  y_coord <- NA
  num_trees <- NA
  total_abh <- NA
  output <- data.frame(x_coord, y_coord, num_trees, total_abh)
  
  # Begin looping through input coordinates
  for(coord_num in 1:length(x)){
    
    focal_stand_summary <- focal_stand %>%
      
      # Calculate distance of each tree from input coordinates
      mutate(dist = sqrt((abs(Xcoord - x[coord_num]) ^ 2) +
                           (abs(Ycoord - y[coord_num]) ^ 2))) %>%
      
      # Subset to trees within input radius of input coordinates
      filter(dist <= nbhd_radius) %>%
      
      # Calculate cross-sectional area at breast height (abh) for each tree
      mutate(abh = pi * ((Dbh_cm / 2) ^ 2)) %>%
      
      # Extract number of trees and total abh within neighborhood
      summarize(x_coord = x[coord_num],
                y_coord = y[coord_num],
                num_trees = n(),
                total_abh = sum(abh))
    
    # Append to output table
    output <- rbind(output, focal_stand_summary)
  }
  
  # Return output with first row (NAs) removed
  output[-1, ]
}