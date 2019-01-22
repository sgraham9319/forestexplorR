
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
      annual_growth = (dbh[2] - dbh[1]) / (year[2] - year[1])
    ) %>% 
    mutate(sqrt_annual_growth = sqrt(annual_growth))
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