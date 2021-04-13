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
      annual_growth = if_else(year[n()] != year[1],
                              (dbh[n()] - dbh[1]) / (year[n()] - year[1]),
                              NA)
    ) %>% 
    mutate(size_corr_growth = sqrt(annual_growth / begin_size))
}