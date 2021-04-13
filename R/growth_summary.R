#======================================================
# Calculate annual growth over entire measurment period
#======================================================

# Function requires a data frame where separate measurements of the same tree
# appear in different rows. The data frame needs to have a column called "treeid"
# that contains unique tree ID values, a column named "year" containing the year
# of the measurement, and a column named "dbh" containing the dbh measurements.

growth_summary <- function(data){
  output <- data %>% 
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
      final_size = dbh[n()]
    ) %>% 
    mutate(
      annual_growth = if_else(
        first_record == last_record,
        NA_real_,
        (final_size - begin_size) / (last_record - first_record)))
  
  # Calculate size corrected growth for trees where defined
  output$size_corr_growth <- NA
  valid_rows <- which(output$annual_growth >= 0)
  output$size_corr_growth[valid_rows] <- 
    sqrt(output$annual_growth[valid_rows] /
           output$begin_size[valid_rows])
  
  # Give warning if some trees showed negative growth
  neg_grow <- length(which(output$annual_growth < 0))
  if(neg_grow > 0){
    print(paste("Warning:", neg_grow, "trees exhibited negative annual growth"))
  }
  
  return(output)
  
}