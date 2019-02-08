
#==========================================
# Function to create a fake mapping dataset
#==========================================

fake_data <- function(){
  
  # Fake stand A
  stand_id <- rep("A", times = 9)
  species <- rep(c("TSHE", "ABAM"), length.out = 9)
  dbh <- rep(2, times = 9)
  x_coord <- rep(c(45, 50, 55), times = 3)
  y_coord <- rep(c(45, 50, 55), each = 3)
  A <- data.frame(stand_id, species, dbh, x_coord, y_coord)
  
  # Combine fake stands into single data frame
  dat <- A
  dat
}

#==================================================================
# Function to create expected nbhd_density output for fake datasets
#==================================================================

density_expected <- function(){
  
  # Fake stand A
  x_coord <- rep(c(45, 50, 55), times = 3)
  y_coord <- rep(c(45, 50, 55), each = 3)
  num_trees <- c(6, 7, 6, 7, 9, 7, 6, 7, 6)
  all_density <- c(6, 7, 6, 7, 9, 7, 6, 7, 6) * 0.01
  tshe_density <- c(4, 3, 4, 3, 5, 3, 4, 3, 4) * 0.01
  abam_density <- c(2, 4, 2, 4, 4, 4, 2, 4, 2) * 0.01
  thpl_density <- rep(0, times = 9)
  tsme_density <- rep(0, times = 9)
  cano_density <- rep(0, times = 9)
  pico_density <- rep(0, times = 9)
  psme_density <- rep(0, times = 9)
  dat <- data.frame(x_coord, y_coord, num_trees, all_density,
                    tshe_density, abam_density, thpl_density, tsme_density,
                    cano_density, pico_density, psme_density)
  dat
}

#============================================
# Function to apply nbhd_density to fake data
#============================================

density_test <- function(data, stand_id){
  nbhd_density(data, stand = stand_id, x = data[, "x_coord"],
               y = data[, "y_coord"], nbhd_radius = 10)
}

