
#===========================
# Create fake growth dataset
#===========================

fake_growth <- function(){
  tree_id <- rep(c("tree1", "tree2", "tree3"), each = 3)
  stand_id <- rep("A", times = 9)
  species <- rep("TSHE", times = 9)
  year <- rep(c(2000, 2005, 2010), times = 3)
  dbh <- c(1, 2, 3, 1, 6, 11, 1, 11, 21)
  dat <- data.frame(tree_id, stand_id, species, year, dbh)
  dat
}

#=======================================================
# Create expected growth_summary output for fake dataset
#=======================================================

growth_summ_exp <- function(){
  tree_id <- c("tree1", "tree2", "tree3")
  stand_id <- rep("A", times = 3)
  species <- rep("TSHE", times = 3)
  first_record <- rep(2000, times = 3)
  last_record <- rep(2010, times = 3)
  begin_size <- rep(1, times = 3)
  mean_size <- c(2, 6, 11)
  midpoint_size <- c(2, 6, 11)
  annual_growth <- c(0.2, 1, 2)
  size_corr_growth <- sqrt(c(0.2, 1, 2))
  dat <- data.frame(tree_id, stand_id, species, first_record, last_record,
                    begin_size, mean_size, midpoint_size, annual_growth,
                    size_corr_growth)
  dat
}

#========================================================
# Create expected detailed_growth output for fake dataset
#========================================================

det_growth_exp <- function(){
  fake_growth() %>%
    mutate(annual_growth = c(0.2, 0.2, NA, 1, 1, NA, 2, 2, NA))
}

#=====================================================================
# Create expected defined_period_annual_growth output for fake dataset
#=====================================================================

def_growth_exp <- function(){
  tree_id <- c("tree1", "tree2", "tree3")
  stand_id <- rep("A", times = 3)
  species <- rep("TSHE", times = 3)
  begin_size <- c(2, 6, 11)
  annual_growth <- c(0.2, 1, 2)
  size_corr_growth <- sqrt(c(0.2, 1, 2) / c(2, 6, 11))
  dat <- data.frame(tree_id, stand_id, species, begin_size, annual_growth,
                    size_corr_growth)
  dat
}

#============================
# Create fake mapping dataset
#============================

fake_map <- function(){
  
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

#=====================================================
# Create expected nbhd_density output for fake dataset
#=====================================================

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

#================================
# Apply nbhd_density to fake data
#================================

density_test <- function(data, stand_id){
  nbhd_density(data, stand = stand_id, x = data[, "x_coord"],
               y = data[, "y_coord"], nbhd_radius = 10)
}

