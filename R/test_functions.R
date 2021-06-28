
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

#===============
# growth_summary
#===============

# Test input
growth_summ_test <- function(data){
  as.data.frame(growth_summary(data))
}

# Expected output
growth_summ_expt <- function(){
  tree_id <- c("tree1", "tree2", "tree3")
  stand_id <- rep("A", times = 3)
  species <- rep("TSHE", times = 3)
  first_record <- rep(2000, times = 3)
  last_record <- rep(2010, times = 3)
  begin_size <- rep(1, times = 3)
  final_size <- c(3, 11, 21)
  mean_size <- c(2, 6, 11)
  midpoint_size <- c(2, 6, 11)
  annual_growth <- c(0.2, 1, 2)
  size_corr_growth <- sqrt(c(0.2, 1, 2))
  dat <- data.frame(tree_id, stand_id, species, first_record, last_record,
                    begin_size, final_size, mean_size, midpoint_size,
                    annual_growth, size_corr_growth)
  dat
}

#================
# detailed_growth
#================

# Test input
det_growth_test <- function(data){
  as.data.frame(detailed_growth(data))
}

# Expected output
det_growth_expt <- function(){
  fake_growth() %>%
    mutate(annual_growth = c(0.2, 0.2, NA, 1, 1, NA, 2, 2, NA))
}

#=============================
# defined_period_annual_growth
#=============================

# Test input
def_growth_test <- function(data){
  as.data.frame(defined_period_annual_growth(data, 2005, 2010))
}

# Expected output
def_growth_expt <- function(){
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
  tree_id <- c("A1", "A2", "A3", "A4", "A5", "A6", "A7", "A8", "A9")
  stand_id <- rep("A", times = 9)
  species <- rep(c("TSHE", "ABAM"), length.out = 9)
  dbh <- rep(2, times = 9)
  x_coord <- rep(c(45, 50, 55), times = 3)
  y_coord <- rep(c(45, 50, 55), each = 3)
  A <- data.frame(tree_id, stand_id, species, dbh, x_coord, y_coord)
  
  # Fake stand B
  tree_id <- c("B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8", "B9")
  stand_id <- rep("B", times = 9)
  species <- rep("TSHE", times = 9)
  dbh <- rep(2, times = 9)
  x_coord <- rep(c(45, 50, 55), times = 3)
  y_coord <- rep(c(45, 50, 55), each = 3)
  B <- data.frame(tree_id, stand_id, species, dbh, x_coord, y_coord)
  
  # Combine fake stands into single data frame
  dat <- rbind(A, B)
  dat
}

#=============
# density_calc
#=============

# Test input
density_calc_test <- function(){
  species <- c("TSHE", "ABAM", "PSME")
  loc_one <- c(3, NA, 1)
  loc_two <- c(NA, 2, 1)
  loc_three <- c(NA, NA, 1)
  data.frame(species, loc_one, loc_two, loc_three)
}

# Expected output
density_calc_expt <- function(){
  all_density <- c(4, 3, 1)
  ABAM_density <- c(0, 2, 0)
  PSME_density <- c(1, 1, 1)
  TSHE_density <- c(3, 0, 0)
  data.frame(all_density, ABAM_density, PSME_density, TSHE_density)
}

#================
# density_summary
#================

# Test input
density_summ_test <- function(data, stand_id){
  output <- density_summary(data, stand_id, radius = 10)
  output[] <- lapply(output, as.character)
  output
}

# Expected output
density_summ_expt <- function(){
  
  # Fake stand A
  tree_id <- c("A1", "A2", "A3", "A4", "A5", "A6", "A7", "A8", "A9")
  x_coord <- rep(c(45, 50, 55), times = 3)
  y_coord <- rep(c(45, 50, 55), each = 3)
  all_density <- c(5, 6, 5, 6, 8, 6, 5, 6, 5) * 0.01
  ABAM_density <- c(2, 3, 2, 3, 4, 3, 2, 3, 2) * 0.01
  TSHE_density <- c(3, 3, 3, 3, 4, 3, 3, 3, 3) * 0.01
  dat <- data.frame(tree_id, x_coord, y_coord, all_density, ABAM_density,
                    TSHE_density)
  dat[] <- lapply(dat, as.character)
  dat
}

#===================
# density_all_stands
#===================

# Test input
density_all_test <- function(data){
  output <- density_all_stands(data, radius = 10)
  output[] <- lapply(output, as.character)
  output
}

# Expected output
density_all_expt <- function(){
  
  # Create expected density_summary output for fake stand B
  tree_id <- c("B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8", "B9")
  all_density <- c(5, 6, 5, 6, 8, 6, 5, 6, 5) * 0.01
  x_coord <- rep(c(45, 50, 55), times = 3)
  y_coord <- rep(c(45, 50, 55), each = 3)
  ABAM_density <- rep(0, times = 9)
  TSHE_density <- c(5, 6, 5, 6, 8, 6, 5, 6, 5) * 0.01
  output <- data.frame(tree_id, all_density, x_coord, y_coord,
                    ABAM_density, TSHE_density)
  
  # Combine with expected density_summary output for fake stand A 
  output <- rbind(density_summ_expt(), output)
  output[] <- lapply(output, as.character)
  output
}

#=================
# density_specific
#=================

# Test input
density_spec_test <- function(data, stand_id){
  x_coord <- c(60, 40, 50)
  y_coord <- c(60, 50, 50)
  focal_coords <- data.frame(x_coord, y_coord)
  output <- density_specific(data, stand_id, 10, focal_coords)
  output
}

# Expected output
density_spec_expt <- function(){
  x_coord <- c(60, 40, 50)
  y_coord <- c(60, 50, 50)
  all_density <- c(1, 4, 9) * 0.01
  ABAM_density <- c(0, 1, 4) * 0.01
  TSHE_density <- c(1, 3, 5) * 0.01
  dat <- data.frame(x_coord, y_coord, all_density, ABAM_density, TSHE_density)
  dat
}

#=============
# graph_matrix
#=============

expt <- function(){
  
  # Fake stand A
  tree_id <- c("A1", "A2", "A3", "A4", "A5", "A6", "A7", "A8", "A9")
  all_density <- c(5, 6, 5, 6, 8, 6, 5, 6, 5) * 0.01
  ABAM_density <- c(2, 3, 2, 3, 4, 3, 2, 3, 2) * 0.01
  TSHE_density <- c(3, 3, 3, 3, 4, 3, 3, 3, 3) * 0.01
  dat <- data.frame(tree_id, all_density, ABAM_density, TSHE_density)
  output <- dat[c(rep(1, times = 5), rep(2, times = 6), rep(3, times = 5),
                  rep(4, times = 6), rep(5, times = 8), rep(6, times = 6),
                  rep(7, times = 5), rep(8, times = 6), rep(9, times = 5)), ]
  rownames(output) <- 1:nrow(output)
  output
}