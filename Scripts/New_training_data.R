######################################################
# Exploring options for new training/test data split #
######################################################

# This script estimates sample size (number of focal trees) for each species in
# the training and test datasets for three potential methods of dividing the 
# data between training and test. The methods are:
# 1) No overlap - any trees that are focals or neighbors in training cannot
#                 be focals or neighbors in test
# 2) Some overlap - trees that are only neighbors in training can be neighbors
#                 but not focals in test
# 3) Strong overlap - trees that are only neighbors in training can be focals
#                 in test

# Load TreeNeighborhood package
devtools::load_all()

#=============
# Loading data
#=============

# Load mapping data
mapping_raw <- read.csv("Data/Cleaned_mapping_2017.csv", stringsAsFactors = F)

# Load tree measurement data
tree_raw <- read.csv("Data/Cleaned_tree_growth_2017.csv", stringsAsFactors = F)

#=================================
# Removing trees with missing data
#=================================

# Calculate annual growth for all trees
growth <- growth_summary(tree_raw)

# Remove trees that were only measured once (growth cannot be calculated)
growth <- growth[growth$first_record != growth$last_record, ]

# Remove trees with negative annual growth
growth <- growth[growth$annual_growth >= 0, ]

# Confirm that growth is defined for all remaining trees
sum(is.nan(growth$annual_growth))
sum(is.nan(growth$size_corr_growth))

# Extract required columns from growth data
growth_cols <- growth[, c("tree_id", "midpoint_size", "size_corr_growth")]

# Join growth data with mapping data, removing trees without mapping
full <- inner_join(growth_cols, mapping_raw, by = "tree_id")

# Remove focals within 20m of stand boundary
potential_focals <- full %>%
  filter(x_coord >= 20 & x_coord <= 80 & y_coord >= 20 & y_coord <= 80)

# View number of potential focals to be divided among test and training
table(potential_focals$species)

#=====================
# Method 1: No overlap
#=====================

# Create training data set
train_sub1 <- potential_focals %>% filter(y_coord <= 30) # Model 1
train_sub2 <- potential_focals %>% filter(x_coord >= 70 & y_coord > 30) # Model 1
#train_sub1 <- potential_focals %>% filter(x_coord <= 30) # Model 2
#train_sub2 <- potential_focals %>% filter(y_coord <= 30 & x_coord > 30) # Model 2
#train_sub1 <- potential_focals %>% filter(y_coord >= 70) # Model 3
#train_sub2 <- potential_focals %>% filter(x_coord <= 30 & y_coord < 70) # Model 3
#train_sub1 <- potential_focals %>% filter(y_coord >= 70) # Model 4
#train_sub2 <- potential_focals %>% filter(x_coord >= 70 & y_coord < 70) # Model 4
train1 <- rbind(train_sub1, train_sub2)

# Create test data set
test1 <- potential_focals %>% filter(x_coord <= 30 & y_coord >= 70) # Model 1
#test1 <- potential_focals %>% filter(x_coord >= 70 & y_coord >= 70) # Model 2
#test1 <- potential_focals %>% filter(x_coord >= 70 & y_coord <= 30) # Model 3
#test1 <- potential_focals %>% filter(x_coord <= 30 & y_coord <= 30) # Model 4

# View number of focals in test and training
table(train1$species)
table(test1$species)

#=======================
# Method 2: Some overlap
#=======================

# Create training data set
train_sub1 <- potential_focals %>% filter(y_coord <= 40) # Model 1
train_sub2 <- potential_focals %>% filter(x_coord >= 60 & y_coord > 40) # Model 1
#train_sub1 <- potential_focals %>% filter(x_coord <= 40) # Model 2
#train_sub2 <- potential_focals %>% filter(y_coord <= 40 & x_coord > 40) # Model 2
#train_sub1 <- potential_focals %>% filter(y_coord >= 60) # Model 3
#train_sub2 <- potential_focals %>% filter(x_coord <= 40 & y_coord < 60) # Model 3
#train_sub1 <- potential_focals %>% filter(y_coord >= 60) # Model 4
#train_sub2 <- potential_focals %>% filter(x_coord >= 60 & y_coord < 60) # Model 4
train2 <- rbind(train_sub1, train_sub2)

# Create test data set
test2 <- potential_focals %>% filter(x_coord <= 40 & y_coord >= 60) # Model 1
#test2 <- potential_focals %>% filter(x_coord >= 60 & y_coord >= 60) # Model 2
#test2 <- potential_focals %>% filter(x_coord >= 60 & y_coord <= 40) # Model 3
#test2 <- potential_focals %>% filter(x_coord <= 40 & y_coord <= 40) # Model 4

# View number of focals in test and training
table(train2$species)
table(test2$species)

#=========================
# Method 3: Strong overlap
#=========================

# Create training data set
train_sub1 <- potential_focals %>% filter(y_coord <= 50) # Model 1
train_sub2 <- potential_focals %>% filter(x_coord >= 50 & y_coord > 50) # Model 1
#train_sub1 <- potential_focals %>% filter(x_coord <= 50) # Model 2
#train_sub2 <- potential_focals %>% filter(y_coord <= 50 & x_coord > 50) # Model 2
#train_sub1 <- potential_focals %>% filter(y_coord >= 50) # Model 3
#train_sub2 <- potential_focals %>% filter(x_coord <= 50 & y_coord < 50) # Model 3
#train_sub1 <- potential_focals %>% filter(y_coord >= 50) # Model 4
#train_sub2 <- potential_focals %>% filter(x_coord >= 50 & y_coord < 50) # Model 4
train3 <- rbind(train_sub1, train_sub2)

# Create test data set
test3 <- potential_focals %>% filter(x_coord < 50 & y_coord > 50) # Model 1
#test3 <- potential_focals %>% filter(x_coord > 50 & y_coord > 50) # Model 2
#test3 <- potential_focals %>% filter(x_coord > 50 & y_coord < 50) # Model 3
#test3 <- potential_focals %>% filter(x_coord < 50 & y_coord < 50) # Model 4

# View number of focals in test and training
table(train3$species)
table(test3$species)
