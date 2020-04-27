#########################################
# Single-species directional graph models
#########################################

# Author: Stuart Graham
# Created: 4/27/2020
# Last edited: 4/27/2020

# Load TreeNeighborhood package
devtools::load_all()

######################
# Part 1. Loading data
######################

# Load mapping data
mapping <- read.csv("../TreeNeighborhood/Data/Cleaned_mapping_2017.csv", stringsAsFactors = F)

# Load tree measurement data
tree <- read.csv("../TreeNeighborhood/Data/Cleaned_tree_growth_2017.csv", stringsAsFactors = F)

#############################
# Part 2. Excluding test data
#############################

# Remove test data from tree dataset (2017 measurements and stands TO04, AE10, and AV02)
tree <- tree %>%
  filter(year != 2017,
         stand_id != "TO04",
         stand_id != "AE10",
         stand_id != "AV02")

# Remove test data from mapping dataset (stands TO04, AE10, and AV02)
mapping <- mapping %>%
  filter(stand_id != "TO04",
         stand_id != "AE10",
         stand_id != "AV02")

################################
# Part 3. Creating neighborhoods
################################

# Obtain all neighborhood data
neighbors <- neighborhoods_all(mapping, 10)

# Remove focals whose neighborhood overlaps stand boundary
neighbors <- neighbors %>%
  filter(x_coord >= 10 & x_coord <= 90 & y_coord >= 10 & y_coord <= 90)

###################################
# Part 4. Calculating annual growth
###################################

# Calculate annual growth for all trees
growth <- growth_summary(tree)

# Remove trees for which annual growth could not be calculated (see likelihood
# model script for details)
growth <- growth[!is.na(growth$size_corr_growth), ]

# Check data one more time
sum(is.nan(growth$size_corr_growth))

# Calculate radial growth to be consistent with Fortunel and Canham
growth$radial_growth <- growth$annual_growth / 2

################################################
# Part 5. Combining growth and neighborhood data
################################################

# Extract required columns from growth data
growth_cols <- growth[, c("tree_id", "midpoint_size", "radial_growth")]

# Join growth and neighborhood data. Use inner join because there will be
# no growth data for focals measured only once or with negative growth
full <- inner_join(neighbors, growth_cols, by = "tree_id")