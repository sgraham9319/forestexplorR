###########################################################
# Neighborhood analysis of MORA data using likelihood model
###########################################################

# Author: Stuart Graham
# Created: 4/3/2020
# Last edited: 6/9/2020

# Load TreeNeighborhood package
devtools::load_all()

######################
# Part 1. Loading data
######################

# Load mapping data
mapping <- read.csv("Data/Cleaned_mapping_2017.csv", stringsAsFactors = F)

# Load tree measurement data
tree <- read.csv("Data/Cleaned_tree_growth_2017.csv", stringsAsFactors = F)

# Load abiotic data
env <- read.csv("Data/stand_abiotic_data.csv", stringsAsFactors = F)

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

# Define neighborhood radius
nb_rad <- 10

# Obtain all neighborhood data
neighbors <- neighborhoods_all(mapping, nb_rad)

# Remove focals whose neighborhood overlaps stand boundary
neighbors <- neighbors %>%
  filter(x_coord >= nb_rad & x_coord <= 100 - nb_rad & y_coord >= nb_rad & y_coord <= 100 - nb_rad)

# Remove small competitors
neighbors <- neighbors %>%
  filter(size_cat_comp == "regular")

###################################
# Part 4. Calculating annual growth
###################################

# Calculate annual growth for all trees
growth <- growth_summary(tree)

# Produces some NA values. Remove those caused by a tree being measured only once,
# which would result in an annual growth measurement of NA
growth <- growth[growth$first_record != growth$last_record, ]

# Check for any remaining NA annual growth values
sum(is.na(growth$annual_growth))

# Investigate NaN values of size_corrected_growth
no_grow <- growth[is.nan(growth$size_corr_growth), ]
max(no_grow$annual_growth)

# NaN values caused by negative annual growth values. These must be caused by
# measurement inaccuracies because trees can't shrink. Converting them to 0s
# would change the distribution of the data in a biased way. Better to treat
# them as missing data and just remove them
growth <- growth[growth$annual_growth >= 0, ]

# Check data one more time
sum(is.nan(growth$size_corr_growth))

# Calculate radial growth to be consistent with Fortunel and Canham
growth$radial_growth <- growth$annual_growth / 2

##########################################################
# Part 5. Combining growth, neighborhood, and abiotic data
##########################################################

# Extract required columns from growth data
growth_cols <- growth[, c("tree_id", "midpoint_size", "radial_growth")]

# Join growth and neighborhood data. Use inner join because there will be
# no growth data for focals measured only once or with negative growth
full <- inner_join(neighbors, growth_cols, by = "tree_id")

# Add abiotic data
full <- left_join(full, env)

################################
# Part 6. Exploring sample sizes
################################

# Find number of focals per species
focal_summ <- full %>% group_by(tree_id) %>% summarize(species = species[1])
table(focal_summ$species)

# Find number of competitors per species for each focal species
comp_summ <- function(data, focal_sps) {
  single_sps <- data[data$species == focal_sps, ]
  output <- table(single_sps$sps_comp)
  output
}
comp_summ(full, "TSME")

####################################
# Part 7. Writing out formatted data
####################################

# Write data for likelihood models
write.csv(full, "Data/lkhd_data.csv", row.names = F)
