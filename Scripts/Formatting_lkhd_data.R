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
# Part 2. Excluding test data - if using original training set
#############################

# Remove original test data from tree dataset (2017 measurements and stands TO04, AE10, and AV02)
#tree <- tree %>%
#  filter(year != 2017,
#         stand_id != "TO04",
#         stand_id != "AE10",
#         stand_id != "AV02")

# Remove original test data from mapping dataset (stands TO04, AE10, and AV02)
#mapping <- mapping %>%
#  filter(stand_id != "TO04",
#         stand_id != "AE10",
#         stand_id != "AV02")

################################
# Part 3. Creating neighborhoods
################################

# Define neighborhood radius
nb_rad <- 20

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

####################################
# Part 6. Defining new training data - unless using original training set
####################################

# Training 1
train_sub1 <- full %>% filter(y_coord <= 40)
train_sub2 <- full %>% filter(x_coord >= 60 & y_coord > 40)
training <- rbind(train_sub1, train_sub2)

# Training 2
train_sub1 <- full %>% filter(x_coord <= 40)
train_sub2 <- full %>% filter(y_coord <= 40 & x_coord > 40)
training <- rbind(train_sub1, train_sub2)

# Training 3
train_sub1 <- full %>% filter(y_coord >= 60)
train_sub2 <- full %>% filter(x_coord <= 40 & y_coord < 60)
training <- rbind(train_sub1, train_sub2)

# Training 4
train_sub1 <- full %>% filter(x_coord >= 60)
train_sub2 <- full %>% filter(y_coord >= 60 & x_coord < 60)
training <- rbind(train_sub1, train_sub2)

################################
# Part 7. Exploring sample sizes
################################

# Find number of focals per species
#focal_summ <- full %>% group_by(tree_id) %>% summarize(species = species[1])
#table(focal_summ$species)
focal_summ <- training %>% group_by(tree_id) %>% summarize(species = species[1])
table(focal_summ$species)

# Find number of competitors per species for each focal species
comp_summ <- function(data, focal_sps) {
  single_sps <- data[data$species == focal_sps, ]
  output <- table(single_sps$sps_comp)
  output
}
comp_summ(full, "TSME")

####################################
# Part 8. Writing out formatted data
####################################

# Write data for likelihood models
#write.csv(full, paste("Data/lkhd_", nb_rad, "m.csv", sep = ""), row.names = F)
write.csv(training, paste("Data/new_train4_", nb_rad, "m.csv", sep = ""), row.names = F)
