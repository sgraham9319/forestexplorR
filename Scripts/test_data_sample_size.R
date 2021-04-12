


# Load TreeNeighborhood package
devtools::load_all()

# Load mapping data
mapping_raw <- read.csv("Data/Cleaned_mapping_2017.csv", stringsAsFactors = F)

# Load tree measurement data
tree_raw <- read.csv("Data/Cleaned_tree_growth_2017.csv", stringsAsFactors = F)

# Load abiotic data
env <- read.csv("Data/stand_abiotic_data.csv", stringsAsFactors = F)

#===============================
# Test stands (TO04, AE10, AV02)
#===============================

# Subset tree data to test stands
tree <- tree_raw %>%
  filter(stand_id == "TO04" | stand_id == "AE10" | stand_id == "AV02")

# Subset mapping data to test stands
mapping <- mapping_raw %>%
  filter(stand_id == "TO04" | stand_id == "AE10" | stand_id == "AV02")

# Calculate annual growth for all trees
growth <- growth_summary(tree)

# Count number of NAs in annual growth
sum(is.nan(growth$annual_growth)) # 126

# Remove those that result from a tree being measured only once
growth <- growth[growth$first_record != growth$last_record, ]
sum(is.nan(growth$annual_growth)) # 0

# Count number remaining NAs in size corrected growth
sum(is.nan(growth$size_corr_growth)) # 46

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

# Obtain all neighborhood data
neighbors <- neighborhoods_all(mapping, 20)

# Remove focals whose neighborhood overlaps stand boundary
neighbors <- neighbors %>%
  filter(x_coord >= 20 & x_coord <= 80 & y_coord >= 20 & y_coord <= 80)

# Remove small competitors
neighbors <- neighbors %>%
  filter(size_cat_comp == "regular")

# Extract required columns from growth data
growth_cols <- growth[, c("tree_id", "midpoint_size", "size_corr_growth")]

# Join growth and neighborhood data. Use inner join because there will be
# no growth data for focals measured only once or with negative growth
full <- inner_join(neighbors, growth_cols, by = "tree_id")

# Add abiotic data
full <- left_join(full, env, by = "stand_id")

# Subset to focal trees and count
focals <- full %>% group_by(tree_id) %>% summarize(species = species[1])
sum(table(focals$species))
# PSME 1, TSHE 154, TSME 1, THPL 9, ABAM 564, CANO 120

#======================
# 2012/2013-2017 growth
#======================

# Confirm tree data collected in 2017 from all stands
sub2017 <- tree_raw[tree_raw$year == 2017, ]
length(table(sub2017$stand_id)) == 15

# Confirm all stands measured in 2012 or 2013 but not both
sub2012 <- tree_raw[tree_raw$year == 2012, ]
sub2013 <- tree_raw[tree_raw$year == 2013, ]
table(sub2012$stand_id)
table(sub2013$stand_id)

# Subset tree data to focal stands and years 2012-2017
tree_new <- tree_raw %>%
  filter(year >= 2012 &
         stand_id != "TO04",
         stand_id != "AE10",
         stand_id != "AV02")

# Subset mapping data to focal stands
map_new <- mapping_raw %>%
  filter(stand_id != "TO04",
         stand_id != "AE10",
         stand_id != "AV02")

# Calculate annual growth for all trees
growth_new <- growth_summary(tree_new)

# Count number of NAs in annual growth
sum(is.nan(growth_new$annual_growth)) # 404 (7%)

# Remove those that result from a tree being measured only once
growth_new <- growth_new[growth_new$first_record != growth_new$last_record, ]
sum(is.nan(growth_new$annual_growth)) # 0

# Count number remaining NAs in size corrected growth
sum(is.nan(growth_new$size_corr_growth)) # 555

# Investigate NaN values of size_corrected_growth
no_grow <- growth_new[is.nan(growth_new$size_corr_growth), ]
max(no_grow$annual_growth)

# NaN values caused by negative annual growth values. These must be caused by
# measurement inaccuracies because trees can't shrink. Converting them to 0s
# would change the distribution of the data in a biased way. Better to treat
# them as missing data and just remove them
growth_new <- growth_new[growth_new$annual_growth >= 0, ]

# Check data one more time
sum(is.nan(growth_new$size_corr_growth))

# Obtain all neighborhood data
neighbors <- neighborhoods_all(map_new, 10)

# Remove focals whose neighborhood overlaps stand boundary
neighbors <- neighbors %>%
  filter(x_coord >= 20 & x_coord <= 80 & y_coord >= 20 & y_coord <= 80)

# Remove small competitors
neighbors <- neighbors %>%
  filter(size_cat_comp == "regular")

# Extract required columns from growth data
growth_cols <- growth_new[, c("tree_id", "midpoint_size", "size_corr_growth")]

# Join growth and neighborhood data. Use inner join because there will be
# no growth data for focals measured only once or with negative growth
full_new <- inner_join(neighbors, growth_cols, by = "tree_id")

# Add abiotic data
full_new <- left_join(full_new, env, by = "stand_id")

# Subset to focal trees and count
focals_new <- full_new %>% group_by(tree_id) %>% summarize(species = species[1])
table(focals_new$species)
# PSME 178, TSHE 637, TSME 74, THPL 81, ABAM 476, CANO 103

#========================
# k-fold cross-validation
#========================

#--------------
# Load raw data
#--------------

# Load mapping data
mapping_raw <- read.csv("Data/Cleaned_mapping_2017.csv", stringsAsFactors = F)

# Load tree measurement data
tree_raw <- read.csv("Data/Cleaned_tree_growth_2017.csv", stringsAsFactors = F)

# Load abiotic data
env <- read.csv("Data/stand_abiotic_data.csv", stringsAsFactors = F)



# Remove 2017 growth from tree data
tree <- tree_raw %>% filter(year != 2017)

# Calculate annual growth for all trees
growth <- growth_summary(tree)

# Count number of NAs in annual growth
sum(is.nan(growth$annual_growth)) # 459

# Remove those that result from a tree being measured only once
growth <- growth[growth$first_record != growth$last_record, ]
sum(is.nan(growth$annual_growth)) # 0

# Count number remaining NAs in size corrected growth
sum(is.nan(growth$size_corr_growth)) # 119

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

# Calculate radial growth for likelihood model
growth$radial_growth <- growth$annual_growth / 2

# Obtain all neighborhood data
neighbors <- neighborhoods_all(mapping_raw, 20)

# Remove small competitors
neighbors <- neighbors %>%
  filter(size_cat_comp == "regular")

# Extract required columns from growth data
growth_cols <- growth[, c("tree_id", "midpoint_size", "radial_growth", "size_corr_growth")]

# Join growth and neighborhood data. Use inner join because there will be
# no growth data for focals measured only once or with negative growth
full <- inner_join(neighbors, growth_cols, by = "tree_id")

# Add abiotic data
full <- left_join(full, env, by = "stand_id")

# Group by tree id to get one row per focal tree
focals <- full %>% group_by(tree_id) %>%
  summarize(species = species[1], x_coord = x_coord[1], y_coord = y_coord[1])

# Create subset for each focal species
PSME <- focals %>% filter(species == "PSME")
TSHE <- focals %>% filter(species == "TSHE")
TSME <- focals %>% filter(species == "TSME")
THPL <- focals %>% filter(species == "THPL")
ABAM <- focals %>% filter(species == "ABAM")
CANO <- focals %>% filter(species == "CANO")

# Create list of radii
rad <- c(18, 20, 15, 17, 20, 17)

# Remove focals whose neighborhood overlaps stand boundary for each species
PSME_sub <- PSME %>%
  filter(x_coord >= rad[1] & x_coord <= 100 - rad[1] & y_coord >= rad[1] & y_coord <= 100 - rad[1])
TSHE_sub <- TSHE %>%
  filter(x_coord >= rad[2] & x_coord <= 100 - rad[2] & y_coord >= rad[2] & y_coord <= 100 - rad[2])
TSME_sub <- TSME %>%
  filter(x_coord >= rad[3] & x_coord <= 100 - rad[3] & y_coord >= rad[3] & y_coord <= 100 - rad[3])
THPL_sub <- THPL %>%
  filter(x_coord >= rad[4] & x_coord <= 100 - rad[4] & y_coord >= rad[4] & y_coord <= 100 - rad[4])
ABAM_sub <- ABAM %>%
  filter(x_coord >= rad[5] & x_coord <= 100 - rad[5] & y_coord >= rad[5] & y_coord <= 100 - rad[5])
CANO_sub <- CANO %>%
  filter(x_coord >= rad[6] & x_coord <= 100 - rad[6] & y_coord >= rad[6] & y_coord <= 100 - rad[6])

# Check expected test and training sample sizes (based on rad list)
dat <- CANO_sub
quad1 <- dat %>% filter(x_coord <= 50, y_coord <= 50)
quad2 <- dat %>% filter(x_coord <= 50, y_coord > 50)
quad3 <- dat %>% filter(x_coord > 50, y_coord > 50)
quad4 <- dat %>% filter(x_coord > 50, y_coord <= 50)
nrow(quad1); nrow(quad2); nrow(quad3); nrow(quad4)
