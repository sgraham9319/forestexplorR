# Load package
devtools::load_all()

#==============
# Cleaning data
#==============

# Load growth data
growth_data <- read.csv("../Data/Tree_growth_2017.csv", stringsAsFactors = F)

# Remove test data (2017 measurements and stands TO04, AE10, and AV02)
growth_data <- growth_data %>%
  filter(year != 2017,
         stand_id != "TO04",
         stand_id != "AE10",
         stand_id != "AV02")

# Remove rows without dbh measurement
growth_data <- growth_data[!is.na(growth_data$dbh), ]

# Identify trees that were never measured as 15 cm dbh or larger
small_trees_data <- growth_data %>% group_by(tree_id) %>% 
  summarize(max_size = max(dbh)) %>%
  filter(max_size < 15)
small_trees <- unique(small_trees_data$tree_id)

# Exclude small trees
growth_data <- growth_data[-which(growth_data$tree_id %in% small_trees), ]

# Load mapping data for individual trees
mapping <- read.csv("../Data/Mapping_2017.csv", stringsAsFactors = F)

# Exclude small trees from mapping data
mapping <- mapping[-which(mapping$tree_id %in% small_trees), ]

#=================
# Summarizing data 
#=================

# Summarize growth for all trees
growth <- growth_summary(growth_data)

# Calculate neighborhood density for all trees
densities <- nbhd_density_all(mapping)

# Add tree ids to densities
densities$tree_id <- mapping$tree_id

# Attach neighborhood information to growth
growth <- left_join(growth, densities)


# SHOULD WE EXCLUDE UNCOMMON TREE SPECIES?

