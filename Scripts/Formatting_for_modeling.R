# Load package
devtools::load_all()

#==============
# Cleaning data
#==============

# Load growth data
growth_data <- read.csv("../Data/Tree_growth_2017.csv", stringsAsFactors = F)

# Remove rows without dbh measurement
growth_data <- growth_data[!is.na(growth_data$dbh), ]

# Identify trees that were never measured as 15 cm dbh or larger
small_trees_data <- growth_data %>% group_by(treeid) %>% 
  summarize(max_size = max(dbh)) %>%
  filter(max_size < 15)
small_trees <- unique(small_trees_data$treeid)

# Exclude small trees
growth_data <- growth_data[-which(growth_data$treeid %in% small_trees), ]

# Load mapping data for individual trees
mapping <- read.csv("../Data/Mapping_2013.csv", stringsAsFactors = F)

# Exclude small trees from mapping data
mapping <- mapping[-which(mapping$TreeID %in% small_trees), ]

#=================
# Summarizing data 
#=================

# Calculate neighborhood density for all trees
densities <- nbhd_density_all(mapping)

# Summarize growth for all trees
growth <- growth_summary(growth_data)
