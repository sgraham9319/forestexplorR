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
mapping <- mapping[-which(mapping$tree_id %in% small_trees), ]

#=================
# Summarizing data 
#=================

# Summarize growth for all trees
growth <- growth_summary(growth_data)

# Calculate neighborhood density for all trees
densities <- nbhd_density_all(mapping)

# Attach neighborhood information to growth
colnames(growth)[1] <- "tree_id"
growth <- inner_join(growth, densities)

# Note - tree id AG05001500036 appears twice with same information apart from 
# mapping (presumably new mapping was added without removing old mapping). Need
# to deal with this in the data cleaning phase
View(growth[growth$tree_id == "AG05001500036", ])
