# Load package
devtools::load_all()

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

ab08 <- mapping[mapping$StandID == "AB08", ]
ab08_density <- nbhd_density(mapping_data = ab08, stand = "AB08", x = ab08$Xcoord, 
                             y = ab08$Ycoord, nbhd_radius = 10)

# Calculating density data
result <- nbhd_density(mapping_data = ab08mapping, stand = "AB08", x = mapping$Xcoord, 
                       y = Ycoord, nbhd_radius = 10)