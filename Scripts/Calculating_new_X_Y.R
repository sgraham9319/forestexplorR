
# Load growth data
growth <- read.csv("../Data/Tree_growth_2017.csv", stringsAsFactors = F)

table(growth$new_mapping)

# Load new mapping data
new_map <- read.csv("../Data/new_mapping_2017.csv", stringsAsFactors = F)


# Load old mapping data
old_map <- read.csv("../Data/Mapping_2013.csv", stringsAsFactors = F)

stand_ids <- unique(new_map$standid)
stand <- new_map[new_map$standid == stand_ids[1], ]

# Calculate degrees to radians conversion factor
cf <- pi / 180

