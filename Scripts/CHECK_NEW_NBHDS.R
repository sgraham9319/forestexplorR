library(tidyr)

# Load old data files
map_old <- read.csv("../../TreeGrowth/Data/Raw_data/mapping_2017.csv", stringsAsFactors = F)
tree_old <- read.csv("../../TreeGrowth/Data/Raw_data/tree_growth_2017.csv", stringsAsFactors = F)
abio <- read.csv("../../TreeGrowth/Data/Raw_data/stand_abiotic_data.csv", stringsAsFactors = F)

# Create new neighborhoods file
neighbors <- neighborhoods2(map_old, radius = 15)
nbhd_summ <- neighborhood_summary(neighbors, radius = 15)
new_nbhds <- neighbors %>%
  left_join(nbhd_summ, by = "tree_id")
nb_rad <- 15
new_nbhds <- new_nbhds %>%
  filter(x_coord >= nb_rad & x_coord <= 100 - nb_rad &
           y_coord >= nb_rad & y_coord <= 100 - nb_rad)
new_nbhds <- new_nbhds %>%
  left_join(map_old %>% select(tree_id, size_cat), by = c("id_comp" = "tree_id"))
new_nbhds <- new_nbhds %>%
  filter(size_cat == "regular")
growth <- growth_summary(tree_old)
growth <- growth %>%
  filter(growth$first_record != growth$last_record &
           annual_growth >= 0)

complete_nbhds <- new_nbhds %>%
  # Inner join here because not all trees in neighbors have growth data
  inner_join(growth %>% select(tree_id, midpoint_size, annual_growth,
                               size_corr_growth),
             by = "tree_id") %>%
  left_join(abio, by = "stand_id")
