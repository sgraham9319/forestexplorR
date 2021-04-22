library(dplyr)
devtools::load_all()

# Load neighborhood data
nbhds <- read.csv("../../TreeGrowth/Data/Output_data/training1.csv")

# Remove unneeded columns from nbhds data
nbhds <- nbhds %>%
  select(species, sps_comp, prox, dbh_comp, grep("density", names(nbhds)),
         pet_mm, size_corr_growth)



test <- growth_model(training = nbhds, outcome_var = "size_corr_growth",
                     focal_sps = "PSME")
test$R_squared

# INTRA AS OPTION FOR FUNCTION???