

# Checking growth model


# Load training data from TreeGrowth
train <- read.csv("../../TreeGrowth/Data/Output_data/training1.csv",
                  stringsAsFactors = F)

# Remove unneeded columns
train <- train %>%
  select(tree_id, species, sps_comp, dbh_comp, prox, all_density, ABAM_density,
         ABPR_density, CANO_density, PSME_density, TABR_density, THPL_density,
         TSHE_density, ABLA_density, TSME_density, ALSI_density, ALVI_density,
         PIMO_density, ALRU_density, PICO_density, PIEN_density, POBA_density,
         ABGR_density, PISI_density, pet_mm, size_corr_growth)

# Run growth model function
results <- growth_model(train, outcome_var = "size_corr_growth",
                        focal_sps = "PSME")

iterations <- 10
