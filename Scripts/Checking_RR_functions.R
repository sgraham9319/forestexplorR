library(dplyr)
devtools::load_all()

# Specify focal species
focal_sps <- "PSME"

# Load training data
train <- read.csv("../../TreeGrowth/Data/Output_data/training1.csv",
                  stringsAsFactors = F)

# Load test data
test <- read.csv("../../TreeGrowth/Data/Output_data/test1.csv",
                 stringsAsFactors = F)

# Remove unneeded columns from training and test
train <- train %>%
  select(tree_id, species, sps_comp, dbh_comp, prox, all_density, ABAM_density,
         ABPR_density, CANO_density, PSME_density, TABR_density, THPL_density,
         TSHE_density, ABLA_density, TSME_density, ALSI_density, ALVI_density,
         PIMO_density, ALRU_density, PICO_density, PIEN_density, POBA_density,
         ABGR_density, PISI_density, pet_mm, size_corr_growth)
test <- test %>%
  select(tree_id, species, sps_comp, dbh_comp, prox, all_density, ABAM_density,
         ABPR_density, CANO_density, PSME_density, TABR_density, THPL_density,
         TSHE_density, ABLA_density, TSME_density, ALSI_density, ALVI_density,
         PIMO_density, ALRU_density, PICO_density, PIEN_density, POBA_density,
         ABGR_density, PISI_density, pet_mm, size_corr_growth)

# Subset to focal species and create additional explanatory variables
sing_sp <- train %>%
  filter(species == focal_sps) %>%
  mutate(
    intra = if_else(sps_comp == focal_sps, 1, 0),
    inter_dens = all_density - get(paste(focal_sps, "density", sep = "_")))
ss_test <- test %>%
  filter(species == focal_sps) %>%
  mutate(
    intra = if_else(sps_comp == focal_sps, 1, 0),
    inter_dens = all_density - get(paste(focal_sps, "density", sep = "_")))

# Load common competitors data and extract for focal species
comm_comp <- read.csv("../../TreeGrowth/Data/Output_data/common_comps.csv",
                      stringsAsFactors = F)
comm_comp <- comm_comp[, focal_sps]

# Get rare competitor density columns
rare_dens <- setdiff(names(sing_sp)[grep("density", names(sing_sp))],
                     c(paste(comm_comp, "density", sep = "_"), "all_density"))

# Convert rare competitors to OTHR and remove unneeded densities
sing_sp <- sing_sp %>%
  mutate(sps_comp = if_else(sps_comp %in% comm_comp, sps_comp, "OTHR"),
         OTHR_density = apply(sing_sp %>% select(rare_dens), 1, sum)) %>%
  select(-rare_dens)
ss_test <- ss_test %>%
  mutate(sps_comp = if_else(sps_comp %in% comm_comp, sps_comp, "OTHR"),
         OTHR_density = apply(ss_test %>% select(rare_dens), 1, sum)) %>%
  select(-rare_dens)




no_test <- growth_model(training = sing_sp, outcome_var = "size_corr_growth",
                     focal_sps = "PSME", iterations = 100)
with_test <- growth_model(training = sing_sp, outcome_var = "size_corr_growth",
                          focal_sps = "PSME", iterations = 100, test = ss_test)



with_test$R_squared
with_test$test_R_squared
View(with_test$mod_coef)
coef(with_test$mod)
