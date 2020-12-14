
# Load required packages
library(dplyr)

# Create data frame to store parameter values of best models
col_names <- c("focal_sps", "training_set",
               "sps_compABAM", "sps_compABLA", "sps_compCANO", "sps_compOTHR",
               "sps_compPICO", "sps_compPIMO", "sps_compPSME", "sps_compTABR",
               "sps_compTHPL", "sps_compTSHE", "sps_compTSME",
               "prox", "size_comp_dbh",
               "ABAM_density", "ABLA_density", "CANO_density", "OTHR_density",
               "PICO_density", "PIMO_density", "PSME_density", "TABR_density",
               "THPL_density", "TSHE_density", "TSME_density", "all_density",
               "intra")
best_mods <- as.data.frame(matrix(NA, nrow = 1, ncol = length(col_names)))
names(best_mods) <- col_names

# Create vectors of focal species and training sets
focal_sps <- c("PSME", "TSHE", "TSME", "THPL", "ABAM", "CANO")
train_set <- 1:4

# Open loop for focal species
for(i in focal_sps){
  
  # Open loop for training sets
  for(j in train_set){
    
    # Load model output
    mods <- read.csv(paste("Data/", i, j, "_no_clim.csv", sep = ""))
    
    # Extract best model
    new_mod <- mods[which.min(mods$MSE), ]
    
    # Add training set and focal species info
    new_mod$training_set <- j; new_mod$focal_sps <- i
    
    # Add to best models data frame
    best_mods <- bind_rows(best_mods, new_mod)
    
  }
}

# Remove empty initializing row and unwanted columns from best models data frame
best_mods <- best_mods[-1, 1:length(col_names)]

# Write best models table to file
write.csv(best_mods, file = "Data/RR_best_mod_coef.csv", row.names = F)
