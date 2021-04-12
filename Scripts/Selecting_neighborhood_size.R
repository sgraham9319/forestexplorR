#############################
# Selecting neighborhood size
#############################

# This script models tree growth as a function of temperature, precipitation, and
# the density of trees in the competitive neighborhood. By using competitor densities
# calculated using a range of neighborhood sizes and comparing fit of the different
# models, this script makes suggestions on the most appropriate neighborhood size to
# use in the analysis

# Load package
devtools::load_all()

# Load required packages
library(glmnet)

#==================================
# Part 1. Loading and cleaning data
#==================================

# Load growth data
growth_data <- read.csv("Data/Cleaned_tree_growth_2017.csv", stringsAsFactors = F)

# Load mapping data for individual trees
mapping <- read.csv("Data/Cleaned_mapping_2017.csv", stringsAsFactors = F)

# Load environmental data
env_raw <- read.csv("Data/stand_abiotic_data.csv", stringsAsFactors = F)

# Remove test data (2017 measurements and stands TO04, AE10, and AV02)
#growth_data <- growth_data %>%
#  filter(year != 2017,
#         stand_id != "TO04",
#         stand_id != "AE10",
#         stand_id != "AV02")

# Remove test data from mapping dataset (stands TO04, AE10, and AV02)
#mapping <- mapping %>%
#  filter(stand_id != "TO04",
#         stand_id != "AE10",
#         stand_id != "AV02")

# Summarize growth for all trees
growth_summ <- growth_summary(growth_data)

# Remove trees for which annual growth could not be calculated (see likelihood
# model script for details)
growth_summ <- growth_summ[!is.na(growth_summ$size_corr_growth), ]

# Join growth and environment data
grow_env <- left_join(growth_summ, env_raw)

# Remove unneeded columns
grow_env <- grow_env[, c("tree_id", "size_corr_growth", "precip_mm", "temp_C")]

#==================
# New training data
#==================

# Create 3D array to store output
output <- array(NA, dim = c(6, 20, 11))

# Define radii and species to investigate
radius_list <- seq(1, 20, 1)
species_list <- c("PSME", "TSHE", "TSME", "THPL", "ABAM", "CANO")

# Loop through neighborhood sizes
for(i in radius_list){
  
  # Create neighborhoods
  neighbors <- graph_mat_all(mapping, radius = radius_list[i])
  
  # Remove focals within 20m of stand boundary
  neighbors <- neighbors %>% filter(x_coord >= 20 & x_coord <= 80 & y_coord >= 20 & y_coord <= 80)
  
  # Remove small competitors
  neighbors <- neighbors %>% filter(size_cat_comp == "regular")
  
  # Join neighborhoods with growth and environmental data
  full <- inner_join(neighbors, grow_env, by = "tree_id")
  
  # Exclude test data
  #full_sub1 <- full %>% filter(y_coord <= 40) # Training 1
  #full_sub2 <- full %>% filter(x_coord >= 60 & y_coord > 40) # Training 1
  #full_sub1 <- full %>% filter(x_coord <= 40) # Training 2
  #full_sub2 <- full %>% filter(y_coord <= 40 & x_coord > 40) # Training 2
  #full_sub1 <- full %>% filter(y_coord >= 60) # Training 3
  #full_sub2 <- full %>% filter(x_coord <= 40 & y_coord < 60) # Training 3
  full_sub1 <- full %>% filter(x_coord >= 60) # Training 4
  full_sub2 <- full %>% filter(y_coord >= 60 & x_coord < 60) # Training 4
  training <- rbind(full_sub1, full_sub2)
  
  # Loop through species
  for(j in 1:length(species_list)){
    
    # Subset to a single species
    sing_sp <- droplevels(training[training$species == species_list[j], ])
    
    # Find common competitors
    comps <- names(which(table(sing_sp$sps_comp) >= 100))
    
    # Change rare competitor species to OTHR
    sing_sp$sps_comp[which(sing_sp$sps_comp %in% comps == F)] <- "OTHR"
    
    # Define densities to include
    comps <- c(comps, "all")
    density_cols <- rep(NA, times = length(comps))
    for(c in 1:length(comps)){
      density_cols[c] <- grep(comps[c], names(sing_sp))
    }
    other_cols <- setdiff(grep("density", names(sing_sp)), density_cols)
    sing_sp$OTHR_density <- apply(sing_sp[, other_cols], 1, sum)
    
    # Convert competitor species to factor
    sing_sp$sps_comp <- as.factor(sing_sp$sps_comp)
    
    # Create matrices of factor variables
    if(length(unique(sing_sp$sps_comp)) > 1){
      dm_fac <- model.matrix(size_corr_growth ~ sps_comp,
                             sing_sp, contrasts.arg = list(sps_comp = contrasts(sing_sp$sps_comp, contrasts = F)))
    } else {
      dm_fac <- cbind(rep(1, times = nrow(sing_sp)), rep(1, times = nrow(sing_sp)))
      colnames(dm_fac) <- c("(Intercept)", paste("sps_comp", sing_sp$sps_comp[1], sep = ""))
    }
    
    # Combine factor predictors with contiunous predictors
    dm <- cbind(dm_fac, as.matrix(sing_sp[, c(9, 12, density_cols)]),
                as.matrix(sing_sp[, c("precip_mm", "temp_C", "OTHR_density")]))
    
    # Standardize variables except for first column (intercept)
    dm[, 2:ncol(dm)] <- apply(dm[, 2:ncol(dm)], 2, z_trans)
    
    # Change columns of NaNs (no variation) to zeros
    dm[, which(is.nan(dm[1, ]))] <- 0
    
    # Check for at least 3 observations per dataset - glmnet will fail otherwise
    if(nrow(dm) >= 30){  
      
      # Fit models
      mod <- cv.glmnet(dm, y = sing_sp$size_corr_growth, family = "gaussian")
      
      # Make predictions
      int_pred <- predict(mod, newx = dm, s = "lambda.1se")
      int_obs <- sing_sp %>% select(tree_id, size_corr_growth)
      comparison <- cbind(int_obs, int_pred)
      colnames(comparison)[2:3] <- c("obs", "pred")
      comparison <- comparison %>%
        group_by(tree_id) %>%
        summarize(observations = obs[1],
                  predictions = mean(pred))
      
      # Add radius to output
      output[j, i, 1] <- radius_list[i]
      
      # Record sample size (number of focal trees)
      output[j, i, 2] <- nrow(comparison)
      
      # Record number of competitor species variables retained
      output[j, i, 3] <- sum(coef(mod)[grep("sps_comp", rownames(coef(mod)))] != 0)
      
      # Record prox coefficient
      output[j, i, 4] <- coef(mod)[rownames((coef(mod))) == "prox"]
      
      # Record size_comp coefficient
      output[j, i, 5] <- coef(mod)[rownames((coef(mod))) == "size_comp_dbh"]
      
      # Record number of density variables retained
      output[j, i, 6] <- sum(coef(mod)[grep("density", rownames(coef(mod)))] != 0)
      
      # Return coefficient of determination
      output[j, i, 7] <- coef_det(comparison)
      
      # Calculate slope of observed growth vs. predicted growth
      slope_fit <- lm(observations ~ 0 + predictions, comparison)
      output[j, i, 8] <- coef(slope_fit)
      
      # Record mean square error of 1se model
      output[j, i, 9] <- mod$cvm[mod$lambda == mod$lambda.1se]
      
      # Record mean square error plus and minus 1 standard error
      output[j, i, 10] <- mod$cvup[mod$lambda == mod$lambda.1se]
      output[j, i, 11] <- mod$cvlo[mod$lambda == mod$lambda.1se]
    }
  }
}

# Divide output among focal species
PSME_output <- as.data.frame(output[1,,])
TSHE_output <- as.data.frame(output[2,,])
TSME_output <- as.data.frame(output[3,,])
THPL_output <- as.data.frame(output[4,,])
ABAM_output <- as.data.frame(output[5,,])
CANO_output <- as.data.frame(output[6,,])
names(PSME_output) <- c("radius", "num_focals", "comp_sps", "prox", "size_comp", "dens", "R_square", "slope",
                        "mse", "mse_plus_se", "mse_min_se")
names(TSHE_output) <- names(PSME_output)
names(TSME_output) <- names(PSME_output)
names(THPL_output) <- names(PSME_output)
names(ABAM_output) <- names(PSME_output)
names(CANO_output) <- names(PSME_output)

# Look at R square and mse patterns for each species
dat <- PSME_output
plot(R_square ~ radius, dat)
which.max(dat$R_square)
max(dat$R_square, na.rm = T)
plot(mse ~ radius, data = dat,
     ylim = c(min(dat$mse_min_se, na.rm = T), max(dat$mse_plus_se, na.rm = T)))
points(dat$radius, dat$mse_plus_se, col = "red")
points(dat$radius, dat$mse_min_se, col = "blue")
which.min(dat$mse)
which(dat$mse_min_se < dat$mse_plus_se[which.min(dat$mse)])

#==================
# Old training data
#==================

# Join growth and environment data
grow_env <- left_join(growth_summ, env_raw)

# Remove unneeded columns
grow_env <- grow_env[, c("tree_id", "size_corr_growth", "precip_mm", "temp_C")]

# Create 3D array to store output
output <- array(NA, dim = c(6, 20, 14))
#output <- array(NA, dim = c(6, 20, 11))

# Define radii and species to investigate
radius_list <- seq(1, 20, 1)
species_list <- c("PSME", "TSHE", "TSME", "THPL", "ABAM", "CANO")

# Loop through neighborhood sizes
for(i in radius_list){
  
  # Create neighborhoods
  neighbors <- graph_mat_all(mapping, radius = radius_list[i])
  
  # Remove focals within 20m of stand boundary
  neighbors <- neighbors %>% filter(x_coord >= 20 & x_coord <= 80 & y_coord >= 20 & y_coord <= 80)
  
  # Remove small competitors
  neighbors <- neighbors %>% filter(size_cat_comp == "regular")
  
  # Join neighborhoods with growth and environmental data
  full <- inner_join(neighbors, grow_env, by = "tree_id")
  
  # Convert competitor species to factor
  full$sps_comp <- as.factor(full$sps_comp)
  
  # INSERT HERE???????????
  
  # Loop through species
  for(j in 1:length(species_list)){
    
    # Subset to a single species
    sing_sp <- droplevels(full[full$species == species_list[j], ])
    
    # Split data into two cross-validation groups
    sing_sp$cv_group <- "A"
    sing_sp$cv_group[which((sing_sp$y_coord > 50 & sing_sp$x_coord <= 50) |
                             (sing_sp$y_coord < 50 & sing_sp$x_coord >= 50))] <- "B"
    
    # Subset data into two cross-validation groups
    a <- sing_sp[sing_sp$cv_group == "A", ]
    b <- sing_sp[sing_sp$cv_group == "B", ]
    
    # Create matrices of factor variables
    dm_fac_a <- model.matrix(size_corr_growth ~ sps_comp, a,
                             contrasts.arg = list(sps_comp = contrasts(a$sps_comp, contrasts = F)))
    dm_fac_b <- model.matrix(size_corr_growth ~ sps_comp, b,
                             contrasts.arg = list(sps_comp = contrasts(b$sps_comp, contrasts = F)))
    #dm_fac <- model.matrix(size_corr_growth ~ sps_comp, sing_sp,
    #                       contrasts.arg = list(sps_comp = contrasts(sing_sp$sps_comp, contrasts = F)))
    # Combine factor predictors with continous predictors (prox, size_comp_dbh, all_density, species-specific densities, precip, temp)
    dm_a <- cbind(dm_fac_a, as.matrix(a[c(9, 12, 14:31, 33:34)]))
    dm_b <- cbind(dm_fac_b, as.matrix(b[c(9, 12, 14:31, 33:34)]))
    #dm <- cbind(dm_fac, as.matrix(sing_sp[c(9, 12, 14:31, 33:34)]))
    
    # Standardize variables except for first column (intercept)
    dm_a[, 2:ncol(dm_a)] <- apply(dm_a[, 2:ncol(dm_a)], 2, z_trans)
    dm_b[, 2:ncol(dm_b)] <- apply(dm_b[, 2:ncol(dm_b)], 2, z_trans)
    #dm[, 2:ncol(dm)] <- apply(dm[, 2:ncol(dm)], 2, z_trans)
    
    # Change columns of NaNs (no variation) to zeros
    dm_a[, which(is.nan(dm_a[1, ]))] <- 0
    dm_b[, which(is.nan(dm_b[1, ]))] <- 0
    #dm[, which(is.nan(dm[1, ]))] <- 0
    
    # Check for at least 3 observations per dataset - glmnet will fail otherwise
    if(nrow(dm_a) >= 30 & nrow(dm_b) >= 30){
    #if(nrow(dm) >= 30){  
      
      # Fit models
      mod_a <- cv.glmnet(dm_a, y = a$size_corr_growth, family = "gaussian")
      mod_b <- cv.glmnet(dm_b, y = b$size_corr_growth, family = "gaussian")
      #mod <- cv.glmnet(dm, y = sing_sp$size_corr_growth, family = "gaussian")
      
      # Make predictions
      int_pred <- c(predict(mod_a, newx = dm_b, s = "lambda.1se"), predict(mod_b, newx = dm_a, s = "lambda.1se"))
      int_obs <- rbind(b %>% select(tree_id, size_corr_growth), a %>% select(tree_id, size_corr_growth))
      #int_pred <- predict(mod, newx = dm, s = "lambda.1se")
      #int_obs <- sing_sp %>% select(tree_id, size_corr_growth)
      comparison <- cbind(int_obs, int_pred)
      colnames(comparison)[2:3] <- c("obs", "pred")
      comparison <- comparison %>%
        group_by(tree_id) %>%
        summarize(observations = obs[1],
                  predictions = mean(pred))
      
      # Add radius to output
      output[j, i, 1] <- radius_list[i]
      
      # Record sample size
      output[j, i, 2] <- nrow(comparison)
      
      # Record number of competitor species variables retained
      output[j, i, 3] <- sum(coef(mod_a)[grep("sps_comp", rownames(coef(mod_a)))] != 0)
      #output[j, i, 3] <- sum(coef(mod)[grep("sps_comp", rownames(coef(mod)))] != 0)
      
      # Record prox coefficient
      output[j, i, 4] <- coef(mod_a)[rownames((coef(mod_a))) == "prox"]
      #output[j, i, 4] <- coef(mod)[rownames((coef(mod))) == "prox"]
      
      # Record size_comp coefficient
      output[j, i, 5] <- coef(mod_a)[rownames((coef(mod_a))) == "size_comp_dbh"]
      #output[j, i, 5] <- coef(mod)[rownames((coef(mod))) == "size_comp_dbh"]
      
      # Record number of density variables retained
      output[j, i, 6] <- sum(coef(mod_a)[grep("density", rownames(coef(mod_a)))] != 0)
      #output[j, i, 6] <- sum(coef(mod)[grep("density", rownames(coef(mod)))] != 0)
      
      # Return coefficient of determination
      output[j, i, 7] <- coef_det(comparison)
      
      # Calculate slope of observed growth vs. predicted growth
      slope_fit <- lm(observations ~ 0 + predictions, comparison)
      output[j, i, 8] <- coef(slope_fit)
      
      # Record mean square error of 1se model
      output[j, i, 9] <- mod_a$cvm[mod_a$lambda == mod_a$lambda.1se]
      output[j, i, 10] <- mod_b$cvm[mod_b$lambda == mod_b$lambda.1se]
      #output[j, i, 9] <- mod$cvm[mod$lambda == mod$lambda.1se]
      
      # Record mean square error plus and minus 1 standard error
      output[j, i, 11] <- mod_a$cvup[mod_a$lambda == mod_a$lambda.1se]
      output[j, i, 12] <- mod_a$cvlo[mod_a$lambda == mod_a$lambda.1se]
      output[j, i, 13] <- mod_b$cvup[mod_b$lambda == mod_b$lambda.1se]
      output[j, i, 14] <- mod_b$cvlo[mod_b$lambda == mod_b$lambda.1se]
      #output[j, i, 10] <- mod$cvup[mod$lambda == mod$lambda.1se]
      #output[j, i, 11] <- mod$cvlo[mod$lambda == mod$lambda.1se]
    }
  }
}

# Divide output among focal species
PSME_output <- as.data.frame(output[1,,])
TSHE_output <- as.data.frame(output[2,,])
TSME_output <- as.data.frame(output[3,,])
THPL_output <- as.data.frame(output[4,,])
ABAM_output <- as.data.frame(output[5,,])
CANO_output <- as.data.frame(output[6,,])
names(PSME_output) <- c("radius", "num_focals", "comp_sps", "prox", "size_comp", "dens", "R_square", "slope",
                        "mse_a", "mse_b", "mse_plus_se_a", "mse_min_se_a", "mse_plus_se_b", "mse_min_se_b")
#names(PSME_output) <- c("radius", "num_focals", "comp_sps", "prox", "size_comp", "dens", "R_square", "slope",
#                        "mse", "mse_plus_se", "mse_min_se")
names(TSHE_output) <- names(PSME_output)
names(TSME_output) <- names(PSME_output)
names(THPL_output) <- names(PSME_output)
names(ABAM_output) <- names(PSME_output)
names(CANO_output) <- names(PSME_output)

# Look at R square and mse patterns for one species
dat <- ABAM_output
plot(R_square ~ radius, dat)
plot(mse_a ~ radius, data = dat,
     ylim = c(min(c(dat$mse_min_se_a, dat$mse_min_se_b), na.rm = T),
              max(c(dat$mse_plus_se_a, dat$mse_plus_se_b), na.rm = T)))
points(dat$radius, dat$mse_plus_se_a, col = "red")
points(dat$radius, dat$mse_min_se_a, col = "blue")
points(dat$radius, dat$mse_b, col = "black", pch = 4)
points(dat$radius, dat$mse_plus_se_b, col = "red", pch = 4)
points(dat$radius, dat$mse_min_se_b, col = "blue", pch = 4)

# Subset to neighborhood sizes that give model with positive R-square and
# that retains at least one neighborhood variable
posr <- dat[dat$R_square > 0 & (dat$comp_sps > 0 | dat$prox > 0 | dat$size_comp > 0 | dat$dens > 0), ]
posr$radius[which.max(posr$R_square)]
max(posr$R_square, na.rm = T)
#posr$radius[which.min(posr$mse_a)]
#posr$radius[which(posr$mse_min_se_a < posr$mse_plus_se_a[which.min(posr$mse_a)])]
#posr$radius[which.min(posr$mse_b)]
#posr$radius[which(posr$mse_min_se_b < posr$mse_plus_se_b[which.min(posr$mse_b)])]

#plot(mse ~ radius, data = dat,
#     ylim = c(min(dat$mse_min_se), max(dat$mse_plus_se)))
#points(dat$radius, dat$mse_plus_se, col = "red")
#points(dat$radius, dat$mse_min_se, col = "blue")
#which.min(dat$mse)
#which(dat$mse_min_se < dat$mse_plus_se[which.min(dat$mse)])


#================================================================
# Old Method: Lasso regression with all densities in single model
#================================================================

# This section makes a glmnet model of size corrected growth predicted by focal
# species identity, stand id, and neighborhood density (all neighbor species 
# combined) measured for neighborhood sizes 1-20m. The coefficients of each of
# the density measurements are plotted, with the assumption being that the
# neighborhood radius producing the density measurement with the highest 
# coefficient has the most explanatory power and should therefore be used.

# Create vector of potential neighborhood radii
radius_list <- seq(1, 20, 1)

# Create data frame to add densities to
density <- mapping[, c("tree_id", "x_coord", "y_coord")]

# Fill matrix with densities by using density_all_stands function with each radius size
for(i in 1:length(radius_list)){
  new_dens <- density_all_stands(mapping, radius = radius_list[i])[ , c("tree_id", "all_density")]
  colnames(new_dens)[2] <- paste("all_density", radius_list[i], sep = "")
  density <- left_join(density, new_dens)
}

# Attach environmental data to density
dens_env <- left_join(density, env_dat)

# Attach density measurements to growth
growth <- left_join(growth_summ, dens_env)

# Subset to single species
growth_sub <- growth[growth$species == "CANO", ]

# Subset to rows where all densities defined (no NAs)
check_no_na <- function(x){all(!is.na(x))}
growth_sub <- growth_sub[apply(growth_sub, 1, check_no_na),]

# Create model matrix of continuous predictor variables
dm <- model.matrix(size_corr_growth ~ precip_mm + temp_C + all_density1 + all_density2 +
                     all_density3 + all_density4 + all_density5 + all_density6 + all_density7 +
                     all_density8 + all_density9 + all_density10 + all_density11 + all_density12 +
                     all_density13 + all_density14 + all_density15 + all_density16 + all_density17 +
                     all_density18 + all_density19 + all_density20, growth_sub)
dm[,2:ncol(dm)] <- apply(dm[,2:ncol(dm)], 2, z_trans)

# Fit and plot model
mod <- cv.glmnet(dm, y = growth_sub$size_corr_growth, family = "gaussian")
plot(mod)

# View model coefficients
coef(mod)

# Return non-zero density coefficients and the one with highest magnitude
dens_coef <- coef(mod)[grep("density", rownames(coef(mod)))]
which(dens_coef != 0)
which.max(abs(dens_coef))