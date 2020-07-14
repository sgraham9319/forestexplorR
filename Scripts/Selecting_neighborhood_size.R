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
growth_data <- read.csv("../Data/Cleaned_tree_growth_2017.csv", stringsAsFactors = F)

# Load mapping data for individual trees
mapping <- read.csv("../Data/Cleaned_mapping_2017.csv", stringsAsFactors = F)

# Load environmental data
env_dat <- read.csv("../Data/stand_abiotic_data.csv", stringsAsFactors = F)

# Create combined mapping and environmental dataset
env_dat <- mapping %>% left_join(env_dat) %>% select(tree_id, precip_mm, temp_C, elev_m)

# Remove test data (2017 measurements and stands TO04, AE10, and AV02)
growth_data <- growth_data %>%
  filter(year != 2017,
         stand_id != "TO04",
         stand_id != "AE10",
         stand_id != "AV02")

# Remove test data from mapping dataset (stands TO04, AE10, and AV02)
mapping <- mapping %>%
  filter(stand_id != "TO04",
         stand_id != "AE10",
         stand_id != "AV02")

# Summarize growth for all trees
growth_summ <- growth_summary(growth_data)

# Remove trees for which annual growth could not be calculated (see likelihood
# model script for details)
growth_summ <- growth_summ[!is.na(growth_summ$size_corr_growth), ]

# Calculate radial growth to be consistent with Fortunel and Canham
growth_summ$radial_growth <- growth_summ$annual_growth / 2

#============================================================
# Part 2. Comparing models using different neighborhood sizes
#============================================================

# The below function takes a list of neighborhood radii to try and creates a glmnet
# model for each neighborhood radius. The dependent variable is size corrected 
# growth and the independent variables are focal species, stand id, and density
# (all neighbor species combined) calculated using the given radius
glmnet_mods <- function(mapping, growth_summ, env_dat, species, radius_list){
  
  # Create output table
  output <- data.frame(matrix(NA, ncol = 7, nrow = 1))
  colnames(output) <- c("radius", "log_lambda", "log_1se_lambda", "mean_sq_err", "mean_plus_sd",
                        "mean_minus_sd", "den_coef")
  
  # Subset growth data to focal species
  growth_sub <- growth_summ[growth_summ$species == species, ]
  
  for(i in 1:length(radius_list)){
    
    # Calculate neighborhood density for all trees
    densities <- density_all_stands(mapping, radius = radius_list[i])
    
    # Join environmental and density data
    dens_env <- left_join(densities, env_dat)
    
    # Attach neighborhood information to growth
    growth <- left_join(growth_sub, dens_env)
    
    # Exclude missing growth and density values
    growth <- growth %>%
      filter(!is.na(size_corr_growth) &
               !is.na(all_density))
    
    # Create model matrix of continuous predictor variables
    dm <- model.matrix(size_corr_growth ~ growth$precip_mm + growth$temp_C + growth$all_density, growth)
    colnames(dm)[2:ncol(dm)] <- c("precip", "temp", "density")
    dm[,2:ncol(dm)] <- apply(dm[,2:ncol(dm)], 2, z_trans)
    
    # Fit model
    mod <- cv.glmnet(dm, y = growth$size_corr_growth, family = "gaussian")
    
    # Check if density retained in lambda.1se model
    if(coef(mod)[rownames(coef(mod)) == "density"] != 0){
      
      # If density reatined, extract model information
      radius <- rep(radius_list[i], times = length(mod$lambda))
      log_lambda <- log(mod$lambda)
      log_1se_lambda <- log(mod$lambda.1se)
      mean_sq_err <- mod$cvm
      mean_plus_sd <- mod$cvup
      mean_minus_sd <- mod$cvlo
      den_coef <- coef(mod)[rownames(coef(mod)) == "density"]
      loop_output <- data.frame(radius, log_lambda, log_1se_lambda, mean_sq_err, mean_plus_sd,
                                mean_minus_sd, den_coef)
      
      # Add loop output to overall output
      output <- rbind(output, loop_output)
      }
    }
  # Return output
  output[-1,]
}

# Create models for a range of neighborhood radii
rad_var <- glmnet_mods(mapping, growth_summ, env_dat, species = "CANO", radius_list = seq(2, 20, 2))

# Plot mean and mean + sd error of 1se model for each radius
mod_summ_1se <- rad_var %>% filter(log_lambda == log_1se_lambda)
plot(mean_sq_err ~ radius, data = mod_summ_1se,
     ylim = c(min(mod_summ_1se$mean_sq_err), max(mod_summ_1se$mean_plus_sd)))
points(mod_summ_1se$radius, mod_summ_1se$mean_plus_sd, col = "red")

# Check for obvious relationships between mean square error and radius size
rad_var$radius_f <- as.factor(rad_var$radius)
ggplot(rad_var, aes(x = log_lambda, y = mean_sq_err, col = radius_f)) +
  geom_point()

# Check for relationship between density model coefficient and radius size
plot(den_coef ~ radius, data = rad_var, xlab = "Neighborhood radius (m)",
     ylab = "Density coefficient")

# Compare pairs of radii more carefully
pair_comp <- function(rad_var, rad_a, rad_b){
  rad_a_sub <- rad_var %>% filter(radius == rad_a)
  rad_b_sub <- rad_var %>% filter(radius == rad_b)
  plot(mean_sq_err ~ log_lambda, data = rad_a_sub, col = "red")
  points(rad_a_sub$log_lambda, rad_a_sub$mean_plus_sd, col = "red")
  points(rad_b_sub$log_lambda, rad_b_sub$mean_sq_err, col = "blue")
  points(rad_b_sub$log_lambda, rad_b_sub$mean_plus_sd, col = "blue")
  legend("topleft", legend = c(rad_a, rad_b), fill = c("red", "blue"))
}
pair_comp(rad_var, 16, 20)

#==================================================
# Part 3. Including all densities in a single model
#==================================================

# This section makes a glmnet model of size corrected growth predicted by focal
# species identity, stand id, and neighborhood density (all neighbor species 
# combined) measured for neighborhood sizes 1-20m. The coefficients of each of
# the density measurements are plotted, with the assumption being that the
# neighborhood radius producing the density measurement with the highest 
# coefficient has the most explanatory power and should therefore be used.

# Create vector of potential neighborhood radii
radius_list <- seq(2, 20, 2)

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
dm <- model.matrix(size_corr_growth ~ precip_mm + temp_C + all_density2 + all_density4 +
                     all_density6 + all_density8 + all_density10 + all_density12 + all_density14 +
                     all_density16 + all_density18 + all_density20, growth_sub)
dm[,2:ncol(dm)] <- apply(dm[,2:ncol(dm)], 2, z_trans)

# Fit and plot model
mod <- cv.glmnet(dm, y = growth_sub$size_corr_growth, family = "gaussian")
plot(mod)
coef(mod) # No density measurements retained in 1se model for any species
# PSME: 10, 14, 20; TSHE: 12*; TSME: 14*, 20*; THPL: 12*, 16*; ABAM: 2*, 4*, 8, 10*, 14*; CANO: none