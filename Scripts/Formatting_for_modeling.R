# Load package
devtools::load_all()

library(moments)
library(glmnet)
library(glmnetUtils)
library(ggplot2)

#==============
# Cleaning data
#==============

# Load growth data
growth_data <- read.csv("../Data/Tree_growth_2017.csv", stringsAsFactors = F)

# Remove test data (2017 measurements and stands TO04, AE10, and AV02)
growth_data <- growth_data %>%
  filter(year != 2017,
         stand_id != "TO04",
         stand_id != "AE10",
         stand_id != "AV02")

# Remove rows without dbh measurement
growth_data <- growth_data[!is.na(growth_data$dbh), ]

# Identify trees that were never measured as 15 cm dbh or larger
small_trees_data <- growth_data %>% group_by(tree_id) %>% 
  summarize(max_size = max(dbh)) %>%
  filter(max_size < 15)
small_trees <- unique(small_trees_data$tree_id)

# Exclude small trees
growth_data <- growth_data[-which(growth_data$tree_id %in% small_trees), ]

# Load mapping data for individual trees
mapping <- read.csv("../Data/Mapping_2017.csv", stringsAsFactors = F)

# Exclude small trees from mapping data
mapping <- mapping[-which(mapping$tree_id %in% small_trees), ]

#=================
# Summarizing data 
#=================

# Summarize growth for all trees
growth_summ <- growth_summary(growth_data)

# Calculate neighborhood density for all trees
densities <- density_all_stands(mapping)

# Attach neighborhood information to growth
growth <- left_join(growth_summ, densities)


# SHOULD WE EXCLUDE UNCOMMON TREE SPECIES?

# Check if Gaussian
kurtosis(growth$size_corr_growth, na.rm = T)

# Exclude missing growth and density values
growth1 <- growth %>%
  filter(!is.na(size_corr_growth) &
           !is.na(all_density))

# Change stand_id and species to factors
growth1$stand_id <- as.factor(growth1$stand_id)
growth1$species <- as.factor(growth1$species)

#=======================
# Creating design matrix
#=======================

# First create a model matrix of factorized predictor variables
dm_factors <- model.matrix(growth1$size_corr_growth ~ growth1$stand_id +
                             growth1$species)[, -1]

# Combine matrix of factors with matrix of continuous predictor variables
# Only all_density
dm <- as.matrix(data.frame(growth1$all_density, dm_factors))
# All species-specific densities
#dm <- as.matrix(data.frame(growth1[, 11:ncol(growth1)], dm_factors))

# Fit model
glmmod <- cv.glmnet(dm, y = growth1$size_corr_growth, family = "gaussian")

# Plot model object - lower mean-squared error means better fit and the numbers
# on the top show the number of variables included. Region delimited by dotted
# vertical lines represents equally good model fit
plot(glmmod)

# For each model I want: log(lambda), cvm, cvup, cvlo

# Get coefficients for minimum lambda model
coef(glmmod, s = "lambda.min")

# Get coefficients of model corresponding to dotted line on right of plot
coef(glmmod, s = "lambda.1se")

#=================================
# Model with tshe-specific density
#=================================

# Create design matrix
dm_factors <- model.matrix(growth1$size_corr_growth ~ growth1$stand_id +
                             growth1$species)[, -1]
dm <- as.matrix(data.frame(growth1$all_density, growth1$tshe_density, dm_factors))

# Fit model
glmmod <- cv.glmnet(dm, y = growth1$size_corr_growth, family = "gaussian")

# Plot model
plot(glmmod)

# Get coefficients for minimum lambda model
coef(glmmod, s = "lambda.min")

# Get coefficients of model corresponding to dotted line on right of plot
coef(glmmod, s = "lambda.1se")

#=============================================
# Try including all species-specific densities
#=============================================

# Create design matrix
dm_factors <- model.matrix(growth1$size_corr_growth ~ growth1$stand_id +
                             growth1$species)[, -1]
dm <- as.matrix(data.frame(growth1$all_density, growth1$tshe_density,
                           growth1$abam_density, growth1$thpl_density,
                           growth1$tsme_density, growth1$cano_density,
                           growth1$pico_density, growth1$psme_density,
                           dm_factors))

# Fit model
glmmod <- cv.glmnet(dm, y = growth1$size_corr_growth, family = "gaussian")

# Plot model
plot(glmmod)

# Get coefficients for minimum lambda model
coef(glmmod, s = "lambda.min")

# Get coefficients of model corresponding to dotted line on right of plot
coef(glmmod, s = "lambda.1se")



# Create design matrix
dm_factors <- model.matrix(growth1$size_corr_growth ~ growth1$stand_id +
                             growth1$species)[, -1]
dm <- as.matrix(data.frame(growth1$tshe_density,
                           growth1$abam_density, growth1$thpl_density,
                           growth1$tsme_density, growth1$cano_density,
                           growth1$pico_density, growth1$psme_density,
                           dm_factors))

# Fit model
glmmod <- cv.glmnet(dm, y = growth1$size_corr_growth, family = "gaussian")

# Plot model
plot(glmmod)

# Get coefficients for minimum lambda model
coef(glmmod, s = "lambda.min")

# Get coefficients of model corresponding to dotted line on right of plot
coef(glmmod, s = "lambda.1se")


# Compare model predictions to observed data
predictions <- predict(glmmod, newx = dm, s = "lambda.1se")
observations <- growth1$size_corr_growth
comparison <- data.frame(predictions, observations)
ggplot(comparison, aes(x = observations, y = X1)) +
  geom_hex() +
  theme_bw() +
  ylim(0, 0.3) +
  geom_abline(intercept = 0, slope = 1)


growth1 %>%
  group_by(stand_id) %>%
  summarize(num_tshe = sum(species == "TSHE"),
            num_abam = sum(species == "ABAM"),
            num_psme = sum(species == "PSME"))


#====================
# Metaparameter sweep
#====================


x <- glmnet_mods(mapping, growth_summ, radius_list = seq(5, 100, 5))

glmnet_mods <- function(mapping, growth_summ, radius_list){
  
  # Create output table
  output <- data.frame(matrix(NA, ncol = 5, nrow = 1))
  colnames(output) <- c("radius", "log_lambda", "mean_sq_err", "mean_plus_sd",
                        "mean_minus_sd")
  
  for(i in 1:length(radius_list)){
    
    # Calculate neighborhood density for all trees
    densities <- density_all_stands(mapping, radius = radius_list[i])
    
    # Attach neighborhood information to growth
    growth <- left_join(growth_summ, densities)
    
    # Exclude missing growth and density values
    growth <- growth %>%
      filter(!is.na(size_corr_growth) &
               !is.na(all_density))
    
    # Change stand_id and species to factors
    growth$stand_id <- as.factor(growth$stand_id)
    growth$species <- as.factor(growth$species)
    
    # First create a model matrix of factorized predictor variables
    dm_factors <- model.matrix(size_corr_growth ~ stand_id + species, growth,
                               contrasts.arg = list(stand_id = contrasts(growth$stand_id, contrasts = F),
                                                    species = contrasts(growth$species, contrasts = F)))
    
    # Combine matrix of factors with matrix of continuous predictor variables
    dm <- as.matrix(data.frame(growth$all_density, dm_factors))
    
    # Fit model
    mod <- cv.glmnet(dm, y = growth$size_corr_growth, family = "gaussian")
    
    # Extract model information
    radius <- rep(radius_list[i], times = length(mod$lambda))
    log_lambda <- log(mod$lambda)
    mean_sq_err <- mod$cvm
    mean_plus_sd <- mod$cvup
    mean_minus_sd <- mod$cvlo
    loop_output <- data.frame(radius, log_lambda, mean_sq_err, mean_plus_sd,
                              mean_minus_sd)
    
    # Add loop output to overall output
    output <- rbind(output, loop_output)
  }
  # Return output
  output[-1,]
}

plot(mean_sq_err ~ log_lambda, col = radius, data = x)
