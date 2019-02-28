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
densities <- density_all_stands(mapping, radius = 10)

# Attach neighborhood information to growth
growth <- left_join(growth_summ, densities)

# Check if Gaussian
kurtosis(growth$size_corr_growth, na.rm = T)

# Exclude missing growth and density values
growth <- growth %>%
  filter(!is.na(size_corr_growth) &
           !is.na(all_density))

# Change stand_id and species to factors
growth$stand_id <- as.factor(growth$stand_id)
growth$species <- as.factor(growth$species)

#=======================================================
# Creating full model to compare with graph matrix model
#=======================================================

# Create model matrix of factorized predictor variables
dm_factors <- model.matrix(size_corr_growth ~ stand_id + species, growth,
                           contrasts.arg = list(stand_id = contrasts(growth$stand_id, contrasts = F),
                                                species = contrasts(growth$species, contrasts = F)))

# Combine factorized predictors with continuous predictors
dm <- cbind(dm_factors, growth$all_density, growth$ABAM_density,
            growth$TSHE_density, growth$PSME_density)
colnames(dm)[(ncol(dm_factors) + 1):ncol(dm)] <- c("all_density", "ABAM_density",
                                                   "TSHE_density", "PSME_density")

# Standardize variables
dm[,2:ncol(dm)] <- apply(dm[,2:ncol(dm)], 2, z_trans)

# Fit model
glmmod <- cv.glmnet(dm, y = growth$size_corr_growth, family = "gaussian")

predictions <- predict(glmmod, newx = dm, s = "lambda.1se")
observations <- growth$size_corr_growth
comparison <- data.frame(predictions, observations)
ggplot(comparison, aes(x = observations, y = X1)) +
  geom_hex() +
  theme_bw() +
  ylim(0, 0.3) +
  geom_abline(intercept = 0, slope = 1)

coef_det(comparison) # 0.260, 0.262, 0.253, 0.253, 0.253

#====================================================================
# Create simple model with no species-specific competitor information
#====================================================================

# Create model matrix of factorized predictor variables
dm_factors <- model.matrix(size_corr_growth ~ stand_id + species, growth,
                           contrasts.arg = list(stand_id = contrasts(growth$stand_id, contrasts = F),
                                                species = contrasts(growth$species, contrasts = F)))

# Combine factorized predictors with continuous predictors
dm <- cbind(dm_factors, growth$all_density)
colnames(dm)[(ncol(dm_factors) + 1):ncol(dm)] <- c("all_density")

# Standardize variables
dm[,2:ncol(dm)] <- apply(dm[,2:ncol(dm)], 2, z_trans)

# Fit model
glmmod <- cv.glmnet(dm, y = growth$size_corr_growth, family = "gaussian")

predictions <- predict(glmmod, newx = dm, s = "lambda.1se")
observations <- growth$size_corr_growth
comparison <- data.frame(predictions, observations)
ggplot(comparison, aes(x = observations, y = X1)) +
  geom_hex() +
  theme_bw() +
  ylim(0, 0.3) +
  geom_abline(intercept = 0, slope = 1)

coef_det(comparison) # 0.247, 0.243, 0.253, 0.239, 0.239


#=============================
# Creating basic model 2/26/19
#=============================

# Subset to common species
comm_sps <- growth %>%
  filter(species %in% c("ABAM", "PSME", "TSHE"))
comm_sps <- droplevels(comm_sps)

# Explore how these species are distributed among stands
table(comm_sps[, c("species", "stand_id")])

# Create model matrix of factorized predictor variables
dm_factors <- model.matrix(size_corr_growth ~ stand_id + species, comm_sps,
                           contrasts.arg = list(stand_id = contrasts(comm_sps$stand_id, contrasts = F),
                                                species = contrasts(comm_sps$species, contrasts = F)))

# Combine factorized predictors with continuous predictors
dm <- cbind(dm_factors, comm_sps$all_density, comm_sps$ABAM_density,
            comm_sps$TSHE_density, comm_sps$PSME_density)
colnames(dm)[(ncol(dm_factors) + 1):ncol(dm)] <- c("all_density", "ABAM_density",
                                                   "TSHE_density", "PSME_density")

# Standardize variables
z_trans <- function(x){
  (x - mean(x)) / sd(x)
}
dm[,2:ncol(dm)] <- apply(dm[,2:ncol(dm)], 2, z_trans)

# Fit model
glmmod <- cv.glmnet(dm, y = comm_sps$size_corr_growth, family = "gaussian")

# Look at coefficients
coef(glmmod, s = "lambda.min")
coef(glmmod, s = "lambda.1se")


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
# Minus the mean and divide by sd
a <- matrix(1:25, ncol = 5, nrow = 5)
a[,2:3] <- apply(a[, 2:3], 2, z_trans)
z_trans <- function(x){
  (x - mean(x)) / sd(x)
}

x <- glmnet_mods(mapping, growth_summ, radius_list = c(1, 10, 20))

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
    dm <- cbind(dm_factors, growth$all_density)
    dm[,2:ncol(dm)] <- apply(dm[,2:ncol(dm)], 2, z_trans)
    
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

x$radius_f <- as.factor(x$radius)
ggplot(x, aes(x = log_lambda, y = mean_sq_err, col = radius_f)) +
  geom_point()

x10 <- x %>% filter(radius == 10)
plot(x10$mean_sq_err ~ x10$log_lambda, col = "red")
points(x = x10$log_lambda, y = x10$mean_plus_sd, col = "red")
x20 <- x %>% filter(radius == 20)
points(x = x20$log_lambda, y = x20$mean_sq_err, col = "blue")
points(x = x20$log_lambda, y = x20$mean_plus_sd, col = "blue")
x1 <- x %>% filter(radius == 1)
points(y = x1$mean_sq_err, x = x1$log_lambda)
points(x = x1$log_lambda, y = x1$mean_plus_sd)


hist(x$mean_plus_sd - x$mean_sq_err)
levels(x$radius_f)

# x = 1, y = 10, nothing = 20
density1 <- density_all_stands(mapping, radius = 1) 
density10 <- density_all_stands(mapping, radius = 10)
density20 <- density_all_stands(mapping, radius = 20)
dens <- merge(density1, density10, by = "tree_id")
dens <- merge(dens, density20, by = "tree_id")

# Attach neighborhood information to growth
growth <- left_join(growth_summ, dens)

# Exclude missing growth and density values
growth <- growth %>%
  filter(!is.na(size_corr_growth) &
           !is.na(all_density) &
           !is.na(all_density.x) &
           !is.na(all_density.y))

# Change stand_id and species to factors
growth$stand_id <- as.factor(growth$stand_id)
growth$species <- as.factor(growth$species)

# First create a model matrix of factorized predictor variables
dm_factors <- model.matrix(size_corr_growth ~ stand_id + species, growth,
                           contrasts.arg = list(stand_id = contrasts(growth$stand_id, contrasts = F),
                                                species = contrasts(growth$species, contrasts = F)))

# Combine matrix of factors with matrix of continuous predictor variables
dm <- cbind(dm_factors, growth$all_density, growth$all_density.x, growth$all_density.y)
dm[,2:ncol(dm)] <- apply(dm[,2:ncol(dm)], 2, z_trans)

# Fit model
mod <- cv.glmnet(dm, y = growth$size_corr_growth, family = "gaussian")
plot(mod)
# Get coefficients for minimum lambda model
coef(mod, s = "lambda.min")
# Get coefficients of model corresponding to dotted line on right of plot
coef(mod, s = "lambda.1se")

#===========================================
# Finding most informative neighborhood size
#===========================================
radius_list <- seq(1, 20, 1)

density <- data.frame(matrix(NA, ncol = length(radius_list), nrow = nrow(mapping)))
colnames(density) <- paste("all_density", radius_list, sep = "")
for(i in 1:length(radius_list)){
  density[, i] <- density_all_stands(mapping, radius = radius_list[i])$all_density
}

# Add tree-id column for joining
density$tree_id <- mapping$tree_id

# Attach neighborhood information to growth
growth <- left_join(growth_summ, density)

# Exclude missing growth and density values
growth <- growth %>%
  filter(!is.na(size_corr_growth) &
           !is.na(all_density1))

# Change stand_id and species to factors
growth$stand_id <- as.factor(growth$stand_id)
growth$species <- as.factor(growth$species)

# Create a model matrix of factorized predictor variables
dm_factors <- model.matrix(size_corr_growth ~ stand_id + species, growth,
                           contrasts.arg = list(stand_id = contrasts(growth$stand_id, contrasts = F),
                                                species = contrasts(growth$species, contrasts = F)))

# Identify columns of growth containing continuous predictor variables
predictors <- (ncol(growth) - (length(radius_list) - 1)):ncol(growth)

# Combine factor and continuous predictor variables in design matrix
dm <- cbind(dm_factors, as.matrix(growth[, predictors]))

# Standardize predictor variables
dm[,2:ncol(dm)] <- apply(dm[,2:ncol(dm)], 2, z_trans)

# Fit and plot model
mod <- cv.glmnet(dm, y = growth$size_corr_growth, family = "gaussian")
plot(mod)

# Get coefficients for minimum lambda model
desired_coef <- (length(coef(mod, s = 0)) - (length(radius_list) - 1)):length(coef(mod, s = 0))
coefficient <- coef(mod, s = 0)[desired_coef]
coef_table <- data.frame(radius_list, coefficient)

# Plot magnitude of density coefficient vs. neighborhood radius
plot(abs(coef_table$coefficient) ~ coef_table$radius_list,
     ylab = "Coefficient magnitude", xlab = "Neighborhood radius")
