#########################################
# Single-species directional graph models
#########################################

# Author: Stuart Graham
# Created: 4/27/2020
# Last edited: 6/9/2020

# Load TreeNeighborhood package
devtools::load_all()

# Load other required packages
library(glmnet)
library(ggplot2)

######################
# Part 1. Loading data
######################

# Load mapping data
mapping <- read.csv("Data/Cleaned_mapping_2017.csv", stringsAsFactors = F)

# Load tree measurement data
tree <- read.csv("Data/Cleaned_tree_growth_2017.csv", stringsAsFactors = F)

#############################
# Part 2. Excluding test data
#############################

# Remove test data from tree dataset (2017 measurements and stands TO04, AE10, and AV02)
tree <- tree %>%
  filter(year != 2017,
         stand_id != "TO04",
         stand_id != "AE10",
         stand_id != "AV02")

# Remove test data from mapping dataset (stands TO04, AE10, and AV02)
mapping <- mapping %>%
  filter(stand_id != "TO04",
         stand_id != "AE10",
         stand_id != "AV02")

################################
# Part 3. Creating neighborhoods
################################

# Obtain all neighborhood data - different radius for each species
#neighbors <- graph_mat_all(mapping, radius = 10)
neighbors_PSME <- graph_mat_all(mapping, radius = 15)
neighbors_TSHE <- graph_mat_all(mapping, radius = 12)
neighbors_TSME <- graph_mat_all(mapping, radius = 15)
neighbors_THPL <- graph_mat_all(mapping, radius = 16)
neighbors_ABAM <- graph_mat_all(mapping, radius = 3)
neighbors_CANO <- graph_mat_all(mapping, radius = 7)

# Subset to single species of focal
neighbors_PSME <- neighbors_PSME %>% filter(species == "PSME")
neighbors_TSHE <- neighbors_TSHE %>% filter(species == "TSHE")
neighbors_TSME <- neighbors_TSME %>% filter(species == "TSME")
neighbors_THPL <- neighbors_THPL %>% filter(species == "THPL")
neighbors_ABAM <- neighbors_ABAM %>% filter(species == "ABAM")
neighbors_CANO <- neighbors_CANO %>% filter(species == "CANO")

# Remove focals whose neighborhood overlaps stand boundary
#neighbors <- neighbors %>% filter(x_coord >= 10 & x_coord <= 90 & y_coord >= 10 & y_coord <= 90)
neighbors_PSME <- neighbors_PSME %>% filter(x_coord >= 15 & x_coord <= 85 & y_coord >= 15 & y_coord <= 85)
neighbors_TSHE <- neighbors_TSHE %>% filter(x_coord >= 12 & x_coord <= 88 & y_coord >= 12 & y_coord <= 88)
neighbors_TSME <- neighbors_TSME %>% filter(x_coord >= 15 & x_coord <= 85 & y_coord >= 15 & y_coord <= 85)
neighbors_THPL <- neighbors_THPL %>% filter(x_coord >= 16 & x_coord <= 84 & y_coord >= 16 & y_coord <= 84)
neighbors_ABAM <- neighbors_ABAM %>% filter(x_coord >= 3 & x_coord <= 97 & y_coord >= 3 & y_coord <= 97)
neighbors_CANO <- neighbors_CANO %>% filter(x_coord >= 7 & x_coord <= 93 & y_coord >= 7 & y_coord <= 93)

# Recombine single species neighborhoods
neighbors <- rbind(neighbors_PSME, neighbors_TSHE, neighbors_TSME, neighbors_THPL, neighbors_ABAM, neighbors_CANO)

# Remove small competitors
neighbors <- neighbors %>% filter(size_cat_comp == "regular")
#neighbors_PSME <- neighbors_PSME %>% filter(size_cat_comp == "regular")
#neighbors_TSHE <- neighbors_TSHE %>% filter(size_cat_comp == "regular")
#neighbors_TSME <- neighbors_TSME %>% filter(size_cat_comp == "regular")
#neighbors_THPL <- neighbors_THPL %>% filter(size_cat_comp == "regular")
#neighbors_ABAM <- neighbors_ABAM %>% filter(size_cat_comp == "regular")
#neighbors_CANO <- neighbors_CANO %>% filter(size_cat_comp == "regular")

###################################
# Part 4. Calculating annual growth
###################################

# Calculate annual growth for all trees
growth <- growth_summary(tree)

# Remove trees for which annual growth could not be calculated (see likelihood
# model script for details)
growth <- growth[!is.na(growth$size_corr_growth), ]

# Check data one more time
sum(is.nan(growth$size_corr_growth))

# Calculate radial growth to be consistent with Fortunel and Canham
growth$radial_growth <- growth$annual_growth / 2

###############################################################
# Part 5. Combining growth, neighborhood and environmental data
###############################################################

# Extract required columns from growth data
growth_cols <- growth[, c("tree_id", "midpoint_size", "radial_growth", "size_corr_growth")]

# Join growth and neighborhood data. Use inner join because there will be
# no growth data for focals measured only once or with negative growth
full <- inner_join(neighbors, growth_cols, by = "tree_id")

# Change stand_id column to factor
full$stand_id <- as.factor(full$stand_id)

# Load environmental data
env_dat <- read.csv("Data/stand_abiotic_data.csv", stringsAsFactors = F)

# Add environmental data to neighborhoods and growth
full <- left_join(full, env_dat)

#################################
# Part 6. Creating model matrices
#################################

# Subset to a single species
sing_sp <- droplevels(full[full$species == "CANO", ])

# Remove growth rate outliers (be sure to remove same data (i.e. size_corr_growth) from lkhd model)
#growth_lim <- mean(sing_sp$size_corr_growth) + (3 * sd(sing_sp$size_corr_growth))
#sing_sp <- sing_sp[sing_sp$size_corr_growth < growth_lim, ]

# Split data into two cross-validation groups
sing_sp$cv_group <- "A"
sing_sp$cv_group[which((sing_sp$y_coord > 50 & sing_sp$x_coord <= 50) |
                         (sing_sp$y_coord < 50 & sing_sp$x_coord >= 50))] <- "B"

# Subset data into two cross-validation groups
a <- sing_sp[sing_sp$cv_group == "A", ]
b <- sing_sp[sing_sp$cv_group == "B", ]

# Create matrices of factor variables
dm_fac_a <- model.matrix(size_corr_growth ~ sps_comp, # + stand_id,
                         a, contrasts.arg = list(sps_comp = contrasts(a$sps_comp, contrasts = F)))#,
                                              #stand_id = contrasts(a$stand_id, contrasts = F)))
dm_fac_b <- model.matrix(size_corr_growth ~ sps_comp, # + stand_id,
                         b, contrasts.arg = list(sps_comp = contrasts(b$sps_comp, contrasts = F)))#,
                                              #stand_id = contrasts(b$stand_id, contrasts = F)))

# Combine factor predictors with continous predictors (prox, size_comp_dbh, all_density, species-specific densities, precip, temp)
dm_a <- cbind(dm_fac_a, as.matrix(a[c(9, 12, 14:31, 35:36)]))
dm_b <- cbind(dm_fac_b, as.matrix(b[c(9, 12, 14:31, 35:36)]))

# Standardize variables except for first column (intercept)
dm_a[, 2:ncol(dm_a)] <- apply(dm_a[, 2:ncol(dm_a)], 2, z_trans)
dm_b[, 2:ncol(dm_b)] <- apply(dm_b[, 2:ncol(dm_b)], 2, z_trans)

# Change columns of NaNs (no variation) to zeros
dm_a[, which(is.nan(dm_a[1, ]))] <- 0
dm_b[, which(is.nan(dm_b[1, ]))] <- 0

# Fit models
mod_a <- cv.glmnet(dm_a, y = a$size_corr_growth, family = "gaussian")
mod_b <- cv.glmnet(dm_b, y = b$size_corr_growth, family = "gaussian")

# Make predictions
int_pred <- c(predict(mod_a, newx = dm_b, s = "lambda.1se"),
              predict(mod_b, newx = dm_a, s = "lambda.1se"))
int_obs <- rbind(b %>% select(tree_id, size_corr_growth),
                 a %>% select(tree_id, size_corr_growth))
comparison <- cbind(int_obs, int_pred)
colnames(comparison)[2:3] <- c("obs", "pred")
comparison <- comparison %>%
  group_by(tree_id) %>%
  summarize(observations = obs[1],
            predictions = mean(pred)) # IS MEAN APPROPRIATE HERE?

# PSME results exploration
#high_growth <- int_obs$tree_id[which(int_pred > 0.075)]
#big_growers <- sing_sp[sing_sp$tree_id %in% high_growth, ]
#table(big_growers$stand_id) # Grows a lot in PP17!

# Extract appropriate axis limits for plot
axis_max <- max(c(max(comparison$predictions), max(comparison$observations))) + 0.01

# Plot predictions vs. observations
ggplot(comparison, aes(x = predictions, y = observations)) +
  geom_hex() +
  theme_bw() +
  ylim(-0.01, axis_max) +
  xlim(-0.01, axis_max) +
  geom_abline(intercept = 0, slope = 1)

# Return coefficient of determination
coef_det(comparison) # psme: 0.61, tshe: 0.28, tsme: -0.24 to -0.27, thpl: 0.03, abam: 0.20, cano: -0.01
# outliers removed; psme: 0.57, tshe: 0.27, tsme: negative still, thpl: 0.06-0.08, abam: 0.19, cano: 0.01
# stand replaced with precip and temp, outliers retained; psme: 0.52, tshe: 0.25, tsme: neg, thpl: 0-0.01, abam: 0.09, cano: 0.04
# best nbhd size: psme 0.60, tshe 0.27, tsme -0.13, thpl 0.16, abam 0.06, cano 0.09

# Calculate slope of observed growth vs. predicted growth
slope_fit <- lm(observations ~ 0 + predictions, comparison)
coef(slope_fit) # TSME: 0.95 meaning predictions increase faster than observations

# View model coefficients
mod_coef <- as.matrix(coef(mod_a))
mod_coef <- cbind(mod_coef, as.matrix(coef(mod_b)))
colnames(mod_coef) <- c("cv_a", "cv_b")

#################################
# Preformance on data used to fit
#################################

# Create matrices of factor variables
dm_fac <- model.matrix(size_corr_growth ~ sps_comp, #+ stand_id,
                       sing_sp, contrasts.arg = list(sps_comp = contrasts(sing_sp$sps_comp, contrasts = F)))#,
                                              #stand_id = contrasts(sing_sp$stand_id, contrasts = F)))

# Combine factor predictors with continous predictors (prox, size_comp_dbh, all_density, species-specific densities)
dm <- cbind(dm_fac, as.matrix(sing_sp[c(9, 12, 14:31, 35:36)]))

# Standardize variables except for first column (intercept)
dm[, 2:ncol(dm)] <- apply(dm[, 2:ncol(dm)], 2, z_trans)

# Change columns of NaNs (no variation) to zeros
dm[, which(is.nan(dm[1, ]))] <- 0

# Fit model
mod <- cv.glmnet(dm, y = sing_sp$size_corr_growth, family = "gaussian")

# Make predictions
int_pred <- predict(mod, newx = dm, s = "lambda.1se")
int_obs <- sing_sp %>% select(tree_id, size_corr_growth)
comparison <- cbind(int_obs, int_pred)
colnames(comparison)[2:3] <- c("obs", "pred")
comparison <- comparison %>%
  group_by(tree_id) %>%
  summarize(observations = obs[1],
            predictions = mean(pred)) # IS MEAN APPROPRIATE HERE?

# Extract appropriate axis limits for plot
axis_max <- max(c(max(comparison$predictions), max(comparison$observations))) + 0.01

# Plot predictions vs. observations
ggplot(comparison, aes(x = predictions, y = observations)) +
  geom_hex() +
  theme_bw() +
  ylim(-0.01, axis_max) +
  xlim(-0.01, axis_max) +
  geom_abline(intercept = 0, slope = 1)

# Return coefficient of determination
coef_det(comparison) # psme: 0.63, tshe: 0.37, tsme: 0.07-0.09, thpl: 0.17, abam: 0.27, cano: 0.23-0.26
# outliers removed; psme: 0.62, tshe: 0.38, tsme: 0.07-0.08, thpl: 0.15, abam: 0.28, cano: 0.12-0.13
# stand replaced with precip and temp, outliers retained; psme: 0.58, tshe: 0.31, tsme: 0.07, thpl: 0.13, abam: 0.15, cano: 0.25
# best nbhd size: psme 0.63, tshe 0.34, tsme 0.18, thpl 0.26, abam 0.10, cano 0.20

# Calculate slope of observed growth vs. predicted growth
slope_fit <- lm(observations ~ 0 + predictions, comparison)
coef(slope_fit) # TSME: 0.95 meaning predictions increase faster than observations

# View model coefficients
mod_coef <- cbind(mod_coef, as.matrix(coef(mod)))
colnames(mod_coef)[3] <- "all_dat"
mod_coef[which(mod_coef == 0)] <- NA
