#########################################
# Single-species directional graph models
#########################################

# Author: Stuart Graham
# Created: 4/27/2020
# Last edited: 4/27/2020

# Load TreeNeighborhood package
devtools::load_all()

# Load other required packages
library(glmnet)
library(ggplot2)

######################
# Part 1. Loading data
######################

# Load mapping data
mapping <- read.csv("../TreeNeighborhood/Data/Cleaned_mapping_2017.csv", stringsAsFactors = F)

# Load tree measurement data
tree <- read.csv("../TreeNeighborhood/Data/Cleaned_tree_growth_2017.csv", stringsAsFactors = F)

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

# Obtain all neighborhood data
neighbors <- graph_mat_all(mapping, radius = 10)

# Remove focals whose neighborhood overlaps stand boundary
neighbors <- neighbors %>%
  filter(x_coord >= 10 & x_coord <= 90 & y_coord >= 10 & y_coord <= 90)

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

################################################
# Part 5. Combining growth and neighborhood data
################################################

# Extract required columns from growth data
growth_cols <- growth[, c("tree_id", "midpoint_size", "radial_growth", "size_corr_growth")]

# Join growth and neighborhood data. Use inner join because there will be
# no growth data for focals measured only once or with negative growth
full <- inner_join(neighbors, growth_cols, by = "tree_id")

# Change stand_id column to factor
full$stand_id <- as.factor(full$stand_id)

#################################
# Part 6. Creating model matrices
#################################

# Subset to a single species
sing_sp <- droplevels(full[full$species == "TSHE", ])

# Split data into two cross-validation groups
sing_sp$cv_group <- "A"
sing_sp$cv_group[which((sing_sp$y_coord > 50 & sing_sp$x_coord <= 50) |
                         (sing_sp$y_coord < 50 & sing_sp$x_coord >= 50))] <- "B"

# Subset data into two cross-validation groups
a <- sing_sp[sing_sp$cv_group == "A", ]
b <- sing_sp[sing_sp$cv_group == "B", ]

# Create matrices of factor variables
dm_fac_a <- model.matrix(size_corr_growth ~ sps_comp + stand_id, a,
                         contrasts.arg = list(sps_comp = contrasts(a$sps_comp, contrasts = F),
                                              stand_id = contrasts(a$stand_id, contrasts = F)))
dm_fac_b <- model.matrix(size_corr_growth ~ sps_comp + stand_id, b,
                         contrasts.arg = list(sps_comp = contrasts(b$sps_comp, contrasts = F),
                                              stand_id = contrasts(b$stand_id, contrasts = F)))

# Combine factor predictors with continous predictors
#dm_a <- cbind(dm_fac_a, as.matrix(a[c(7, 9, 10, 11, 21, 24)]))
#dm_b <- cbind(dm_fac_b, as.matrix(b[c(7, 9, 10, 11, 21, 24)]))
# Try excluding species-specific densities to see if more variation is attributed to
# species-specific interaction coefficients
dm_a <- cbind(dm_fac_a, as.matrix(a[c(7, 9, 10)]))
dm_b <- cbind(dm_fac_b, as.matrix(b[c(7, 9, 10)]))

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
            predictions = mean(pred))

# Extract appropriate axis limits for plot
axis_max <- max(c(max(comparison$predictions), max(comparison$observations)))

# Plot predictions vs. observations
ggplot(comparison, aes(x = observations, y = predictions)) +
  geom_hex() +
  theme_bw() +
  ylim(0, axis_max) +
  xlim(0, axis_max) +
  geom_abline(intercept = 0, slope = 1)

# Return coefficient of determination
coef_det(comparison) # tsme: -0.146, psme: 0.62 (0.61), abam: 0.27, cano: -0.007, tshe: 0.31 (0.25)
coef(mod_b)
