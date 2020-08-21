############################
# Checking for overfitting #
############################

# This script checks for overfitting in the graphical model under three
# different methods for defining the training and test data. This helps
# us determine the best method for splitting the data into training and
# test because we want to use as much of the data as possible without
# causing any overfitting

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

# Load environmental data
env_dat <- read.csv("Data/stand_abiotic_data.csv", stringsAsFactors = F)

################################
# Part 2. Creating neighborhoods
################################

# Obtain all neighborhood data - different radius for each species
neighbors <- graph_mat_all(mapping, radius = 20)

# Remove focals whose neighborhood overlaps stand boundary
neighbors <- neighbors %>% filter(x_coord >= 20 & x_coord <= 80 & y_coord >= 20 & y_coord <= 80)

# Remove small competitors
neighbors <- neighbors %>% filter(size_cat_comp == "regular")

###################################
# Part 3. Calculating annual growth
###################################

# Calculate annual growth for all trees
growth <- growth_summary(tree)

# Remove trees for which annual growth could not be calculated (see likelihood
# model script for details)
growth <- growth[!is.na(growth$size_corr_growth), ]

# Check data one more time
sum(is.nan(growth$size_corr_growth))

###############################################################
# Part 4. Combining growth, neighborhood and environmental data
###############################################################

# Extract required columns from growth data
growth_cols <- growth[, c("tree_id", "midpoint_size", "size_corr_growth")]

# Join growth and neighborhood data. Use inner join because there will be
# no growth data for focals measured only once or with negative growth
full <- inner_join(neighbors, growth_cols, by = "tree_id")

# Add environmental data to neighborhoods and growth
full <- left_join(full, env_dat)

######################################
# Part 5. Formatting data for modeling
######################################

# Exclude test data
full_sub1 <- full %>% filter(y_coord <= 40)
full_sub2 <- full %>% filter(x_coord >= 60 & y_coord > 40)
no_test <- rbind(full_sub1, full_sub2)

# Subset to a single species
sing_sp <- droplevels(no_test[no_test$species == "TSHE", ])

# Find common competitors
comps <- names(which(table(sing_sp$sps_comp) >= 100))

# Change rare competitor species to OTHR
sing_sp$sps_comp[which(sing_sp$sps_comp %in% comps == F)] <- "OTHR"

# Define densities to include
comps <- c(comps, "all")
density_cols <- rep(NA, times = length(comps))
for(i in 1:length(comps)){
  density_cols[i] <- grep(comps[i], names(sing_sp))
}
other_cols <- setdiff(grep("density", names(sing_sp)), density_cols)
sing_sp$OTHR_density <- apply(sing_sp[, other_cols], 1, sum)

# Divide into training and validation data sets
#validation <- sing_sp %>% filter(x_coord < 30)
#training <- sing_sp %>% filter(x_coord >= 70 | y_coord > 67)
validation <- sing_sp %>% filter(x_coord < 40)
training <- sing_sp %>% filter(x_coord >= 60)
#validation <- sing_sp %>% filter(x_coord < 50)
#training <- sing_sp %>% filter(x_coord >= 50)

# Randomize growth measurements among focals in training set
focals <- training %>% group_by(tree_id) %>% summarize(size_corr_growth = size_corr_growth[1])
new_order <- sample(1:nrow(focals), nrow(focals))
focals$size_corr_growth_rand <- focals$size_corr_growth[new_order]
training <- inner_join(training, focals[, c(1, 3)], by = "tree_id")

# Convert competitor species to factor
training$sps_comp <- as.factor(training$sps_comp)
validation$sps_comp <- as.factor(validation$sps_comp)

#####################################
# Part 6. Building tree growth models
#####################################

# Create matrices of factor variables
#dm_fac <- model.matrix(size_corr_growth_rand ~ sps_comp, # randomized growth
#                      training, contrasts.arg = list(sps_comp = contrasts(training$sps_comp, contrasts = F)))
dm_fac <- model.matrix(size_corr_growth ~ sps_comp, # true growth
                       training, contrasts.arg = list(sps_comp = contrasts(training$sps_comp, contrasts = F)))

# Combine factor predictors with continuous predictors (prox, size_comp_dbh, all_density, species-specific densities)
dm <- cbind(dm_fac, as.matrix(training[c(9, 12, density_cols, 35:36, 38)]))

# Standardize variables except for first column (intercept)
dm[, 2:ncol(dm)] <- apply(dm[, 2:ncol(dm)], 2, z_trans)

# Change columns of NaNs (no variation) to zeros
dm[, which(is.nan(dm[1, ]))] <- 0

# Fit model
#mod <- cv.glmnet(dm, y = training$size_corr_growth_rand, family = "gaussian") # randomized growth
mod <- cv.glmnet(dm, y = training$size_corr_growth, family = "gaussian") # true growth

# Create model matrix of validation data
dm_fac_v <- model.matrix(size_corr_growth ~ sps_comp,
                         validation, contrasts.arg = list(sps_comp = contrasts(validation$sps_comp, contrasts = F)))
dm_v <- cbind(dm_fac_v, as.matrix(validation[c(9, 12, density_cols, 35:36, 38)]))
dm_v[, 2:ncol(dm_v)] <- apply(dm_v[, 2:ncol(dm_v)], 2, z_trans)
dm_v[, which(is.nan(dm_v[1, ]))] <- 0

##################################
# Part 7. Checking for overfitting
##################################

# Add missing columns to dm_v
#new_cols <- rep(0, times = nrow(dm_v))
#colnames(new_cols) <- c("sps_compOTHR", "sps_compPICO", "sps_compPSME")
#dm_v <- cbind(dm_v[, 1], new_cols, dm_v[, 2:ncol(dm_v)])
#colnames(dm_v)[2] <- "sps_compABAM"

# Make predictions
int_pred <- predict(mod, newx = dm_v, s = "lambda.1se")
int_obs <- validation %>% select(tree_id, size_corr_growth)
comparison <- cbind(int_obs, int_pred)
colnames(comparison)[2:3] <- c("obs", "pred")
comparison <- comparison %>%
  group_by(tree_id) %>%
  summarize(observations = obs[1],
            predictions = mean(pred))

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
coef_det(comparison)

# Calculate slope of observed growth vs. predicted growth
slope_fit <- lm(observations ~ 0 + predictions, comparison)
coef(slope_fit)
