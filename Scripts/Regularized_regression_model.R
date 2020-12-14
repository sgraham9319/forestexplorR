###############################################
# Regularized regression model of tree growth #
###############################################

# Author: Stuart Graham
# Created: 10/3/2020
# Last edited: 11/10/2020

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

################################
# Part 2. Creating neighborhoods
################################

# Obtain all neighborhood data
neighbors <- graph_mat_all(mapping, radius = 15)

# Remove focals whose neighborhood overlaps stand boundary
neighbors <- neighbors %>% filter(x_coord >= 15 & x_coord <= 85 & y_coord >= 15 & y_coord <= 85)

# Remove small competitors
neighbors <- neighbors %>% filter(size_cat_comp == "regular")

###################################
# Part 3. Calculating annual growth
###################################

# Calculate annual growth for all trees
growth <- growth_summary(tree)

# Produces some NA values. Remove those caused by a tree being measured only once,
# which would result in an annual growth measurement of NA
growth <- growth[growth$first_record != growth$last_record, ]

# Check for any remaining NA annual growth values
sum(is.na(growth$annual_growth))

# Check for remaining NaNs in size corrected growth
sum(is.nan(growth$size_corr_growth))

# NaNs could be caused by negative annual growth values. Remove these
growth <- growth[growth$annual_growth >= 0, ]

# Confirm no remaining NaNs in size corrected growth
sum(is.nan(growth$size_corr_growth))

################################################
# Part 4. Combining growth and neighborhood data
################################################

# Extract required columns from growth data
growth_cols <- growth[, c("tree_id", "midpoint_size", "size_corr_growth")]

# Join growth and neighborhood data. Use inner join because there will be
# no growth data for focals measured only once or with negative growth
full <- inner_join(neighbors, growth_cols, by = "tree_id")

# Change stand_id column to factor
full$stand_id <- as.factor(full$stand_id)

##############################################
# Part 5. Dividing into training and test data
##############################################

# Specify training data subsets
#full_sub1 <- full %>% filter(y_coord <= 50 - (15 / 2)) # Training 1
#full_sub2 <- full %>% filter(x_coord >= 50 + (15 / 2) & y_coord > 50 - (15 / 2)) # Training 1
#full_sub1 <- full %>% filter(x_coord <= 50 - (15 / 2)) # Training 2
#full_sub2 <- full %>% filter(y_coord <= 50 - (15 / 2) & x_coord > 50 - (15 / 2)) # Training 2
#full_sub1 <- full %>% filter(y_coord >= 50 + (15 / 2)) # Training 3
#full_sub2 <- full %>% filter(x_coord <= 50 - (15 / 2) & y_coord < 50 + (15 / 2)) # Training 3
full_sub1 <- full %>% filter(x_coord >= 50 + (15 / 2)) # Training 4
full_sub2 <- full %>% filter(y_coord >= 50 + (15 / 2) & x_coord < 50 + (15 / 2)) # Training 4

# Define full training set
training <- rbind(full_sub1, full_sub2)

# Define test data
#test <- full %>% filter(y_coord >= 50 + (15 / 2) & x_coord <= 50 - (15 / 2)) # Test 1
#test <- full %>% filter(y_coord >= 50 + (15 / 2) & x_coord >= 50 + (15 / 2)) # Test 2
#test <- full %>% filter(y_coord <= 50 - (15 / 2) & x_coord >= 50 + (15 / 2)) # Test 3
test <- full %>% filter(y_coord <= 50 - (15 / 2) & x_coord <= 50 - (15 / 2)) # Test 4


######################################
# Part 6. Creating model design matrix
######################################

# Subset to a single species
sing_sp <- droplevels(training[training$species == "CANO", ])

# Find common competitors
comps <- names(which(table(sing_sp$sps_comp) >= 100))

# Change rare competitor species to OTHR
sing_sp$sps_comp[which(sing_sp$sps_comp %in% comps == F)] <- "OTHR"

# Create OTHR density column by summing densities over rare competitors
comps <- c(comps, "all")
density_cols <- rep(NA, times = length(comps))
for(i in 1:length(comps)){
  density_cols[i] <- grep(comps[i], names(sing_sp))
}
other_cols <- setdiff(grep("density", names(sing_sp)), density_cols)
sing_sp$OTHR_density <- apply(sing_sp[, other_cols], 1, sum)

# Add intra vs. inter column
sing_sp$intra <- 0
sing_sp$intra[sing_sp$sps_comp == sing_sp$species] <- 1

# Convert competitor species to factor
sing_sp$sps_comp <- as.factor(sing_sp$sps_comp)

# Create matrix of factor variables
dm_fac <- model.matrix(size_corr_growth ~ sps_comp + stand_id, sing_sp,
                       contrasts.arg = list(sps_comp = contrasts(sing_sp$sps_comp, contrasts = F),
                       stand_id = contrasts(sing_sp$stand_id, contrasts = F)))

# Combine factor predictors with continuous predictors (prox, size_comp_dbh, all_density, species-specific densities)
dm <- cbind(dm_fac, as.matrix(sing_sp[c(9, 12, density_cols, 35, 36)])) # including intra

# Standardize variables except for first column (intercept)
dm[, 2:ncol(dm)] <- apply(dm[, 2:ncol(dm)], 2, z_trans)

# Change columns of NaNs (no variation) to zeros
dm[, which(is.nan(dm[1, ]))] <- 0

###########################
# Part 7. Fitting the model
###########################

# Fit the model 100 times
all_fits_list <- mod_iter(des_mat = dm, growth_dat = sing_sp, runs = 100)

# Separate coefficients table from model in output
all_fits <- all_fits_list[[1]]
mod <- all_fits_list[[2]]

# Identify parameter columns needed to test hypotheses
comp_id_cols <- grep("sps_comp", names(all_fits))
prox_col <- which(names(all_fits) == "prox")
size_comp_col <- which(names(all_fits) == "size_comp_dbh")
intra_col <- which(names(all_fits) == "intra")
dens_cols <- grep("density", names(all_fits))

# Create new summation columns needed to test hypotheses
all_fits$nbhds <- abs(apply(all_fits[, c(comp_id_cols, prox_col, size_comp_col,
                                     intra_col, dens_cols)], 1, sum)) 
all_fits$compID <- abs(apply(all_fits[, c(comp_id_cols, intra_col)], 1, sum))

# In what percentage of model runs is competitive neighborhood important?
pct_dif(all_fits$nbhds, "greater", 0)

# In what percentage of model runs is species identity of competitors important?
pct_dif(all_fits$compID, "greater", 0)

# In what percentage of model runs is conspecific competition strongest?
pct_dif(all_fits$intra, "less", 0)

# In what percentage of model runs is heterospecific competition strongest?
pct_dif(all_fits$intra, "greater", 0)

# In what percentage of models is each species a strong competitor?
apply(all_fits[, comp_id_cols], 2, pct_dif, direction = "less", value = 0) 

# In what percentage of models is each species a weak competitor?
apply(all_fits[, comp_id_cols], 2, pct_dif, direction = "greater", value = 0) 

# Make predictions of training data using best model
int_pred <- predict(mod, newx = dm, s = "lambda.1se")
int_obs <- sing_sp %>% select(tree_id, size_corr_growth)
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

# Compute coefficient of determination
coef_det(comparison)

# Calculate slope of observed growth vs. predicted growth
slope_fit <- lm(observations ~ 0 + predictions, comparison)
coef(slope_fit)

# Extract model coefficients
#mod_coef <- as.matrix(coef(mod))
#mod_coef <- cbind(mod_coef, as.matrix(coef(mod)))

# Save model coefficients table
write.csv(all_fits, "Data/CANO4_no_clim.csv")

##############################
# Part 8. Predicting test data
##############################

# Subset to a single species
test_ss <- droplevels(test[test$species == "CANO", ])

# Change rare competitor species to OTHR
test_ss$sps_comp[which(test_ss$sps_comp %in% comps == F)] <- "OTHR"

# Create OTHR density column by summing densities over rare competitors
test_ss$OTHR_density <- apply(test_ss[, other_cols], 1, sum)

# Add intra vs. inter column
test_ss$intra <- 0
test_ss$intra[test_ss$sps_comp == test_ss$species] <- 1

# Convert competitor species to factor
test_ss$sps_comp <- as.factor(test_ss$sps_comp)

# Create matrix of factor variables
dm_fac_test <- model.matrix(size_corr_growth ~ sps_comp + stand_id, test_ss,
                       contrasts.arg = list(sps_comp = contrasts(test_ss$sps_comp, contrasts = F),
                                            stand_id = contrasts(test_ss$stand_id, contrasts = F)))

# Combine factor predictors with continuous predictors (prox, size_comp_dbh, all_density, species-specific densities)
dm_test <- cbind(dm_fac_test, as.matrix(test_ss[c(9, 12, density_cols, 35, 36)])) # including intra

# Add any columns to dm_test that are missing (but present in dm)
missing_cols <- colnames(dm)[colnames(dm) %in% colnames(dm_test) == F]
for(i in 1:length(missing_cols)){
  dm_test <- cbind(dm_test, rep(0, times = nrow(dm_test)))
  colnames(dm_test)[ncol(dm_test)] <- missing_cols[i]
}
dm_test <- dm_test[, match(colnames(dm), colnames(dm_test))]

# Standardize variables except for first column (intercept)
dm_test[, 2:ncol(dm_test)] <- apply(dm_test[, 2:ncol(dm_test)], 2, z_trans)

# Change columns of NaNs (no variation) to zeros
dm_test[, which(is.nan(dm_test[1, ]))] <- 0

# Make predictions
test_pred <- predict(mod, newx = dm_test, s = "lambda.1se")
test_obs <- test_ss %>% select(tree_id, size_corr_growth)
comparison_test <- cbind(test_obs, test_pred)
colnames(comparison_test)[2:3] <- c("obs", "pred")
comparison_test <- comparison_test %>%
  group_by(tree_id) %>%
  summarize(observations = obs[1],
            predictions = mean(pred))

# Extract appropriate axis limits for plot
axis_max <- max(c(max(comparison_test$predictions), max(comparison_test$observations))) + 0.01

# Plot predictions vs. observations
ggplot(comparison_test, aes(x = predictions, y = observations)) +
  geom_hex() +
  theme_bw() +
  ylim(-0.01, axis_max) +
  xlim(-0.01, axis_max) +
  geom_abline(intercept = 0, slope = 1)

# Compute coefficient of determination
coef_det(comparison_test)

# Calculate slope of observed growth vs. predicted growth
slope_fit <- lm(observations ~ 0 + predictions, comparison_test)
coef(slope_fit)
