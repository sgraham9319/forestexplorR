
# Load mapping and growth data
mapping <- read.csv("../Data/Mapping_2017.csv", stringsAsFactors = F)
growth <- read.csv("../Data/Tree_growth_2017.csv", stringsAsFactors = F)

# Remove test data (2017 measurements and stands TO04, AE10, and AV02)
growth <- growth %>%
  filter(year != 2017,
         stand_id != "TO04",
         stand_id != "AE10",
         stand_id != "AV02")

# Remove rows with no dbh measurement
growth <- growth %>%
  filter(!is.na(dbh))

# Identify trees larger than minimum size threshold
trees <- growth %>%
  group_by(tree_id) %>% 
  summarize(max_size = max(dbh)) %>%
  filter(max_size >= 15) %>%
  select(tree_id)

# Exclude small trees
growth <- growth %>%
  filter(tree_id %in% trees$tree_id)

# Summarize growth for all trees
growth_summ <- growth_summary(growth)

# Calculate neighborhood density for all trees
densities <- density_all_stands(mapping, radius = 10)

# Attach neighborhood information to growth
no_int_mat <- left_join(growth_summ, densities)

# Exclude missing growth and density values and zero density values
no_int_mat <- no_int_mat %>%
  filter(!is.na(size_corr_growth),
         !is.na(all_density),
         all_density > 0)

# Remove trees whose neighborhood extends beyond stand edge
no_int_mat <- no_int_mat %>%
  filter(x_coord >= 10 & x_coord <= 90 & y_coord >= 10 & y_coord <= 90)

# Change stand_id and species to factors
no_int_mat$stand_id <- as.factor(no_int_mat$stand_id)
no_int_mat$species <- as.factor(no_int_mat$species)

# Divide int_mat data into two cross-validation groups
no_int_mat$cv_group <- "A"
no_int_mat$cv_group[which((no_int_mat$y_coord > 50 & no_int_mat$x_coord <= 50) |
                         (no_int_mat$y_coord < 50 & no_int_mat$x_coord >= 50))] <- "B"

#====================================================
# Model 1: no species-specific competitor information
#====================================================

# Subset data into two cross-validation groups
a <- no_int_mat[no_int_mat$cv_group == "A", ]
b <- no_int_mat[no_int_mat$cv_group == "B", ]

# Create matrices of factor variables
dm_fac_a <- model.matrix(size_corr_growth ~ species + stand_id, a,
                         contrasts.arg = list(species = contrasts(a$species, contrasts = F),
                                              stand_id = contrasts(a$stand_id, contrasts = F)))
dm_fac_b <- model.matrix(size_corr_growth ~ species + stand_id, b,
                         contrasts.arg = list(species = contrasts(b$species, contrasts = F),
                                              stand_id = contrasts(b$stand_id, contrasts = F)))

# Combine factor predictors with continous predictors
dm_a <- cbind(dm_fac_a, a$all_density)
colnames(dm_a)[(ncol(dm_fac_a) + 1):ncol(dm_a)] <- "all_density"
dm_b <- cbind(dm_fac_b, b$all_density)
colnames(dm_b)[(ncol(dm_fac_b) + 1):ncol(dm_b)] <- "all_density"

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
predictions <- c(predict(mod_a, newx = dm_b, s = "lambda.1se"),
                 predict(mod_b, newx = dm_a, s = "lambda.1se"))
observations <- c(b$size_corr_growth, a$size_corr_growth)
comparison <- data.frame(predictions, observations)

# Extract appropriate axis limits for plot
axis_max <- max(c(max(predictions), max(observations)))

# Plot predictions vs. observations
ggplot(comparison, aes(x = observations, y = predictions)) +
  geom_hex() +
  theme_bw() +
  ylim(0, axis_max) +
  xlim(0, axis_max) +
  geom_abline(intercept = 0, slope = 1)

# Return coefficient of determination
coef_det(comparison) # 0.231

#=============================================
# Model 2: species-specific densities included
#=============================================

# Combine factor predictors with continous predictors
dm_a <- cbind(dm_fac_a, a$all_density, a$ABAM_density,
              a$TSHE_density, a$PSME_density)
colnames(dm_a)[(ncol(dm_fac_a) + 1):ncol(dm_a)] <- c("all_density", "ABAM_density",
                                                     "TSHE_density", "PSME_density")
dm_b <- cbind(dm_fac_b, b$all_density, b$ABAM_density,
              b$TSHE_density, b$PSME_density)
colnames(dm_b)[(ncol(dm_fac_b) + 1):ncol(dm_b)] <- c("all_density", "ABAM_density",
                                                     "TSHE_density", "PSME_density")

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
predictions <- c(predict(mod_a, newx = dm_b, s = "lambda.1se"),
                 predict(mod_b, newx = dm_a, s = "lambda.1se"))
observations <- c(b$size_corr_growth, a$size_corr_growth)
comparison <- data.frame(predictions, observations)

# Extract appropriate axis limits for plot
axis_max <- max(c(max(predictions), max(observations)))

# Plot predictions vs. observations
ggplot(comparison, aes(x = observations, y = predictions)) +
  geom_hex() +
  theme_bw() +
  ylim(0, axis_max) +
  xlim(0, axis_max) +
  geom_abline(intercept = 0, slope = 1)

# Return coefficient of determination
coef_det(comparison) # 0.239
