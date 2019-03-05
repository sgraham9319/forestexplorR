devtools::load_all()
library(glmnet)


# Load mapping data
mapping <- read.csv("../Data/Mapping_2017.csv", stringsAsFactors = F)

# Create interaction matrix
int_mat <- graph_mat_all(mapping, 10)

# Remove trees whose neighborhood extends beyond stand edge
int_mat <- int_mat %>%
  filter(x_coord >= 10 & x_coord <= 90 & y_coord >= 10 & y_coord <= 90)

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

# Calculate annual growth
ann_growth <- growth_summary(growth_data)

# Add growth data to interaction data
int_mat <- left_join(int_mat, ann_growth[, c(1, ncol(ann_growth))])

# Exclude missing growth and density values
int_mat <- int_mat %>%
  filter(!is.na(size_corr_growth) &
           !is.na(all_density))

# Change focal species and stand_id columns to factors
int_mat$species <- as.factor(int_mat$species)
int_mat$stand_id <- as.factor(int_mat$stand_id)

# Divide int_mat data into two cross-validation groups
int_mat$cv_group <- "A"
int_mat$cv_group[which((int_mat$y_coord > 50 & int_mat$x_coord <= 50) |
  (int_mat$y_coord < 50 & int_mat$x_coord >= 50))] <- "B"

# Subset data into two cross-validation groups
a <- int_mat[int_mat$cv_group == "A", ]
b <- int_mat[int_mat$cv_group == "B", ]
  
# Create matrices of factor variables
dm_fac_a <- model.matrix(size_corr_growth ~ species + sps_comp + stand_id, a,
                         contrasts.arg = list(species = contrasts(a$species, contrasts = F),
                                              sps_comp = contrasts(a$sps_comp, contrasts = F),
                                              stand_id = contrasts(a$stand_id, contrasts = F)))
dm_fac_b <- model.matrix(size_corr_growth ~ species + sps_comp + stand_id, b,
                         contrasts.arg = list(species = contrasts(b$species, contrasts = F),
                                              sps_comp = contrasts(b$sps_comp, contrasts = F),
                                              stand_id = contrasts(b$stand_id, contrasts = F)))

# Combine factor predictors with continous predictors
dm_a <- cbind(dm_fac_a, a$prox, a$size_comp,
              a$all_density, a$ABAM_density,
              a$TSHE_density, a$PSME_density)
colnames(dm_a)[(ncol(dm_fac_a) + 1):ncol(dm_a)] <- c("proximity", "comp_size",
                                                     "all_density", "ABAM_density",
                                                     "TSHE_density", "PSME_density")
dm_b <- cbind(dm_fac_b, b$prox, b$size_comp,
              b$all_density, b$ABAM_density,
              b$TSHE_density, b$PSME_density)
colnames(dm_b)[(ncol(dm_fac_b) + 1):ncol(dm_b)] <- c("proximity", "comp_size",
                                                     "all_density", "ABAM_density",
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
coef_det(comparison) # 0.257
