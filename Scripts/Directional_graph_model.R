library(glmnet)

# For each tree need one row for each of it's neighbors. On each row I need the
# growth and species identity of the focal tree (will be the same for each row
# referring to that tree), species identity, size, and proximity of the neighbor
# tree that row refers to, and the species-specific densities of the neighborhood
# (also the same across all rows referring to the same focal tree)

# Before creating dm_factors, need to get a dataset that has same number of rows
# as matrix described above, but only one column for focal species and one column
# for neighbor species.

# GET THIS WORKING FOR A SINGLE STAND FIRST!


# Load mapping data
mapping <- read.csv("../Data/Mapping_2017.csv", stringsAsFactors = F)

# Create interaction matrix
int_mat <- graph_mat_all(mapping, 10)

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

# Create model matrix of factorized predictor variables
dm_factors <- model.matrix(size_corr_growth ~ species + sps_comp + stand_id, int_mat,
                           contrasts.arg = list(species = contrasts(int_mat$species, contrasts = F),
                                                sps_comp = contrasts(int_mat$sps_comp, contrasts = F),
                                                stand_id = contrasts(int_mat$stand_id, contrasts = F)))

# Combine factorized predictors with continuous predictors
dm <- cbind(dm_factors, int_mat$prox, int_mat$size_comp,
            int_mat$all_density, int_mat$ABAM_density,
            int_mat$TSHE_density, int_mat$PSME_density)
colnames(dm)[(ncol(dm_factors) + 1):ncol(dm)] <- c("proximity", "comp_size",
                                                   "all_density", "ABAM_density",
                                                   "TSHE_density", "PSME_density")

# Remove factor levels always equal to zero
dm <- dm[, -c(which(apply(dm, 2, max) == 0))]

# Standardize variables
dm[,2:ncol(dm)] <- apply(dm[,2:ncol(dm)], 2, z_trans)

# Fit model
glmmod <- cv.glmnet(dm, y = int_mat$size_corr_growth, family = "gaussian")

# Plot model
plot(glmmod)

# Look at coefficients
coef(glmmod, s = "lambda.min")
coef(glmmod, s = "lambda.1se")


#========================
# Fitting model to 1 tree
#========================

# rows 1 through 7
dm_few_tree <- dm[1:90, ]
mod <- cv.glmnet(dm_few_tree, y = int_mat$size_corr_growth[1:90], family = "gaussian")
predictions <- predict(mod, newx = dm_few_tree, s = "lambda.1se")
observations <- int_mat$size_corr_growth[1:90]
comparison <- data.frame(predictions, observations)
ggplot(comparison, aes(x = observations, y = X1)) +
  geom_hex() +
  theme_bw() +
  geom_abline(intercept = 0, slope = 1)
length(unique(predictions))

coef(mod, s = "lambda.1se")
# 47 through 57

#==================================================
# Compare observations to predicitons in full model
#==================================================

predictions <- predict(glmmod, newx = dm, s = "lambda.1se")
observations <- int_mat$size_corr_growth
comparison <- data.frame(predictions, observations)
ggplot(comparison, aes(x = observations, y = X1)) +
  geom_hex() +
  theme_bw() +
  ylim(0, 0.3) +
  geom_abline(intercept = 0, slope = 1)
plot(glmmod)
coef_det(comparison) # 0.294, 0.298, 0.298, 0.296, 0.297
coef(glmmod, s = "lambda.1se")
# 5 competitor species identities retained and size of competitor retained

predict_eval <- function(mod, dm, y){
  predictions <- predict(mod, newx = dm, s = "lambda.1se")
  observations <- y
  comparison <- data.frame(predictions, observations)
  print(coef_det(comparison))
  ggplot(comparison, aes(x = observations, y = X1)) +
    geom_hex() +
    theme_bw() +
    ylim(0, 0.3) +
    geom_abline(intercept = 0, slope = 1)
}
predict_eval(glmmod, dm, int_mat$size_corr_growth)
