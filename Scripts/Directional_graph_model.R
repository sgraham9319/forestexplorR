devtools::load_all()
library(glmnet)
library(ggplot2)
library(tidyr)

#===========================
# Loading and formating data
#===========================

# Load mapping data
mapping <- read.csv("Data/Mapping_2017.csv", stringsAsFactors = F)

# Create interaction matrix
int_mat <- graph_mat_all(mapping, radius = 10)

# Remove trees whose neighborhood extends beyond stand edge
int_mat <- int_mat %>%
  filter(x_coord >= 10 & x_coord <= 90 & y_coord >= 10 & y_coord <= 90)

# Load growth data
growth_data <- read.csv("Data/Tree_growth_2017.csv", stringsAsFactors = F)

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

#==================
# Create full model
#==================

# Create matrices of factor variables
dm_fac <- model.matrix(size_corr_growth ~ species + sps_comp + stand_id, int_mat,
                         contrasts.arg = list(species = contrasts(int_mat$species, contrasts = F),
                                              sps_comp = contrasts(int_mat$sps_comp, contrasts = F),
                                              stand_id = contrasts(int_mat$stand_id, contrasts = F)))

# Combine factor predictors with continous predictors
dm <- cbind(dm_fac, as.matrix(int_mat[c(7, 9, 10, 11, 24, 27)]))

# Standardize variables except for first column (intercept)
dm[, 2:ncol(dm)] <- apply(dm[, 2:ncol(dm)], 2, z_trans)

# Change columns of NaNs (no variation) to zeros
dm[, which(is.nan(dm[1, ]))] <- 0

# Fit models
mod <- cv.glmnet(dm, y = int_mat$size_corr_growth, family = "gaussian")

# Make predictions
pred <- predict(mod, newx = dm, s = "lambda.1se")
obs <- int_mat %>% select(tree_id, species, sps_comp, size_corr_growth)
all_pred <- cbind(obs, pred)
colnames(all_pred)[4:5] <- c("obs", "pred")
plot(mod)
coef(mod)
#=====================
# Cross-validate model
#=====================

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
dm_a <- cbind(dm_fac_a, as.matrix(a[c(7, 9, 10, 11, 24, 27)]))
dm_b <- cbind(dm_fac_b, as.matrix(b[c(7, 9, 10, 11, 24, 27)]))

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

#=================================
# Quantifying species interactions
#=================================

# Get species list
all_sps <- unique(c(as.character(all_pred$species),
                    as.character(all_pred$sps_comp)))

# Create output matrix
result <- matrix(NA, nrow = length(all_sps), ncol = length(all_sps))
rownames(result) <- all_sps
colnames(result) <- all_sps

# Create se and sample size output matrices
se_mat <- result
samp_size <- result

for(focal in 1:length(all_sps)){
  
  # Obtain single prediction for each focal by averaging over interactions
  complete_nbhd <- all_pred %>%
    filter(species == all_sps[focal]) %>%
    group_by(tree_id) %>%
    summarize(pred = mean(pred))
  
  for(comp in 1:length(all_sps)){
    
    # Calculate average predictions with one competitor species excluded
    reduced_nbhd <- all_pred %>%
      filter(species == all_sps[focal] & sps_comp != all_sps[comp]) %>%
      group_by(tree_id) %>%
      summarize(pred_sub = mean(pred))
    
    # Match and drop those with no prediction (those that only had this competitor
    # species as a competitor)
    pred_comparison <- left_join(reduced_nbhd, complete_nbhd)
    
    # Obtain average effect of this competitor species on this focal species
    result[focal, comp] <- mean(pred_comparison$pred - pred_comparison$pred_sub)
    se_mat[focal, comp] <- st_err(pred_comparison$pred - pred_comparison$pred_sub)
    samp_size[focal, comp] <- sum(pred_comparison$pred != pred_comparison$pred_sub)
  }
}

# Each cell is the effect of a neighbor of species (column name) on the growth
# of species (row name)
gathered_result <- gather(data.frame(result), key = "comp_sps",
                          value = "growth_effect") %>%
  mutate(focal_sps = rep(rownames(result), times = ncol(result)))
gathered_result <- gathered_result[, c(3, 1, 2)]

# Remove rare species (appear less than 100 times as competitor or focal in
# interaction matrix)
rare_sps <- unique(c(names(which(table(int_mat$species) < 100)),
                     names(which(table(int_mat$sps_comp) < 100))))
gathered_result <- gathered_result %>%
  filter(!(focal_sps %in% rare_sps) &
           !(comp_sps %in% rare_sps))

# Create heat map
ggplot(data = gathered_result, aes(x = comp_sps, y = focal_sps,
                                   fill = growth_effect)) +
  labs(x = "Competitor species", y = "Focal species") +
  geom_tile() +
  coord_equal() +
  scale_fill_gradient2(breaks = c(min(gathered_result$growth_effect), 0, 
                                  max(gathered_result$growth_effect)),
                       labels = c("Negative", "No effect", "Positive"), 
                       low = "purple", high = "green", mid = "black",
                       midpoint = 0, name = "Effect on annual growth") +
  theme_minimal()


# Removing those interactions whose standard error overlaps with zero
result2 <- result
result2[which((abs(result2) - se_mat) <= 0)] <- NA
gathered_result2 <- gather(data.frame(result2), key = "comp_sps",
                          value = "growth_effect") %>%
  mutate(focal_sps = rep(rownames(result2), times = ncol(result2)))
gathered_result2 <- gathered_result2[, c(3, 1, 2)]
rare_sps <- unique(c(names(which(table(int_mat$species) < 100)),
                     names(which(table(int_mat$sps_comp) < 100))))
gathered_result2 <- gathered_result2 %>%
  filter(!(focal_sps %in% rare_sps) &
           !(comp_sps %in% rare_sps))
ggplot(data = gathered_result2, aes(x = comp_sps, y = focal_sps,
                                   fill = growth_effect)) +
  labs(x = "Competitor species", y = "Focal species") +
  geom_tile() +
  coord_equal() +
  scale_fill_gradient2(breaks = c(min(gathered_result2$growth_effect), 0, 
                                  max(gathered_result2$growth_effect)),
                       labels = c("Negative", "No effect", "Positive"), 
                       low = "purple", high = "green", mid = "black",
                       midpoint = 0, name = "Effect on annual growth") +
  theme_minimal()

#===============================================================
# Comparing no competitor growth to growth with focal competitor
#===============================================================

# Create output matrix
result1 <- matrix(NA, nrow = length(all_sps), ncol = length(all_sps))
rownames(result1) <- all_sps
colnames(result1) <- all_sps

# Create se and sample size output matrices
se_mat1 <- result1
samp_size1 <- result1

for(focal in 1:length(all_sps)){
  
  # Obtain single prediction for each focal by averaging over interactions
  complete_nbhd <- all_pred %>%
    filter(species == all_sps[focal]) %>%
    group_by(tree_id) %>%
    summarize(pred = mean(pred))
  
  for(comp in 1:length(all_sps)){
    
    # Calculate average predictions with one competitor species excluded
    reduced_nbhd <- all_pred %>%
      filter(species == all_sps[focal] & sps_comp != all_sps[comp]) %>%
      group_by(tree_id) %>%
      summarize(pred_sub = mean(pred))
    
    # Match and drop those with no prediction (those that only had this competitor
    # species as a competitor)
    pred_comparison <- left_join(reduced_nbhd, complete_nbhd)
    
    # Obtain average effect of this competitor species on this focal species
    result[focal, comp] <- mean(pred_comparison$pred - pred_comparison$pred_sub)
    se_mat[focal, comp] <- st_err(pred_comparison$pred - pred_comparison$pred_sub)
    samp_size[focal, comp] <- sum(pred_comparison$pred != pred_comparison$pred_sub)
  }
}

# Each cell is the effect of a neighbor of species (column name) on the growth
# of species (row name)
gathered_result <- gather(data.frame(result), key = "comp_sps",
                          value = "growth_effect") %>%
  mutate(focal_sps = rep(rownames(result), times = ncol(result)))
gathered_result <- gathered_result[, c(3, 1, 2)]

# Remove rare species (appear less than 100 times as competitor or focal in
# interaction matrix)
rare_sps <- unique(c(names(which(table(int_mat$species) < 100)),
                     names(which(table(int_mat$sps_comp) < 100))))
gathered_result <- gathered_result %>%
  filter(!(focal_sps %in% rare_sps) &
           !(comp_sps %in% rare_sps))

# Create heat map
ggplot(data = gathered_result, aes(x = comp_sps, y = focal_sps,
                                   fill = growth_effect)) +
  labs(x = "Competitor species", y = "Focal species") +
  geom_tile() +
  scale_fill_gradient2(breaks = c(min(gathered_result$growth_effect), 0, 
                                  max(gathered_result$growth_effect)),
                       labels = c("Negative", "No effect", "Positive"), 
                       low = "purple", high = "green", mid = "black",
                       midpoint = 0, name = "Effect on annual growth") +
  theme_minimal()