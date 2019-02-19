# Load package
devtools::load_all()

library(moments)
library(glmnet)
library(glmnetUtils)

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
growth <- growth_summary(growth_data)

# Calculate neighborhood density for all trees
densities <- density_all_stands(mapping)

# Attach neighborhood information to growth
growth <- left_join(growth, densities)


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
dm <- as.matrix(data.frame(growth1$all_density, dm_factors))

# Fit model
glmmod <- cv.glmnet(dm, y = growth1$size_corr_growth, family = "gaussian")

# Plot model object - lower mean-squared error means better fit and the numbers
# on the top show the number of variables included. Region delimited by dotted
# vertical lines represents equally good model fit
plot(glmmod)

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

levels(growth1$stand_id)

growth1 %>%
  group_by(stand_id) %>%
  summarize(num_tshe = sum(species == "TSHE"),
            num_abam = sum(species == "ABAM"),
            num_psme = sum(species == "PSME"))
