# Load package
devtools::load_all()

#==============
# Cleaning data
#==============

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

# Load mapping data for individual trees
mapping <- read.csv("Data/Mapping_2017.csv", stringsAsFactors = F)

# Exclude small trees from mapping data
mapping <- mapping[-which(mapping$tree_id %in% small_trees), ]

# Summarize growth for all trees
growth_summ <- growth_summary(growth_data)

# Calculate neighborhood density for all trees
densities <- density_all_stands(mapping, radius = 10)

# Attach neighborhood information to growth
growth <- left_join(growth_summ, densities)

# Exclude missing growth and density values
growth <- growth %>%
  filter(!is.na(size_corr_growth) &
           !is.na(all_density))

# Change stand_id and species to factors
growth$stand_id <- as.factor(growth$stand_id)
growth$species <- as.factor(growth$species)

#====================
# Metaparameter sweep
#====================

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
