###########################################################
# Neighborhood analysis of MORA data using likelihood model
###########################################################

# Author: Stuart Graham
# Created: 4/3/2020
# Last edited: 4/24/2020

# Load TreeNeighborhood package
devtools::load_all()

# Load package to time functions
library(tictoc)

# Source general functions from TreeNeigborhood package
#source("R/utils.R")

# Source data rearrangement functions for likelihood model
#source("R/likelihood_model_functions.R")

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
neighbors <- neighborhoods_all(mapping, 10)

# Remove focals whose neighborhood overlaps stand boundary
neighbors <- neighbors %>%
  filter(x_coord >= 10 & x_coord <= 90 & y_coord >= 10 & y_coord <= 90)

###################################
# Part 4. Calculating annual growth
###################################

# Calculate annual growth for all trees
growth <- growth_summary(tree)

# Produces some NA values. Remove those caused by a tree being measured only once,
# which would result in an annual growth measurement of NA
growth <- growth[growth$first_record != growth$last_record, ]

# Check for any remaining NA annual growth values
sum(is.na(growth$annual_growth))

# Investigate NaN values of size_corrected_growth
no_grow <- growth[is.nan(growth$size_corr_growth), ]
max(no_grow$annual_growth)

# NaN values caused by negative annual growth values. These must be caused by
# measurement inaccuracies because trees can't shrink. Converting them to 0s
# would change the distribution of the data in a biased way. Better to treat
# them as missing data and just remove them
growth <- growth[growth$annual_growth >= 0, ]

# Check data one more time
sum(is.nan(growth$size_corr_growth))

# Calculate radial growth to be consistent with Fortunel and Canham
growth$radial_growth <- growth$annual_growth / 2

################################################
# Part 5. Combining growth and neighborhood data
################################################

# Extract required columns from growth data
growth_cols <- growth[, c("tree_id", "midpoint_size", "radial_growth")]

# Join growth and neighborhood data. Use inner join because there will be
# no growth data for focals measured only once or with negative growth
full <- inner_join(neighbors, growth_cols, by = "tree_id")

################################
# Part 6. Exploring sample sizes
################################

# Find number of focals per species
focal_summ <- full %>% group_by(tree_id) %>% summarize(species = species[1])
table(focal_summ$species)

# Find number of competitors per species for each focal species
comp_summ <- function(data, focal_sps) {
  single_sps <- data[data$species == focal_sps, ]
  output <- table(single_sps$sps_comp)
  output
}
comp_summ(full, "TSME")

# TSME HAS ONLY 4 COMPETITOR SPECIES - TRY THIS FIRST!!!

###############################
# Part 7. Likelihood model TSME
###############################

# Subset to TSME focals
tsme <- droplevels(full[full$species == "TSME", ])

# Check size of dataset
table(tsme$sps_comp)
length(unique(tsme$tree_id))

# Create competitor species list
sps_list <- unique(tsme$sps_comp) # lambdas will be in order of sps_list

# Create NCI (neighborhood crowding index) function
NCI_calc <- function(neighbors, alpha, beta, lambda1, lambda2, lambda3, lambda4) {
  # Put lambdas into a vector to help ensure we use the right one
  lambda_list <- c(lambda1, lambda2, lambda3, lambda4)
  # Determine which species are competitors
  comp_sps <- which(sps_list %in% neighbors$sps_comp)
  # Initiate nci value
  nci_tot <- 0
  # Loop through competitor species
  for(sps in comp_sps) {
    ngb_sub <- neighbors[neighbors$sps_comp == sps_list[sps], ]
    # Calculate nci for this competitor species
    nci_sub <- sum(lambda_list[sps] * ((ngb_sub$size_comp ^ alpha) /
                                     (ngb_sub$prox ^ beta)))
    # Add to overall nci value
    nci_tot <- nci_tot + nci_sub
  }
  return(nci_tot)
}

# Create growth prediction function
growth_pred <- function(nbhd_data, gmax, X0, Xb, alpha, beta,
                        lambda1, lambda2, lambda3, lambda4, C, D){
  ids <- unique(nbhd_data$tree_id)
  growth <- rep(NA, times = length(ids))
  for(id in 1:length(ids)){
    neighbors <- nbhd_data[nbhd_data$tree_id == ids[id], ]
    growth[id] <- gmax * exp((-1/2) * ((log(neighbors$dbh[1] / X0) / Xb) ^ 2)) *
      exp(-C * (NCI_calc(neighbors, alpha, beta, lambda1, lambda2, lambda3, lambda4) ^ D))
  }
  results <- data.frame(ids, growth)
  return(results)
}

# Create negative log-likelihood function to optimize
lkhd <- function(par){
  gmax <- par[1]
  X0 <- par[2]
  Xb <- par[3]
  alpha <- par[4]
  beta <- par[5]
  C <- par[6]
  D <- par[7]
  sigma <- par[8]
  #lambda1 <- par[9]
  #lambda2 <- par[10]
  #lambda3 <- par[11]
  #lambda4 <- par[12]
  if(sigma<0) {return(Inf)}
  if(gmax<0 | X0<0) {return(Inf)}
  if(Xb<0 | Xb>20) {return(Inf)}           # from Uriarte: Xb[0,20]
  if(C<0 | C>10) {return(Inf)}              # from Uriarte: C [0,10] 
  if(D<1 | D>5) {return(Inf)}              # from Uriarte: D [1,5]   
  if(alpha<0 | alpha>4) {return(Inf)}      # from Uriarte: alpha [0,4]
  if(beta<0 | beta>4) {return(Inf)}
  #if(lambda1<0 | lambda1>1) {return(Inf)}
  #if(lambda2<0 | lambda2>1) {return(Inf)}
  #if(lambda3<0 | lambda3>1) {return(Inf)}
  #if(lambda4<0 | lambda4>1) {return(Inf)}
  #pred <- growth_pred(tsme, gmax, X0, Xb, alpha, beta, lambda1, lambda2, lambda3, lambda4, C, D)[, 2]
  pred <- growth_pred(tsme, gmax, X0, Xb, alpha, beta, lambda1=1, lambda2=1, lambda3=1, lambda4=1, C, D)[, 2]
  NLL <- (-sum(dnorm(tsme$radial_growth, mean = pred, sd = sigma, log = T)))
  return(NLL)
}

# Fit the model - TAKES ABOUT 10 MINUTES TO RUN
tic("model fitting")
fit <- optim(par = c(5, 100, 8, 3, 3, 0.5, 2, 5, 0.5, 0.5, 0.5, 0.5), fn = lkhd, method = "SANN")
toc()

# View parameters and compare to set parameter values
par_names <- c("gmax", "X0", "Xb", "alpha", "beta", "lambda1", "lambda2", "lambda3", "lambda4", 
               "C", "D", "sigma")
initial_val <- c(5, 100, 8, 3, 3, 0.5, 0.5, 0.5, 0.5, 0.5, 2, 5)
optim_val <- fit$par
result <- data.frame(par_names, initial_val, optim_val)
View(result)

# No lambdas - TAKES ABOUT 14 MINUTES TO RUN
tic("model fitting")
fit <- optim(par = c(5, 100, 8, 3, 3, 0.5, 2, 5), fn = lkhd, method = "SANN")
toc()

par_names <- c("gmax", "X0", "Xb", "alpha", "beta", "C", "D", "sigma")
initial_val <- c(5, 100, 8, 3, 3, 0.5, 2, 5)
optim_val <- fit$par
result <- data.frame(par_names, initial_val, optim_val)
View(result)
result$attempt3 <- fit$par

# Parameter values still change every time

# Save optimized parameters
params <- fit$par[1:7]

# Predict growth using model parameters
predictions <- growth_pred(tsme, params[1], params[2], params[3], params[4], params[5],
                    lambda1=1, lambda2=1, lambda3=1, lambda4=1, params[6], params[7])
names(predictions) <- c("tree_id", "predictions")

# Extract observations
observations <- tsme %>% group_by(tree_id) %>% summarize(growth = radial_growth[1])
names(observations)[2] <- "observations"

# Combine into table
comparison <- inner_join(predictions, observations, by = "tree_id")

# Plot predictions vs. observations
ggplot(comparison, aes(x = observations, y = predictions)) +
  geom_hex() +
  theme_bw() +
  ylim(0, max(c(max(comparison$predictions), max(comparison$observations)))) +
  xlim(0, max(c(max(comparison$predictions), max(comparison$observations)))) +
  geom_abline(intercept = 0, slope = 1)

# Calculate R^2
coef_det(comparison) # -4
plot(observations~predictions, data = comparison)

# Explore prediction for a single tree
one_tree2 <- tsme[tsme$tree_id == unique(tsme$tree_id)[1], ]
one_tree <- tsme[tsme$tree_id == comparison$tree_id[42], ]
growth_pred(one_tree, params[1], params[2], params[3], params[4], params[5],
            lambda1=1, lambda2=1, lambda3=1, lambda4=1, params[6], params[7])
params[1] # 5.64
exp((-1/2) * ((log(one_tree$dbh[1] / params[2]) / params[3]) ^ 2)) # 0.98
exp(-params[6] * (NCI_calc(one_tree2, params[4], params[5], 1, 1, 1, 1) ^ params[7])) # 9.3e-82
# neighborhood effect is definitely the problem here

#######################################################
# Part 8. Applying simple models to single-species data
#######################################################

# Subset to TSME focals
sing_sp <- droplevels(full[full$species == "TSME", ])

# Extract growth of each individual
focals <- sing_sp %>% group_by(tree_id) %>% summarize(dbh = dbh[1], radial_growth = radial_growth[1])

# Visualize relationship between focal size and growth
plot(radial_growth ~ dbh, data = focals, xlim = c(0, max(focals$dbh)))

# Define negative log likelihood function to be minimized
simp_NLL <- function(par) {
  
  # Define parameters
  gmax <- par[1]    # Maximum annual growth
  X0 <- par[2]      # DBH at which maximum growth occurs
  Xb <- par[3]      # Breadth of size-effect function
  sigma <- par[4]
  
  # Prevent parameter values from becoming nonsensical (i.e. negative)
  if(sigma <= 0 | gmax < 0 | X0 < 0| Xb < 0) {NLL <- 10000000} else {
    
    # Calculate negative log likelihood
    NLL <- -sum(dnorm(focals$radial_growth,
                      # LINE BELOW IS MODEL FORMULA
                      mean = gmax * exp((-1/2) * ((log(focals$dbh / X0) / Xb) ^ 2)),
                      sd = sigma,
                      log = T))
  }
  return(NLL)
}

# Create function to predict growth using model parameters
simp_pred <- function(data, par){
  par[1] * exp((-1/2) * ((log(data$dbh / par[2]) / par[3]) ^ 2))
}

# Create set of starting values
gmax_vals <- seq(0.03, 0.15, 0.03)
X0_vals <- seq(450, 1000, 50)
Xb_vals <- seq(0.5, 5.5, 1)
gmax <- rep(gmax_vals, each = length(X0_vals) * length(Xb_vals))
X0 <- rep(rep(X0_vals, each = length(Xb_vals)), times = length(gmax_vals))
Xb <- rep(Xb_vals, times = length(gmax_vals) * length(X0_vals))
starting_vals <- data.frame(gmax, X0, Xb)

# Create empty matrix for optimized values and minimized NLL
optim_vals <- as.data.frame(matrix(NA, nrow = nrow(starting_vals), ncol = 5))
names(optim_vals) <- c("gmax_opt", "X0_opt", "Xb_opt", "sigma_opt", "NLL")

# Loop through sets of starting values, optimizing and extracting results
for(i in 1:nrow(starting_vals)) {
  fit <- optim(par = c(starting_vals[i, 1], starting_vals[i, 2], starting_vals[i, 3], 5),
               fn = simp_NLL, method = "SANN")
  optim_vals[i, 1:4] <- fit$par
  optim_vals[i, 5] <- fit$value
}

# Combine starting and optimized values
output <- cbind(starting_vals, optim_vals)

#=============================================
# Exploring results of varying starting values
#=============================================
# Some repeated optimizations were run on a different machine - results can be loaded
# with this line
#output <- read.csv("Data/abam_3_param.csv")

# Exploring results
hist(output$NLL)
plot(NLL~gmax_opt, data = output)
plot(NLL~X0_opt, data = output)
plot(NLL~Xb_opt, data = output)
plot(gmax_opt~gmax, data = output)
plot(X0_opt~X0, data = output)
plot(Xb_opt~Xb, data = output)

# Order by increasing NLL
output <- output %>% arrange(NLL)

# Extract parameters of best fitting model
params <- c(output[1, 4], output[1, 5], output[1, 6])

# Calculate coefficient of determination
1 - (sum((focals$radial_growth - simp_pred(focals, params)) ^ 2) / 
       sum((focals$radial_growth - mean(focals$radial_growth)) ^ 2))
# psme: 0.07, tsme: 0.007, thpl: 0.16, tshe: 0.16, cano: 0.06, abam: 0.02
# with X0 <- seq(200, 600, 50), psme R^2 = 0.08
# with X0 <- seq(450, 1000, 50), psme R^2 = 0.083

# Plot optimized function over data
plot(radial_growth ~ dbh, data = focals, xlim = c(0, max(focals$dbh)))
curve(params[1]*exp((-0.5)*((log(x/params[2])/params[3])^2)), from = 0, to = max(focals$dbh), add = T)

#################################
# Part 9. Trying no lambda models
#################################
# Use THPL data because it has a smaller sample size but got good R^2 in 
# size effect model

