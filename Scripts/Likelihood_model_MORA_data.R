###########################################################
# Neighborhood analysis of MORA data using likelihood model
###########################################################

# Author: Stuart Graham
# Created: 4/3/2020
# Last edited: 4/23/2020

# Load TreeNeighborhood package
devtools::load_all()

# Load package to time functions
library(tictoc)

# Source general functions from TreeNeigborhood package
source("R/utils.R")

# Source data rearrangement functions for likelihood model
source("R/likelihood_model_functions.R")

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
  lambda1 <- par[6]
  lambda2 <- par[7]
  lambda3 <- par[8]
  lambda4 <- par[9]
  C <- par[10]
  D <- par[11]
  sigma <- par[12]
  if(sigma<0) {return(Inf)}
  if(gmax<0 | X0<0) {return(Inf)}
  if(Xb<0 | Xb>20) {return(Inf)}           # from Uriarte: Xb[0,20]
  if(C<0 | C>1) {return(Inf)}              # from Uriarte: C [0,10] 
  if(D<1 | D>5) {return(Inf)}              # from Uriarte: D [1,5]   
  if(alpha<0 | alpha>4) {return(Inf)}      # from Uriarte: alpha [0,4]
  if(beta<0 | beta>4) {return(Inf)}
  if(lambda1<0 | lambda1>1) {return(Inf)}
  if(lambda2<0 | lambda2>1) {return(Inf)}
  if(lambda3<0 | lambda3>1) {return(Inf)}
  if(lambda4<0 | lambda4>1) {return(Inf)}
  pred <- growth_pred(tsme, gmax, X0, Xb, alpha, beta, lambda1, lambda2, lambda3, lambda4, C, D)[, 2]
  NLL <- (-sum(dnorm(tsme$radial_growth, mean = pred, sd = sigma, log = T)))
  return(NLL)
}

# Fit the model - TAKES ABOUT 10 MINUTES TO RUN
tic("model fitting")
fit <- optim(par = c(5, 100, 8, 3, 3, 0.5, 0.5, 0.5, 0.5, 0.5, 2, 5), fn = lkhd, method = "SANN")
toc()

# View parameters and compare to set parameter values
par_names <- c("gmax", "X0", "Xb", "alpha", "beta", "lambda1", "lambda2", "lambda3", "lambda4", 
               "C", "D", "sigma")
initial_val <- c(5, 100, 8, 3, 3, 0.5, 0.5, 0.5, 0.5, 0.5, 2, 5)
optim_val <- fit$par
result <- data.frame(par_names, initial_val, optim_val)
View(result)
