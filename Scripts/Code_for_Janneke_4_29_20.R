#####################################
# Modeling tree size effect on growth
#####################################

# Author: Stuart Graham
# Created: 4/29/2020

# Load required packages
library(dplyr)

# Load neighborhood data
dat <- read.csv("Data/lkhd_data.csv", stringsAsFactors = F)

# Order data by tree_id to avoid problems later
dat <- dat %>% arrange(tree_id)

# Subset to a single focal species
sing_sp <- dat[dat$species == "THPL", ]

# Extract annual growth of each focal individual
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

# Optimize - values provided are gmax, X0, Xb, and sigma in that order
parameter_fits <- optim(par = c(0.5, 50, 10, 5), fn = simp_NLL, method = "SANN")
parameter_fits$par # View fitted values
parameter_fits$value #-1454

# Save optimized parameters
params <- parameter_fits$par[1:3]

# Create function to predict growth using model parameters
simp_pred <- function(data, par){
  par[1] * exp((-1/2) * ((log(data$dbh / par[2]) / par[3]) ^ 2))
}

# Plot predictions vs. observations
plot(y = focals$radial_growth, x = simp_pred(focals, params),
     ylab = "Observations", xlab = "Predictions",
     ylim = c(0, max(c(max(simp_pred(focals, params)), max(focals$radial_growth)))),
     xlim = c(0, max(c(max(simp_pred(focals, params)), max(focals$radial_growth)))))
abline(a = 0, b = 1)

# Calculate coefficient of determination
1 - (sum((focals$radial_growth - simp_pred(focals, params)) ^ 2) / 
       sum((focals$radial_growth - mean(focals$radial_growth)) ^ 2))

# Plot optimized function over data
plot(radial_growth ~ dbh, data = focals)
curve(0.11*exp((-0.5)*((log(x/257)/1.99)^2)), from = 0, to = max(focals$dbh), add = T)

#############################################
# Janneke's advice - try many starting values
#############################################

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

##############################################
# Exploring results of varying starting values
##############################################
# Some repeated optimizations were run on a different machine - results loaded below

# Load psme data
output <- read.csv("Data/thpl_3_param.csv")

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

################
# C and D models
################

# Subset to a single focal species
sing_sp <- dat[dat$species == "THPL", ]

# Extract annual growth of each focal individual
focals <- sing_sp %>% group_by(tree_id) %>% summarize(dbh = dbh[1], radial_growth = radial_growth[1])

# Define NCI function
nci <- function(neighbors){
  sum(neighbors$size_comp / neighbors$prox)
}

# Create growth prediction function
growth_pred <- function(nbhd_data, gmax, X0, Xb, C, D){
  ids <- unique(nbhd_data$tree_id)
  pred_grow <- rep(NA, times = length(ids))
  for(id in 1:length(ids)){
    neighbors <- nbhd_data[nbhd_data$tree_id == ids[id], ]
    pred_grow[id] <- gmax * exp((-1/2) * ((log(neighbors$dbh[1] / X0) / Xb) ^ 2)) *
      exp(-C * (nci(neighbors) ^ D))
  }
  results <- data.frame(ids, pred_grow, stringsAsFactors = F)
  return(results)
}

# Define negative log likelihood function to be minimized
c_and_d_NLL <- function(par) {
  
  # Define parameters
  gmax <- par[1]    # Maximum annual growth
  X0 <- par[2]      # DBH at which maximum growth occurs
  Xb <- par[3]      # Breadth of size-effect function
  C <- par[4]
  D <- par[5]
  sigma <- par[6]
  
  # Prevent parameter values from becoming nonsensical
  if(sigma < 0) {return(Inf)}
  if(gmax < 0 | X0 < 0) {return(Inf)}
  if(Xb < 0 | Xb > 20) {return(Inf)}
  if(C < 0 | C > 10) {return(Inf)}
  if(D < 1 | D > 5) {return(Inf)}
  
  # Make growth predictions
  pred <- growth_pred(sing_sp, gmax, X0, Xb, C, D)
  
  # Join predictions to observations by tree_id
  combined <- inner_join(focals, pred, by = c("tree_id" = "ids"))
  
  # Calculate negative log likelihood
  NLL <- -sum(dnorm(combined$radial_growth, mean = combined$pred_grow, sd = sigma, log = T))
  
  # Return value
  return(NLL)
}

# Optimize - values provided are gmax, X0, Xb, C, D, and sigma in that order
tic("model fitting")
fit <- optim(par = c(0.137, 256, 1.95, 2, 2, 5), fn = c_and_d_NLL, method = "SANN")
toc()

fit$par # View fitted values
fit$value #-1454
pred <- growth_pred(sing_sp, 0.137, 256, 1.95, 0.0002, 1)[, 2]
-sum(dnorm(focals$radial_growth, mean = pred, sd = 0.1, log = T))
# Extract parameters of best fitting model
params <- fit$par

# Calculate coefficient of determination
1 - (sum((focals$radial_growth - growth_pred(sing_sp, params[1], params[2], params[3], params[4], params[5])[, 2]) ^ 2) / 
       sum((focals$radial_growth - mean(focals$radial_growth)) ^ 2))

1 - (sum((focals$radial_growth - growth_pred(sing_sp, 0.137, 256, 1.95, 0.0002, 1)[, 2]) ^ 2) / 
       sum((focals$radial_growth - mean(focals$radial_growth)) ^ 2))

# Calculate coefficient of determination
params <- c(0.137, 256, 1.95)
1 - (sum((focals$radial_growth - simp_pred(focals, params)) ^ 2) / 
       sum((focals$radial_growth - mean(focals$radial_growth)) ^ 2))
x <- growth_pred(sing_sp, 0.137, 256, 1.95, 0, 1)[, 2] - simp_pred(focals, params)

problems <- which(x != 0)
prob_ids <- focals$tree_id[problems]
prob_dat <- sing_sp[sing_sp$tree_id %in% prob_ids, ]

growth_pred <- function(nbhd_data, gmax, X0, Xb, C, D){
  ids <- unique(nbhd_data$tree_id)
  growth <- rep(NA, times = length(ids))
  for(id in 1:length(ids)){
    neighbors <- nbhd_data[nbhd_data$tree_id == ids[id], ]
    growth[id] <- gmax * exp((-1/2) * ((log(neighbors$dbh[1] / X0) / Xb) ^ 2)) *
      exp(-C * (nci(neighbors) ^ D))
  }
  results <- data.frame(ids, growth)
  return(results)
}
simp_pred <- function(data, par){
  par[1] * exp((-1/2) * ((log(data$dbh / par[2]) / par[3]) ^ 2))
}

focals2 <- prob_dat %>% group_by(tree_id) %>% summarize(dbh = dbh[1], radial_growth = radial_growth[1])
