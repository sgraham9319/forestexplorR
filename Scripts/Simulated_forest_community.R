

# Simulate mapping data
num_trees <- 700 # Could be around 500 if just large trees, up to 1300 if also small trees
x_coord <- runif(num_trees, min = 0, max = 100)
y_coord <- runif(num_trees, min = 0, max = 100)
tree_id <- paste("A", 1:num_trees, sep = "")
stand_id <- rep("STANDA", times = num_trees)
species <- sample(c("A", "B", "C"), num_trees, replace = T)
dbh <- runif(num_trees, min = 5, max = 300)
mapping <- data.frame(tree_id, stand_id, species, x_coord, y_coord, dbh)
mapping <- mapping %>%
  mutate(size_cat = if_else(dbh >= 15, "regular", "small"))

# Plot stand map
plot(y_coord ~ x_coord, data = mapping)

# Create neighborhoods
neighbors <- neighborhoods(mapping, "STANDA", 10)

# Remove focals whose neighborhood overlaps stand boundary
neighbors <- neighbors %>%
  filter(x_coord >=10 & x_coord <= 90 & y_coord >=10 & y_coord <= 90 &
           size_cat == "regular")

#-------------------------
# Simulating annual growth
#-------------------------

# Define size effect parameters
gmax <- c(3, 1, 6)
X0 <- c(100, 250, 50)
Xb <- c(1, 1, 1)
growth_params <- data.frame(gmax, X0, Xb)
growth_params$species <- c("A", "B", "C")

# Plot size-growth relationships to check they are sensible
size_function <- function(x){
  gmax * exp((-1/2) * (log(x/X0) / Xb) ^ 2)
}
plot(NA, ylim = c(0, 6), xlim = c(0, 300))
for(i in 1:nrow(growth_params)){
  gmax <- growth_params$gmax[i]
  X0 <- growth_params$X0[i]
  Xb <- growth_params$Xb[i]
  curve(size_function(x), from = 1, to = 300, add = T)
}

# Define competition parameters
alpha <- c(1, 1, 1, 1, 1, 1, 1, 1, 1)
beta <- c(1, 1, 1, 1, 1, 1, 1, 1, 1)
lmd <- c(0.6, 0.9, 0.2, 0.3, 0.9, 0, 0.7, 1, 0.1)
comp_params <- data.frame(alpha, beta, lmd)
comp_params$species <- rep(c("A", "B", "C"), each = 3)
comp_params$sps_comp <- rep(c("A", "B", "C"), times = 3)

# Calculate NCI
nci_measurement <- neighbors %>%
  left_join(comp_params, by = c("species", "sps_comp")) %>%
  mutate(partial_nci = lmd * ((size_comp ^ alpha) / ((prox * 100) ^ beta))) %>%
  group_by(tree_id) %>%
  summarize(nci = sum(partial_nci))
  
# Make annual growth data
growth_estimates <- mapping %>%
  filter(x_coord >=10 & x_coord <= 90 & y_coord >=10 & y_coord <= 90 &
           size_cat == "regular") %>%
  left_join(growth_params, by = "species") %>%
  left_join(nci_measurement, by = "tree_id") %>%
  mutate(noise = rnorm(nrow(nci_measurement), sd = 0.01),
         size_only = gmax * exp((-1/2) * (log(dbh/X0) / Xb) ^ 2),
         comp_effect = exp((-0.5) * nci),
         annual_growth = (size_only * comp_effect) + noise)

# Check range of values seems reasonable
hist(growth_estimates$annual_growth)

# See how annual growth relates to size effects
plot(annual_growth ~ dbh, data = growth_estimates)
for(i in 1:nrow(growth_params)){
  gmax <- growth_params$gmax[i]
  X0 <- growth_params$X0[i]
  Xb <- growth_params$Xb[i]
  curve(size_function(x), from = 1, to = 300, add = T)
}

# See how annual growth differs between species
ggplot(data = growth_estimates, aes(x = dbh, y = annual_growth, col = species)) +
  geom_point()

#-------------------------------
# Creating tree growth data file
#-------------------------------

# Simulate growth data
census_years <- seq(1970, 2020, 5)
first_record <- sample(census_years, num_trees, replace = T, 
                       prob = c(0.8, rep((0.2 / (length(census_years) - 1)),
                                         times = length(census_years) - 1)))


0.1*50
hist(mapping$dbh)
x <- mapping %>% filter(dbh < 15 & dbh > 10)
