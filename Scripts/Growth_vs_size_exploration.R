


# Load TreeNeighborhood package
devtools::load_all()

# Load tree measurement data
tree <- read.csv("Data/Cleaned_tree_growth_2017.csv", stringsAsFactors = F)

# How quickly do these trees grow?
total_growth <- growth_summary(tree)
tapply(total_growth$annual_growth, total_growth$species, mean, na.rm = T)
tapply(total_growth$annual_growth, total_growth$species, st_err)

# Calculate growth between each consecutive pair of censuses
tree_growth <- detailed_growth(tree)

# Calculate number censuses each tree id appears in
censuses <- tree_growth %>%
  group_by(tree_id) %>%
  summarize(num_censuses = n())

# Join number of censuses to tree growth
tree_growth <- left_join(tree_growth, censuses, by = "tree_id")

# Calculate slope of growth rate over time
growth_slopes <- tree_growth %>%
  filter(num_censuses > 2 & !is.na(annual_growth)) %>%
  group_by(tree_id) %>%
  summarize(species = species[1], start_dbh = dbh[1],
            slope = lm(annual_growth ~ year)$coefficients[2])

# Subset to a single species and plot against dbh
plot(NULL, xlim = c(0, max(growth_slopes$start_dbh)),
     ylim = c(-0.02, 0.02),
     ylab = "Slope", xlab = "Starting DBH",)
sps <- c("PSME", "TSHE", "TSME", "THPL", "ABAM", "CANO")
cols <- c("red", "pink", "light blue", "light green", "dark blue", "orange")
slope_results <- matrix(NA, ncol = 2, nrow = 6)

for(i in 1:length(sps)){
  sing_sp <- growth_slopes %>% filter(species == sps[i])
  mod <- lm(sing_sp$slope ~ sing_sp$start_dbh)
  slope_results[i, 1] <- summary(mod)$coefficients[2, 1]
  slope_results[i, 2] <- summary(mod)$coefficients[2, 4]
  size_mod <- function(x){(mod$coefficients[2] * x) + mod$coefficients[1]}
  curve(size_mod(x), from = min(sing_sp$start_dbh), to = max(sing_sp$start_dbh),
        add = T, col = cols[i])
}

# Add legend
legend("topleft", legend = c("PSME", "TSHE", "TSME", "THPL", "ABAM", "CANO"),
       fill = cols, cex = 0.8, bty = "n")

slope_results <- cbind(sps, as.data.frame(slope_results))
names(slope_results)[2:3] <- c("slope", "p-value")


sing_sp <- growth_slopes %>% filter(species == "CANO")
plot(sing_sp$slope ~ sing_sp$start_dbh)

tapply(growth_slopes$slope, growth_slopes$species, st_err)




# Load mapping data
mapping <- read.csv("Data/Cleaned_mapping_2017.csv", stringsAsFactors = F)

# Join mapping data to growth slopes
sub_gslopes <- left_join(growth_slopes, mapping, by = "tree_id")

# Remove tree within 15m of stand boundary
sub_gslopes <- sub_gslopes %>%
  filter(x_coord >= 15 & x_coord <= 85 & y_coord >= 15 & y_coord <= 85)

# Restrict to training 1
train_sub1 <- sub_gslopes %>% filter(y_coord <= 50 - (15 / 2)) # Training 1
train_sub2 <- sub_gslopes %>% filter(x_coord >= 50 + (15 / 2) & y_coord > 50 - (15 / 2)) # Training 1
#train_sub1 <- sub_gslopes %>% filter(x_coord <= 50 - (15 / 2)) # Training 2
#train_sub2 <- sub_gslopes %>% filter(y_coord <= 50 - (15 / 2) & x_coord > 50 - (15 / 2)) # Training 2
#train_sub1 <- sub_gslopes %>% filter(y_coord >= 50 + (15 / 2)) # Training 3
#train_sub2 <- sub_gslopes %>% filter(x_coord <= 50 - (15 / 2) & y_coord < 50 + (15 / 2)) # Training 3
#train_sub1 <- sub_gslopes %>% filter(x_coord >= 50 + (15 / 2)) # Training 4
#train_sub2 <- sub_gslopes %>% filter(y_coord >= 50 + (15 / 2) & x_coord < 50 + (15 / 2)) # Training 4
training <- rbind(train_sub1, train_sub2)

tapply(training$slope, training$species.x, mean)
tapply(training$slope, training$species.x, st_err)
table(training$species.x)
