#################################
# Exploring 2013 tree growth data
#################################

# This script:
# 1) Loads the 2017 tree growth data
# 2) Investigates cases of NA dbh measurement to confirm they can be excluded
# from the analysis
# 3) Explores m2017 mapping data - 83 trees in growth data without mapping
# 4) Calculates annual growth rates and explores NA values
# 5) Explores how annual growth rate over the entire measurement period
# relate to tree size - finds expected positive relationship
# 6) Explores how annual growth rate between each census relates to size within
# individual trees - finds too much noise to detect this relationship

# Load package
devtools::load_all()

#=================================
# Loading and checking growth data
#=================================

# Load growth data
growth <- read.csv("../Data/Tree_growth_2017.csv", stringsAsFactors = F)

# Remove columns not required for exploratory analysis
growth <- growth %>%
  select(tree_id, stand_id, plot, species, tag, year, tree_status, dbh, 
         dbh_code, new_mapping)

# Identify trees that were never measured as 15 cm dbh or larger
small_trees_data <- growth %>%
  group_by(tree_id) %>% 
  summarize(max_size = max(dbh)) %>%
  filter(max_size < 15)

# Exclude small trees
growth <- growth[-which(growth$tree_id %in% small_trees_data$tree_id), ]

# Check for missing size data
sum(is.na(growth$dbh))
# 1930 cases of missing dbh

#=========================================
# Investigate rows with no dbh measurement
#=========================================

# Isolate missing dbh records
no_dbh <- growth %>%
  filter(is.na(dbh))

# Do all these rows have the dbh code for "missing"
sum(no_dbh$dbh_code == "M") == nrow(no_dbh)

# Do any rows have dbh code for missing but have a dbh measurement?
sum(growth$dbh_code == "M") == nrow(no_dbh)

# Were trees marked as missing actually missing or were they dead?
table(no_dbh$tree_status)
# 1911 dead (tree status = 6), 19 not found (tree status = 9)

# Confirm that no tree was recorded dead twice or came back to life
dead_trees <- growth[growth$tree_status == 6, "tree_id"]
dead_tree_summary <- growth %>%
  filter(tree_id %in% dead_trees) %>%
  group_by(tree_id) %>%
  summarize(times_rec_dead = sum(tree_status == 6),
            year_dead = year[which(tree_status == 6)],
            last_record = max(year))
max(dead_tree_summary$times_rec_dead)
# No trees recorded dead twice

any(dead_tree_summary$year_dead != dead_tree_summary$last_record)
# No dead trees come back to life

# Remove all rows with no dbh measurement from further analysis
growth <- growth %>%
  filter(!is.na(dbh))

#=======================
# Exploring mapping data
#=======================

# Load mapping data
mapping <- read.csv("../Data/Mapping_2017.csv", stringsAsFactors = F)

# Determine how many trees we do not have mapping for
tree_ids <- unique(growth$tree_id)
table(tree_ids %in% mapping$tree_id)
# Have mapping for 6892 trees, missing for 54 trees

# Isolate trees with no mapping data
no_map <- growth %>%
  filter(tree_id %in% mapping$tree_id == F)

#================================
# Calculating annual growth rates
#================================

# Calculate annual growth over entire measurement period for each tree
growth_summ <- growth_summary(growth)

# How many trees were measured only once?
sum(growth_summ$first_record == growth_summ$last_record) # 277

# Is this the only cause of NA in the growth columns?
sum(is.na(growth_summ$annual_growth)) # Yes
sum(is.na(growth_summ$size_corr_growth)) # No

# How many cases of negative growth are there?
sum(growth_summ$annual_growth < 0, na.rm = T) # 106

# Are being measured only once and negative growth the only causes of NA?
sum(is.na(growth_summ$size_corr_growth)) == 277 + 106 # Yes

# Calculate annual growth between 1977 and 1980
specific_annual_growth <- defined_period_annual_growth(data = growth,
                                                       begin = 1977, end = 1980)

#===============================================
# Exploring how growth rate relates to tree size
#===============================================

# Check for relationship when all species combined
plot(growth_summ$annual_growth, growth_summ$begin_size)

# Check for relationship in western hemlock
tshe <- growth_summ %>%
  filter(species == "TSHE")
plot(tshe$annual_growth, tshe$begin_size)

# Check for relationship among western hemlock growing in the same stand
tshe_ta01 <- tshe %>%
  filter(stand_id == "TA01")
plot(tshe_ta01$annual_growth, tshe_ta01$begin_size)

#=====================================================
# Exploring growth rate versus size within individuals
#=====================================================

# Calculate annual growth rate between each pair of consecutive censuses
# for each tree
partitioned_growth <- detailed_growth(growth)

# Identify unique tree IDs
tree_ids <- unique(partitioned_growth$tree_id)

# Randomly sample 10 tree IDs
rand_ids <- sample(tree_ids, size = 10)

# Subset to these 10 tree IDs
ten_trees <- partitioned_growth[partitioned_growth$tree_id %in% rand_ids,]

# Plot abh vs. annual growth
ggplot(data = ten_trees, aes(x = dbh, y = annual_growth, group = tree_id)) +
  geom_line() +
  geom_point()

# May be too much noise among measurements of the same tree. Check this by
# extracting slopes from models
growth_vs_size <- partitioned_growth %>%
  filter(!is.na(annual_growth)) %>%
  group_by(tree_id) %>%
  summarize(slope = unname(lm(annual_growth ~ dbh)["coefficients"][[1]][2]),
            av_dbh = mean(dbh))

# Plot slope of model vs average size (expect large positive slope at small size 
# but small positive to zero slope at large size)
plot(growth_vs_size$slope ~ growth_vs_size$av_dbh,
     ylim = c(0, max(growth_vs_size$slope, na.rm = T)))
summary(lm(slope ~ av_dbh, data = growth_vs_size))
# Relationship overall is positive - opposite to expectation