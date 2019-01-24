#################################
# Exploring 2013 tree growth data
#################################

# This script loads the cleaned 2013 tree growth data (already includes mapping 
# information), checks and excludes missing data, explores the quantity of data,
#......

# Load package
devtools::load_all()

#==========================
# Loading and checking data
#==========================

# Load cleaned 2013 data
cleanData <- read.csv("../Data/Tree_growth_2017.csv")

# Remove columns not required for exploratory analysis to keep things tidy
unneededCols <- which(names(cleanData) %in% c("dbcode", "entity", "psp_studyid", "canopy_class", 
                                              "tree_vigor", "crown_ratio", "main_stem", "rooting", 
                                              "crown_pct", "tree_pct", "lean_angle")) 
cleanData <- cleanData[,-unneededCols]

# Determine number of trees without mapping data
#missingCoords <- droplevels(cleanData[is.na(cleanData$Xcoord) | 
#                                        is.na(cleanData$Ycoord), ])
#length(unique(missingCoords$TreeID))

# Remove these 162 trees from the dataset
#cleanData <- droplevels(cleanData[!is.na(cleanData$Xcoord) &
#                                    !is.na(cleanData$Ycoord), ])

# Final check for gaps in required data
sum(is.na(cleanData$treeid))
sum(is.na(cleanData$standid))
sum(is.na(cleanData$species))
sum(is.na(cleanData$year))
sum(is.na(cleanData$dbh))

#=======================
# Attaching mapping data
#=======================

# Load 2013 mapping data
mapping <- read.csv("../Data/Mapping_2013.csv")

# Extract unique tree IDs from growth data
treeIDs <- unique(cleanData$treeid) # 8803 tree IDs

# How many of these trees do we have mapping data for?
sum(treeIDs %in% mapping$TreeID) # 8416

# Add mapping data to growth data
cleanData$x_coord <- mapping[match(cleanData$treeid, mapping$TreeID), "Xcoord"]
cleanData$y_coord <- mapping[match(cleanData$treeid, mapping$TreeID), "Ycoord"]

#===========================================
# Explore new mapping data collected in 2017
#===========================================

# Determine how many remappings there are
table(cleanData$tree_status, cleanData$new_mapping) # 230 ingrowth and 166 old of
                                                    # a total of 7198 measured trees

# Isolate new mapping data
newMap <- cleanData[cleanData$new_mapping == "Y",]
update <- newMap[!is.na(newMap$x_coord), ]

table(cleanData$year)
ingrowth2017 <- cleanData[cleanData$tree_status == 2 & cleanData$new_mapping == "Y", "treeid"]
oldgrowth2017 <- cleanData[cleanData$tree_status == 1 & cleanData$new_mapping == "Y", "treeid"]
sum(ingrowth2017 %in% mapping$TreeID)
sum(oldgrowth2017 %in% mapping$TreeID)
table(mapping$Year)

#=========================================
# Investigate rows with no dbh measurement
#=========================================

# How many rows have NA for dbh measurement?
sum(is.na(cleanData$dbh)) # 1930

# Do all these rows have the dbh code for "missing"
any(cleanData$dbh_code[which(is.na(cleanData$dbh))] != "M") # Yes

# Do any rows have the dbh code for missing but have a dbh measurement?
sum(is.na(cleanData$dbh)) == sum(cleanData$dbh_code == "M") # No

# All of the trees without dbh measurements are said to be missing, but does this
# mean they were not found or that they were dead?
nodbh <- cleanData[cleanData$dbh_code == "M",]
table(nodbh$tree_status)  # Of 1930 trees, 1911 were dead (tree status = 6), and
                          # 19 were not found (tree status = 9)

# Check that dead trees were not recorded twice or ever came back to life
deadTrees <- cleanData[cleanData$tree_status == 6, "treeid"]
deadTreeDat <- cleanData[cleanData$treeid %in% deadTrees,]
deadTreeSumm <- deadTreeDat %>% group_by(treeid) %>% 
  summarize(numDead = sum(tree_status == 6),
            yearDead = year[which(tree_status == 6)],
            lastYearRecorded = max(year))
table(deadTreeSumm$numDead) # No trees recorded dead twice
any(deadTreeSumm$yearDead != deadTreeSumm$lastYearRecorded) # no dead trees come back to life

# All rows with no dbh measurement can be excluded from the analysis
cleanData <- cleanData[!is.na(cleanData$dbh), ]

#================================
# Calculating annual growth rates
#================================

# Calculate annual growth over entire measurement period for each tree
annual_growth <- overall_annual_growth(cleanData)

# How many cases of negative growth are there?
sum(growth$annual_growth < 0, na.rm = T) # 134
sum(is.na(growth$annual_growth))

# Calculate annual growth between 1977 and 1980
specific_annual_growth <- defined_period_annual_growth(data = cleanData,
                                                       begin = 1977, end = 1980)

#==========================================
# Exploring how growth rate relates to size
#==========================================

# Remove trees that were only measured once
onceMeasured <- names(which(table(cleanData$treeid) < 2))
cleanData <- cleanData[-which(cleanData$treeid %in% onceMeasured),]

# Sort by year
cleanData <- cleanData[order(cleanData$treeid, cleanData$year),]

# Convert DBH to area at breast height (ABH) to get more meaningful growth
cleanData$abh <- pi*((cleanData$dbh / 2) ^ 2)

# Calculate ABH difference between each row and the row below and add as new column
growth <- diff(cleanData$abh)
period <- diff(cleanData$year)
annGrowth <- growth/period
cleanData$annGrowth <- c(annGrowth, NA)

growthD <- diff(cleanData$dbh)
annGrowthD <- growthD/period
cleanData$annGrowthD <- c(annGrowthD, NA)

# Remove growth estimates that are comparing different trees
meaninglessGrowth <- cumsum(as.vector(table(cleanData$treeid)))
cleanData$annGrowth[meaninglessGrowth] <- NA
cleanData$annGrowthD[meaninglessGrowth] <- NA

# Check for relationship when all species combined
plot(cleanData$annGrowth ~ cleanData$abh, ylim = c(0, 300))

# Check for relationship in western hemlock
TSHE <- cleanData[cleanData$species == "TSHE",]
plot(TSHE$annGrowth ~ TSHE$abh, ylim = c(0, 300), ylab = "Annual Growth (cm2)",
     xlab = "Cross-sectional area at breast height (cm2)", 
     main = "Western hemlock")

# Check for relationship among western hemlock growing in the same stand
data <- TSHE[TSHE$standid == "TA01",]
plot(data$annGrowth ~ data$abh, ylim = c(0, 100), ylab = "Annual Growth (cm2)",
     xlab = "Cross-sectional area at breast height (cm2)", 
     main = "Western hemlock in TA01")

data <- TSHE[TSHE$standid == "PP17",]
plot(data$annGrowth ~ data$abh, ylim = c(0, 100), ylab = "Annual Growth (cm2)",
     xlab = "Cross-sectional area at breast height (cm2)", 
     main = "Western hemlock in PP17")

############################################################
# Exploring changes in growth rate with size for individuals
############################################################

# Create subset of one species in one stand
TA01TSHE <- droplevels(TSHE[TSHE$standid == "TA01",])

# Remove rows with NA for annual growth
TA01TSHE <- TA01TSHE[!is.na(TA01TSHE$annGrowth), ]

# Create overall plot
plot(TA01TSHE$annGrowth ~ TA01TSHE$abh, ylim = c(0, 100), ylab = "Annual Growth (cm2)",
     xlab = "Cross-sectional area at breast height (cm2)", 
     main = "Western hemlock in TA01")
plot(TA01TSHE$annGrowthD ~ TA01TSHE$dbh, ylab = "Annual Growth (cm)",
     xlab = "Cross-sectional area at breast height (cm)", 
     main = "Western hemlock in TA01")

# Identify unique tree IDs
uniqueIDs <- unique(TA01TSHE$treeid)

# Randomly sample 10 tree IDs
randIDs <- sample(uniqueIDs, size = 10)

# Subset to these 10 tree IDs
tenTrees <- TA01TSHE[TA01TSHE$treeid %in% randIDs,]

# Plot abh vs. annual growth
ggplot(data=tenTrees, aes(x=abh, y=annGrowth, group=treeid)) +
  geom_line() +
  geom_point()
ggplot(data=tenTrees, aes(x=dbh, y=annGrowthD, group=treeid)) +
  geom_line() +
  geom_point()

# Extract y-intercept and slope for linear model of annual growth ~ size for each
# tree as well as average size of the tree during the course of monitoring
y <- TA01TSHE %>% group_by(treeid) %>% summarize(y_intercept = unname(lm(annGrowth ~ abh)["coefficients"][[1]][1]),
                                                 slope = unname(lm(annGrowth ~ abh)["coefficients"][[1]][2]),
                                                 av_abh = mean(abh))

# Plot slope of model vs average size (expect large positive slope at small size 
# but small positive to zero slope at large size)
plot(y$slope ~ y$av_abh, ylim = c(0, max(y$slope, na.rm = T)))

# Should we be looking at dbh instead of abh?
plot(annGrowth ~ abh, data = TA01TSHE)
plot(annGrowthD ~ dbh, data = TA01TSHE)
