#################################
# Exploring 2013 tree growth data
#################################

# This script loads the cleaned 2013 tree growth data (already includes mapping 
# information), checks and excludes missing data, explores the quantity of data,
#......

# Load required packages
library(dplyr)

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

# Remove 1930 lines with NA for DBH
cleanData <- cleanData[!is.na(cleanData$dbh), ]

#=============================
# Summarizing quantity of data
#=============================

# Determine how many trees were measured at least twice (allows growth 
# estimation) and how these are distributed among species
indsPerSpecies <- cleanData %>% group_by(species) %>% 
  summarize(numTrees = sum(table(treeid) >= 2))
sum(indsPerSpecies$numTrees) # 8376 trees

# Determine how many times trees have been measured
hist(table(cleanData$treeid)) # most have been measured at least 7 times

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

# Remove growth estimates that are comparing different trees
meaninglessGrowth <- cumsum(as.vector(table(cleanData$treeid)))
cleanData$annGrowth[meaninglessGrowth] <- NA

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
