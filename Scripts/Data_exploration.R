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
cleanData <- read.csv("2013_cleaned.csv")

# Remove columns not required for exploratory analysis to keep things tidy
unneededCols <- which(names(cleanData) %in% c("Canopy_class", "Tree_vigor", 
                                              "Crown_ratio", "Main_stem", "Rooting", 
                                              "Crown_pct", "Tree_pct", "Lean_angle")) 
cleanData <- cleanData[,-unneededCols]

# Determine number of trees without mapping data
missingCoords <- droplevels(cleanData[is.na(cleanData$Xcoord) | 
                                        is.na(cleanData$Ycoord), ])
length(unique(missingCoords$TreeID))

# Remove these 162 trees from the dataset
cleanData <- droplevels(cleanData[!is.na(cleanData$Xcoord) &
                                    !is.na(cleanData$Ycoord), ])

# Final check for gaps in required data
sum(is.na(cleanData$TreeID))
sum(is.na(cleanData$StandID))
sum(is.na(cleanData$Species))
sum(is.na(cleanData$Year))
sum(is.na(cleanData$DBH_cm))

# All required data are present

#=============================
# Summarizing quantity of data
#=============================

# Determine how many trees we were measured at least twice (allows growth 
# estimation) and how these are distributed among species
indsPerSpecies <- cleanData %>% group_by(Species) %>% 
  summarize(numTrees = sum(table(TreeID) >= 2))
sum(indsPerSpecies$numTrees) # 8547 trees
#View(indsPerSpecies) # 10 species with at least 75 individual trees

# Determine how many times trees have been measured
hist(table(cleanData$TreeID)) # most have been measured 7 times

#==========================================
# Exploring how growth rate relates to size
#==========================================

# Remove trees that were only measured once
onceMeasured <- names(which(table(cleanData$TreeID) < 2))
cleanData <- cleanData[-which(cleanData$TreeID %in% onceMeasured),]

# Sort by year
cleanData <- cleanData[order(cleanData$TreeID, cleanData$Year),]

# Convert DBH to area at breast height (ABH) to get more meaningful growth
cleanData$ABH_cm2 <- pi*((cleanData$DBH_cm / 2) ^ 2)

# Calculate ABH difference between each row and the row below and add as new column
growth <- diff(cleanData$ABH_cm)
period <- diff(cleanData$Year)
annGrowth <- growth/period
cleanData$annGrowth <- c(annGrowth, NA)

# Remove growth estimates that are comparing different trees
meaninglessGrowth <- cumsum(as.vector(table(cleanData$TreeID)))
cleanData$annGrowth[meaninglessGrowth] <- NA

# Check for relationship when all species combined
plot(cleanData$annGrowth ~ cleanData$ABH_cm, ylim = c(0, 300))

# Check for relationship in western hemlock
TSHE <- cleanData[cleanData$Species == "TSHE",]
plot(TSHE$annGrowth ~ TSHE$ABH_cm, ylim = c(0, 300), ylab = "Annual Growth (cm2)",
     xlab = "Cross-sectional area at breast height (cm2)", 
     main = "Western hemlock")

# Check for relationship among western hemlock growing in the same stand
data <- TSHE[TSHE$StandID == "TA01",]
plot(data$annGrowth ~ data$ABH_cm, ylim = c(0, 100), ylab = "Annual Growth (cm2)",
     xlab = "Cross-sectional area at breast height (cm2)", 
     main = "Western hemlock in TA01")

data <- TSHE[TSHE$StandID == "PP17",]
plot(data$annGrowth ~ data$ABH_cm, ylim = c(0, 100), ylab = "Annual Growth (cm2)",
     xlab = "Cross-sectional area at breast height (cm2)", 
     main = "Western hemlock in PP17")
