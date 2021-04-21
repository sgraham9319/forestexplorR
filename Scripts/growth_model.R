# Load neighborhood data
nbhds <- read.csv("../../TreeGrowth/Data/Output_data/training1.csv")

# Remove unneeded columns from nbhds data
nbhds <- nbhds %>%
  select(species, sps_comp, prox, dbh_comp, grep("density", names(nbhds)),
         pet_mm, size_corr_growth)

# User inputs
training <- nbhds
focal_sps <- "PSME"
outcome_var <- "size_corr_growth"
# INTRA AS OPTION FOR FUNCTION???

# Extract focal species data and convert character variables to factors
sing_sp <- training %>%
  filter(species == focal_sps) %>%
  select(-species) %>%
  mutate(across(where(is.character), as.factor))

# Create model formula object
mod_form <- as.formula(
  paste(outcome_var, "~",
        paste(setdiff(names(sing_sp), outcome_var), collapse= "+")))

# Find factor column names
fctr_cols <- names(sing_sp)[sapply(sing_sp, is.factor)]

# Create list of factor columns
fctr_list <- list()
for(i in 1:length(fctr_cols)){
  fctr_list[[i]] <- sing_sp[, fctr_cols[i]]
  names(fctr_list)[i] <- fctr_cols[i]
}

# Create design matrix
dm <- model.matrix(mod_form, sing_sp,
                   contrasts.arg = lapply(fctr_list, contrasts, contrasts = F))

# Standardize variables except for first column (intercept)
dm[, 2:ncol(dm)] <- apply(dm[, 2:ncol(dm)], 2, z_trans)

# Change columns of NaNs (no variation) to zeros
dm[, which(is.nan(dm[1, ]))] <- 0

# Fit glmnet model
mod <- glmnet::cv.glmnet(x = dm, y = sing_sp$size_corr_growth,
                         family = "gaussian")

