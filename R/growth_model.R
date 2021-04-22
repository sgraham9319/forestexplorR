#' Model tree growth
#' 
#' Fits a regularized linear regression model of tree growth for a single
#' species, returning the model object, growth predictions, the coefficient
#' of determination, and the fitted model coefficients.
#' 
#' All variables other than the indicated outcome variable are included as
#' explanatory variables in the model. Regularized regression model is fitted
#' using the glmnet package. Predictions, R-squared, and model coefficients
#' returned are all based on the \code{"lambda.1se"} model but all models with
#' different lambda values are included in the returned mod object - see glmnet 
#' documentation for details.
#' 
#' @param training A dataframe containing a separate row for each focal x
#' neighbor tree pair and a column for each variable to include in the
#' model, including the outcome variable. Outcome variable must be numeric,
#' other variables can be of any type, including character. Character variables
#' will be converted to factors before fitting the model.
#' @param outcome_var Name of column in \code{training} containing the outcome
#' variable, provided as a string.
#' @param focal_sps Name of species for which growth should be modeled,
#' provided as a string.
#' @return A list containing four elements: \code{mod} is the fitted glmnet
#' object, \code{obs_pred} is a dataframe containing observed and predicted 
#' growth, \code{R_squared} is the coefficient of determination, \code{mod_coef}
#' is a matrix containing the explanatory variables and their fitted
#' coefficients. These elements can be called with e.g. \code{$R_squared}.
#' 


growth_model <- function(training, outcome_var, focal_sps){
  
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
  
  # Calculate model predictions for training data
  preds <- predict(mod, newx = dm, s = "lambda.1se")
  
  # Combine with observations
  obs_pred <- cbind(sing_sp %>% select(size_corr_growth), preds)
  names(obs_pred) <- c("observations", "predictions")
  
  # Calculate coefficient of determination
  R_squared <- coef_det(obs_pred)
  
  # Extract model coefficients
  mod_coef <- as.matrix(coef(mod, s = "lambda.1se"))
  
  # Create output list
  output <- list(mod = mod, obs_pred = obs_pred, R_squared = R_squared,
                 mod_coef = mod_coef)
  
  # Return output
  output
}