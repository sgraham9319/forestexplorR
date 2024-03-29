#' Model tree growth
#' 
#' Fits regularized linear regression models of tree growth for a single
#' species, returning the model object, its growth predictions, its coefficient
#' of determination when applied to the training set (and optionally, a
#' user-provided test set), and the fitted model coefficients.
#' 
#' All variables in the user-provided training data other than tree_id and the
#' indicated outcome variable are included as explanatory variables in the
#' model. The regularized regression model is fitted using the glmnet package.
#' Predictions, R-squared, and model coefficients returned are all based on the
#' \code{"lambda.1se"} model but all models with different lambda values are
#' included in the returned mod object - see glmnet documentation for details.
#' 
#' As the glmnet model fitting process is stochastic, the fitted model can
#' differ with each run. The \code{iterations} argument allows the user to
#' specify how many times the model should be fitted. If the model is fitted
#' more than once, the fitted coefficients for all models will be returned
#' but only the model object for the best model (lowest cross-validated mean
#' square error) will be returned.
#' 
#' Rare competitor species can be grouped together as "RARE" for modeling if
#' desired. The optional argument \code{rare_comps} is a number indicating the
#' minimum number of interactions a species must appear in to remain separate
#' from the "RARE" category. If handling of rare competitor species is requested
#' and species-specific densities are included in \code{training}, the densities
#' of rare species can be summed together under "RARE_density" using the
#' optional argument \code{density_suffix}.
#' 
#' @param training A dataframe containing a separate row for each focal x
#' neighbor tree pair and a column for each variable to include in the
#' model, including the outcome variable. Outcome variable must be numeric.
#' Other variables can be of any type, including character. Any character
#' variables will be converted to factors before fitting the model.
#' @param outcome_var Name of column in \code{training} containing the outcome
#' variable, provided as a string.
#' @param iterations Number of times the model should be fitted.
#' @param rare_comps Minimum number of interactions a competitor species must
#' appear in to remain separate from the "RARE" category (see details). If not
#' specified, no "RARE" category will be created and all competitor species will
#' remain separate.
#' @param density_suffix Suffix of columns containing species-specific densities
#' e.g. if density columns are of the form \code{speciesA_dens}, this argument
#' should have the value \code{"_dens"}. This optional argument is only used
#' if the argument \code{rare_comps} is specified.
#' @param test An optional dataframe of test data that must be in exactly the
#' same format as the training data i.e. all the same columns with the same
#' names.
#' @return A list containing four or five elements:
#' \itemize{
#' \item \code{mod} is a fitted glmnet model
#' object - if \code{iterations > 1} this will be the model with the lowest
#' cross-validated mean square error
#' \item \code{obs_pred} is a dataframe containing
#' observed and predicted growth - if \code{iterations > 1} this will 
#' correspond to the best model
#' \item \code{R_squared} is the coefficient of
#' determination - if \code{iterations > 1} this will correspond to the best
#' model
#' \item \code{test_R_squared} is the coefficient of determination of the best
#' model on the test data (this element will not appear if no test data are
#' provided)
#' \item \code{mod_coef} is a data frame containing the fitted coefficients, 
#' cross-validated mean square error (\code{mse}), and coefficient of 
#' determination for each fitted model with rows in ascending order of 
#' \code{mse}.
#' }
#' @examples
#' # See vignette "Modeling tree growth and mortality"
#' @export
#' @importFrom magrittr %>%
#' @import dplyr

growth_model <- function(training, outcome_var, iterations = 1,
                         rare_comps = "none", density_suffix = "none",
                         test = NULL){
  
  # Order training by tree id
  sing_sp <- training %>%
    arrange(tree_id)
  
  # Create vector of tree ids then remove this column
  tree_ids <- sing_sp$tree_id
  sing_sp <- sing_sp %>%
    select(-tree_id)
  
  # Handle rare competitors if requested
  if(rare_comps != "none"){
    
    # Make list of rare competitor species
    rare <- names(which(table(sing_sp$sps_comp) < rare_comps))
    
    # Change rare competitor species to "RARE"
    sing_sp <- sing_sp %>%
      mutate(sps_comp = if_else(sps_comp %in% rare, "RARE",
                                as.character(sps_comp)))
    
    if(density_suffix != "none"){
      
      # Create rare density column
      rare_dens_cols <- paste(rare, density_suffix, sep = "")
      sing_sp <- sing_sp %>%
        mutate(RARE_density = apply(sing_sp %>%
                                      select(all_of(rare_dens_cols)), 1, sum)) %>%
        select(-all_of(rare_dens_cols))
       
    }
  }
  
  # Convert character variables to factors
  sing_sp <- sing_sp %>%
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
  
  # Iterate the model fitting
  for(i in 1:iterations){
    
    # Fit glmnet model and make predictions for training data
    mod <- glmnet::cv.glmnet(x = dm, y = sing_sp[, outcome_var],
                             family = "gaussian")
    preds <- predict(mod, newx = dm, s = "lambda.1se")
    
    # Combine with observations
    obs_pred <- cbind(tree_ids, sing_sp %>% select(all_of(outcome_var)), preds)
    names(obs_pred) <- c("tree_id", "observations", "predictions")
    
    # Get single prediction for each tree
    obs_pred <- obs_pred %>%
      group_by(tree_id) %>%
      summarize(observations = observations[1],
                predictions = mean(as.numeric(predictions)))
    
    # Calculate coefficient of determination
    R_squared <- coef_det(obs_pred)
    
    # Extract model coefficients
    new_coef <- as.matrix(coef(mod, s = "lambda.1se"))
    
    # Add model mean square error and R squared to coefficients
    mse <- mod$cvm[mod$lambda == mod$lambda.1se]
    new_coef <- rbind(new_coef, mse, R_squared)
    
    # Update best model output if new model is the best
    if(i == 1){
      best_mod <- mod
      final_obs_pred <- obs_pred
      best_R_squared <- R_squared
    } else if(mse < min(mod_coef["mse", ])){
      best_mod <- mod
      final_obs_pred <- obs_pred
      best_R_squared <- R_squared
    }
    
    # Add to model coefficients matrix
    if(i == 1){
      mod_coef <- new_coef
    } else{
      mod_coef <- cbind(mod_coef, new_coef)
    }
  }
  
  # Clean model coefficients table
  mod_coef <- as.data.frame(t(mod_coef))
  mod_coef <- mod_coef %>%
    select(-2) %>%
    arrange(mse)
  rownames(mod_coef) <- 1:nrow(mod_coef)
  
  # Evaluate fit to test data if requested
  if(!is.null(test)){
    
    # Order test by tree id
    ss_test <- test %>%
      arrange(tree_id)
    
    # Create vector of tree ids then remove this column
    test_ids <- ss_test$tree_id
    ss_test <- ss_test %>%
      select(-tree_id)
    
    # Handle rare competitors if requested
    if(rare_comps != "none"){
      
      # Change rare competitor species to "RARE"
      ss_test <- ss_test %>%
        mutate(sps_comp = if_else(sps_comp %in% rare, "RARE",
                                  as.character(sps_comp)))
      
      if(density_suffix != "none"){
        
        # Create OTHR density column
        ss_test <- ss_test %>%
          mutate(RARE_density = apply(ss_test %>% select(all_of(rare_dens_cols)),
                                      1, sum)) %>%
          select(-all_of(rare_dens_cols))
      }
    }
    
    # Convert character variables to factors
    ss_test <- ss_test %>%
      mutate(across(where(is.character), as.factor))
    
    # Find factor column names
    fctr_cols <- names(ss_test)[sapply(ss_test, is.factor)]
    
    # Create list of factor columns
    fctr_list <- list()
    for(i in 1:length(fctr_cols)){
      fctr_list[[i]] <- ss_test[, fctr_cols[i]]
      names(fctr_list)[i] <- fctr_cols[i]
    }
    
    # Create design matrix
    dm_test <- model.matrix(mod_form, ss_test,
                            contrasts.arg = lapply(fctr_list, contrasts,
                                                   contrasts = F))
    
    # Find any columns in dm missing from dm_test
    missing_cols <- setdiff(colnames(dm), colnames(dm_test))
    
    # Add these columns to dm_test
    for(column in missing_cols){
      dm_test <- cbind(dm_test, rep(0, times = nrow(dm_test)))
      colnames(dm_test)[ncol(dm_test)] <- column
    }
    
    # Reorder dm_test columns to match order of dm
    dm_test <- dm_test[, match(colnames(dm), colnames(dm_test))]
    
    # Standardize variables except for first column (intercept)
    dm_test[, 2:ncol(dm_test)] <- apply(dm_test[, 2:ncol(dm_test)], 2, z_trans)
    
    # Change columns of NaNs (no variation) to zeros
    dm_test[, which(is.nan(dm_test[1, ]))] <- 0
    
    # Calculate model predictions for test data
    test_preds <- predict(best_mod, newx = dm_test, s = "lambda.1se")
    
    # Combine with observations
    test_obs_pred <- cbind(test_ids, ss_test %>% select(all_of(outcome_var)),
                           test_preds)
    names(test_obs_pred) <- c("tree_id", "observations", "predictions")
    
    # Get single prediction for each tree
    test_obs_pred <- test_obs_pred %>%
      group_by(tree_id) %>%
      summarize(observations = observations[1],
                predictions = mean(as.numeric(predictions)))
    
    # Calculate coefficient of determination
    test_R_squared <- coef_det(test_obs_pred)
    
  }
  
  # Create output list
  if(is.null(test)){
    output <- list(mod = best_mod, obs_pred = final_obs_pred,
                   R_squared = best_R_squared, mod_coef = mod_coef)
  } else {
    output <- list(mod = best_mod, obs_pred = final_obs_pred,
                   R_squared = best_R_squared, 
                   test_R_squared = test_R_squared, mod_coef = mod_coef)
  }
  
  # Return output
  output
}
