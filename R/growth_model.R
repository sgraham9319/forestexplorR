#' Model tree growth
#' 
#' Fits regularized linear regression models of tree growth for a single
#' species, returning the model object, its growth predictions, its coefficient
#' of determination when applied to the training set (and optionally a provided
#' test set), and the fitted model coefficients.
#' 
#' All variables in the user-provided training data other than species,
#' tree_id, and the indicated outcome variable are included as explanatory
#' variables in the model. Regularized regression model is fitted
#' using the glmnet package. Predictions, R-squared, and model coefficients
#' returned are all based on the \code{"lambda.1se"} model but all models with
#' different lambda values are included in the returned mod object - see glmnet 
#' documentation for details.
#' 
#' As the glmnet model fitting process is stochastic, the fitted model can
#' differ with each run. The \code{iterations} argument allows the user to
#' specify how many times the model should be fitted. If the model is fitted
#' more than once, the fitted coefficients for all models will be returned in
#' but only the model object for the best model (lowest cross-validated mean
#' square error) will be returned.
#' 
#' Rare competitor species can be grouped together as "OTHR" for modeling if
#' desired. The optional argument \code{rare_comps} is a number indicating the
#' minimum number of interactions a species must appear in to remain separate
#' from the "OTHR" category. If species-specific densities are included in
#' \code{training} those of rare competitors will be combined into a new 
#' "OTHR_density" variable.
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
#' @param iterations Number of times the model should be fitted.
#' @param rare_comps Minimum number of interactions a competitor species must
#' appear in to remain separate from the "OTHR" category (see details). If not
#' specified, no "OTHR" category will be created and all competitor species will
#' remain separate.
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

growth_model <- function(training, outcome_var, focal_sps, iterations = 1,
                         rare_comps = "none", test = NULL){
  
    # Extract focal species data
  sing_sp <- training %>%
    arrange(tree_id) %>%
    filter(species == focal_sps) %>%
    select(-species)
  
  # Create vector of tree ids then remove this column
  tree_ids <- sing_sp$tree_id
  sing_sp <- sing_sp %>%
    select(-tree_id)
  
  # Remove any density columns containing only zeros
  sing_sp <- sing_sp %>%
    select(-names(
      which(colSums(sing_sp[, grep("density", names(sing_sp))]) == 0)))
  
  # Handle rare competitors if requested
  if(rare_comps != "none"){
    
    # Make list of rare competitor species
    rare <- names(which(table(sing_sp$sps_comp) < rare_comps))
    
    # Change rare competitor species to OTHR
    sing_sp <- sing_sp %>%
      mutate(sps_comp = if_else(sps_comp %in% rare, "OTHR",
                                as.character(sps_comp)))
    
    # Create OTHR density column
    rare_dens_cols <- paste(rare, "density", sep = "_")
    sing_sp <- sing_sp %>%
      mutate(OTHR_density = apply(sing_sp %>%
                                    select(all_of(rare_dens_cols)), 1, sum)) %>%
      select(-all_of(rare_dens_cols))
    
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
    
    # Fit glmnet model
    mod <- glmnet::cv.glmnet(x = dm, y = sing_sp$size_corr_growth,
                             family = "gaussian")
    
    # Calculate model predictions for training data
    preds <- predict(mod, newx = dm, s = "lambda.1se")
    
    # Combine with observations
    obs_pred <- cbind(tree_ids, sing_sp %>% select(size_corr_growth), preds)
    names(obs_pred) <- c("tree_id", "observations", "predictions")
    
    # Get single prediction for each tree
    obs_pred <- obs_pred %>%
      group_by(tree_id) %>%
      summarize(observations = observations[1],
                predictions = mean(predictions))
    
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
    
    # Subset test to focal species
    ss_test <- test %>%
      arrange(tree_id) %>%
      filter(species == focal_sps) %>%
      select(-species)
    
    # Create vector of tree ids then remove this column
    test_ids <- ss_test$tree_id
    ss_test <- ss_test %>%
      select(-tree_id)
    
    # Handle rare competitors if requested
    if(rare_comps != "none"){
      
      # Change rare competitor species to OTHR
      ss_test <- ss_test %>%
        mutate(sps_comp = if_else(sps_comp %in% rare, "OTHR",
                                  as.character(sps_comp)))
      
      # Create OTHR density column
      ss_test <- ss_test %>%
        mutate(OTHR_density = apply(ss_test %>% select(all_of(rare_dens_cols)),
                                    1, sum)) %>%
        select(-all_of(rare_dens_cols))
      
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
    test_obs_pred <- cbind(test_ids, ss_test %>% select(size_corr_growth),
                           test_preds)
    names(test_obs_pred) <- c("tree_id", "observations", "predictions")
    
    # Get single prediction for each tree
    test_obs_pred <- test_obs_pred %>%
      group_by(tree_id) %>%
      summarize(observations = observations[1],
                predictions = mean(predictions))
    
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