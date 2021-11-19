#' Select neighborhood size
#' 
#' Fits tree growth models using neighborhood metrics based on differing
#' neighborhood sizes, calculating mean square error of each model, thereby
#' providing useful data on the appropriate neighborhood size for analysis. It
#' is strongly recommended that this sensitivity analysis or a similar one is
#' used to select a neighborhood size for modeling because the appropriate 
#' neighborhood size varies between study systems and ecological processes.
#' 
#' This function uses a data frame of tree coordinates (\code{map_data}) to
#' calculate a series of metrics to describe the neighborhoods around each tree.
#' It does this by applying the \code{neighborhoods} and then the
#' \code{neighborhood_summary} function. It is therefore necessary to provide,
#' under the argument \code{dens_type}, instruction for the method that
#' \code{neighborhood_summary} should use to calculate densities. This process
#' is repeated for each of the neighborhood sizes specified by the argument
#' \code{radii}. To prevent edge effects, trees whose neighborhood overlaps the
#' stand boundary are excluded from modeling (hence the need to provide stand
#' sizes with the arguments \code{max_x} and \code{max_y}). Boundary overlap is
#' determined using the largest neighborhood size in \code{radii} to ensure
#' the number of focal trees is constant across neighborhood sizes.
#' 
#' The neighborhood metrics are then combined with the growth rates of their
#' focal tree, specifically the \code{size_corr_growth} output by 
#' \code{growth_summary} (any trees with illogical negative annual growth rates
#' are excluded as focal trees). If the user also provides abiotic data for the
#' stands using the argument \code{abiotic_data}, these are also joined to the
#' neighborhood metric and growth data. The result is a design matrix that is
#' given to the \code{growth_model} function. This entire process is repeated
#' for each neighborhood size (\code{radii}) and \code{focal_sps}, with the 
#' mean square error of each resulting model being recorded. The function
#' outputs a list containing a data frame of all mean square error values and
#' a plot of mean square error vs. neighborhood size for each \code{focal_sps}.
#' 
#' @param radii Vector of the neighborhood radii to try.
#' @param map_data Data frame containing tree coordinates.
#' @param growth_data Data frame containing repeated measurements of DBH for
#' trees in \code{map_data}
#' @param abiotic_data Data frame containing abiotic data for each of the stands
#' in \code{map_data}. Must contain a stand_id column. All other columns will be
#' included as explanatory variables in the fitted models. This argument is
#' optional.
#' @param focal_sps Vector of species names (i.e. character strings) for which
#' a neighborhood size should be selected.
#' @param dens_type Method for calculating tree densities in each neighborhood.
#' Must take one of the values: \code{"raw"}, \code{"proportional"}, or
#' \code{"angular"}. See \code{?neighborhood_summary} for details.
#' @param max_x Maximum expected x coordinate (i.e. should be 100 if the stands
#' are 100 x 100 m).
#' @param max_y Maximum expected y coordinate.
#' @param rare_sps Minimum number of interactions a competitor species must
#' appear in to remain separate from the "RARE" category in the growth models
#' (see \code{?growth_model} for details). If not specified, this argument will
#' default to a value of 100 interactions.
#' @return A list containing a number of named elements (exact number will be
#' the number of \code{focal_sps} plus one:
#' \itemize{
#' \item \code{mse_vals} is a data frame containing the mean square error of
#' the fitted model for each \code{focal_sps} using each neighborhood size
#' listed in \code{radii}.
#' \item each additional element of the list will be a plot of model mean square
#' error vs. neighborhood size for one of the \code{focal_sps}.
#' }
#' @examples
#' See vignette "Selecting neighborhood size"
#' @export
#' @importFrom magrittr %>%
#' @import dplyr

select_nbhd_size <- function(radii, map_data, growth_data,
                             abiotic_data = NULL, focal_sps, dens_type,
                             max_x, max_y, rare_sps = 100){
  
  # Create data frame to store results
  nb_rad_comp <- as.data.frame(matrix(NA, nrow = length(radii),
                                      ncol = length(focal_sps) + 1))
  names(nb_rad_comp) <- c("radius", paste0(focal_sps, "_mse"))
  nb_rad_comp$radius <- radii
  
  # Loop through the neighborhood radii
  for(i in 1:length(radii)){
    
    # Select neighborhood radius
    nb_rad <- radii[i]
    
    # Construct neighborhoods
    nbhds <- neighborhoods(map_data, radius = nb_rad)
    
    # Describe neighborhoods
    nbhd_summ <- neighborhood_summary(nbhds, id_column = "tree_id",
                                      radius = nb_rad, densities = dens_type)
    
    # Combine neighborhoods with their summaries
    nbhds <- nbhds %>%
      left_join(nbhd_summ, by = "tree_id")
    
    # Remove trees whose max neighborhood overlaps a plot boundary
    nbhds <- nbhds %>%
      filter(x_coord >= max(radii) & x_coord <= max_x - max(radii) &
               y_coord >= max(radii) & y_coord <= max_y - max(radii))
    
    # Calculate annual growth for all trees
    growth <- growth_summary(growth_data)
    
    # Remove trees that were only measured once and/or had negative growth
    growth <- growth %>%
      filter(first_record != last_record & annual_growth >= 0)
    
    # Add growth data to neighborhoods
    nbhds <- nbhds %>%
      inner_join(growth %>% select(tree_id, size_corr_growth),
                 by = "tree_id")
    
    # Add plot abiotic data, if provided
    if(!is.null(abiotic_data)){
      nbhds <- nbhds %>%
        left_join(abiotic_data, by = "stand_id")
    }
    
    # Loop through focal species
    for(sps in focal_sps){
      
      # Subset nbhds to focal species
      one_sps <- nbhds %>%
        filter(species == sps)
      
      # Drop columns not needed for the model
      one_sps <- one_sps %>%
        select(-c(stand_id, species, dbh, abh, x_coord, y_coord,
                  id_comp, abh_comp))
      
      # Run model
      if(dens_type == "angular"){
        mod <- growth_model(one_sps, outcome_var = "size_corr_growth",
                            rare_comps = rare_sps,
                            density_suffix = "_angle_sum")
      } else{
        mod <- growth_model(one_sps, outcome_var = "size_corr_growth",
                            rare_comps = rare_sps,
                            density_suffix = "_density")
      }
      
      # Store mean square error
      nb_rad_comp[i, paste0(sps, "_mse")] <- mod$mod_coef$mse[1] 
    }
  }
  
  # Create list to store overall output
  output_list <- as.list(rep(NA, times = length(focal_sps) + 1))
  names(output_list) <- c("mse_vals", paste0(focal_sps, "_plot"))
  
  # Make table of mse values the first element of output list
  output_list$mse_vals <- nb_rad_comp
  
  # Make a graph for each species, storing in output list
  for(a in 1:length(focal_sps)){
    
    # the local() function is used so that ggplot code is forced to evaluate
    # locally (i.e. within the loop)
    output_list[[a + 1]] <- local({
      local_a <- a
      sps_plot <- ggplot(data = nb_rad_comp,
                         aes(x = radius,
                             y = get(paste0(focal_sps[local_a], "_mse")))) +
        geom_line(col = "green") +
        labs(x = "Neighborhood radius (m)", y = "Mean square error") +
        theme_classic()
    })
  }
  
  # Return the output list
  return(output_list)
}
