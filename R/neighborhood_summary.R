#' Summarize neighborhoods
#'
#' Calculates species-specific densities and species richness for each 
#' neighborhood in the provided data frame. Densities can be provided as
#' area covered in the neighborhood, area covered per hectare, or the sum of
#' horizontal angles.
#' 
#' This function contains an optional argument for handling edge effects. When
#' specified, the density and species richness values for neighborhoods that
#' overlap the stand boundary will be multiplied according to the fraction of 
#' their neighborhood that is missing from the stand. Species richness values
#' are rounded to the nearest whole number.
#'
#' @param neighbors Data frame of neighborhoods outputted by
#' \code{neighborhoods} function.
#' @param id_column Name of column in \code{neighbors} containing site names as
#' a string.
#' @param radius Numeric vector describing neighborhood radius in meters.
#' @param densities Character specifying the type of density measurements to
#' calculate - raw (m^2 per hectare), proportional (as proportions of overall
#' tree density), angular (angular size of trees).
#' @param edge_correction Boolean indicating whether edge correction should be
#' used.
#' @param x_limit maximum possible x-coordinate in the stand (only required if
#' \code{edge_correction = T})
#' @param y_limit maximum possible y-coordinate in the stand (only required if
#' \code{edge_correction = T})
#' @return Data frame containing a summary for each neighborhood in
#' \code{neighbors}.
#' @examples
#' # Create a neighborhoods object
#' nbhds <- neighborhoods(mapping, stands = "AB08", radius = 10)
#' 
#' # Summarize neighborhoods using angular densities
#' nbhd_summ <- neighborhood_summary(nbhds, id_column = "tree_id", radius = 10,
#'                                   densities = "angular")
#' 
#' # Using raw densities with edge correction
#' nbhd_summ <- neighborhood_summary(nbhds, id_column = "tree_id", radius = 10,
#'                                   densities = "raw", edge_correction = T,
#'                                   x_limit = 100, y_limit = 100)
#' @export
#' @importFrom magrittr %>%
#' @import dplyr

neighborhood_summary <- function(neighbors, id_column, radius,
                                 densities = "raw", edge_correction = F,
                                 x_limit = NULL, y_limit = NULL){
  
  # Create unique ID for each row
  neighbors$id <- 1:nrow(neighbors)
  
  # Calculate species richness of each neighborhood
  sps_rich <- neighbors %>%
    group_by(get(id_column)) %>%
    summarize(sps_richness = length(unique(sps_comp)),
              x_coord = x_coord[1],
              y_coord = y_coord[1])
  
  # Calculate proportion of neighborhood contained in stand if requested
  if(edge_correction == T){
    sps_rich <- suppressWarnings(nbhd_captured(sps_rich, x_limit, y_limit,
                                               radius))
  }
  
  # Remove x and y coordinate columns
  sps_rich <- sps_rich %>%
    select(-c(x_coord, y_coord))
  
  # Extract required columns for density calculations
  dens <- neighbors %>%
    select(id, id_column, prox, dbh_comp, sps_comp, abh_comp)
  
  # Convert to wide format
  if(densities == "angular"){
    dens <- dens %>%
      mutate(sum_angle = atan(dbh_comp / (prox * 100))) %>%
      tidyr::pivot_wider(id_cols = id,
                         names_from = sps_comp, names_sort = T,
                         values_from = sum_angle, values_fill = 0) %>%
      left_join(neighbors %>% select(id, id_column), by = "id") %>%
      select(-id) %>%
      group_by(get(id_column)) %>%
      select(-id_column) %>%
      summarize(across(.fn = sum))
  } else{
    dens <- dens %>%
      tidyr::pivot_wider(id_cols = id,
                         names_from = sps_comp, names_sort = T,
                         values_from = abh_comp, values_fill = 0) %>%
      left_join(neighbors %>% select(id, id_column), by = "id") %>%
      select(-id) %>%
      group_by(get(id_column)) %>%
      select(-id_column) %>%
      summarize(across(.fn = sum))
    
    # Convert densities to units of m^2 / hectare
    dens <- dens %>%
      rowwise() %>%
      mutate(across(2:ncol(dens), ~ ./circ_area(radius)))
  }
  
  # Calculate all_density
  dens <- dens %>%
    ungroup() %>%
    mutate(all = rowSums(dens[, 2:ncol(dens)]))
  
  # Add "angle_sum" or "density" to column names
  if(densities == "angular"){
    names(dens)[-1] <- paste(names(dens)[-1], "angle_sum", sep = "_")
  } else{
    names(dens)[-1] <- paste(names(dens)[-1], "density", sep = "_")
  }
  
  # Convert densities to proportions if requested
  if(densities == "proportional"){
    dens <- dens %>%
      rowwise() %>%
      mutate(across(2:(ncol(dens) - 1), ~ ./all_density))
  }
  
  # Join species richness to densities
  output <- sps_rich %>%
    left_join(dens, by = "get(id_column)")
  
  # Multiply species richness and densities for edge neighborhoods
  if(edge_correction == T & densities != "proportional"){
    output <- output %>%
      rowwise() %>%
      mutate(across(c(2, 4:ncol(output)), ~ .*(1 / prop_area_inc))) %>%
      mutate(sps_richness = round(sps_richness)) %>%
      select(-(prop_area_inc))
  }
  
  # Reorder rows to match input data frame
  output <- output[match(unique(neighbors[, id_column]),
                         output$`get(id_column)`), ]
  
  # Correct name of id column
  names(output)[1] <- id_column
  
  # Return output
  return(output)
  
}
