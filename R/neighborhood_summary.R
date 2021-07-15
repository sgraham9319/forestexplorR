#' Summarize neighborhoods
#'
#' Calculates species-specific densities and species richness for each 
#' neighborhood in the provided data frame. Densities can be provided as
#' area covered in the neighborhood, area covered per hectare, or the sum of
#' horizontal angles.
#'
#' @param neighbors Data frame of neighborhoods outputted by
#' \code{neighborhoods} function.
#' @param id_column Name of column in \code{neighbors} containing site names as
#' a string.
#' @param radius Numeric vector describing neighborhood radius in meters.
#' @param densities Character specifying the type of density measurements to
#' calculate - raw (m^2 per hectare), proportional (as proportions of overall
#' tree density), angular (angular size of trees).
#' @return Data frame containing a summary for each neighborhood in
#' \code{neighbors}.

neighborhood_summary <- function(neighbors, id_column, radius,
                                 densities = "raw"){
  
  # Create unique ID for each row
  neighbors$id <- 1:nrow(neighbors)
  
  # Calculate species richness of each neighborhood
  sps_rich <- neighbors %>%
    group_by(get(id_column)) %>%
    summarize(sps_richness = length(unique(sps_comp)))
  
  # Extract required columns
  dens <- neighbors %>%
    select(id, id_column, prox, dbh_comp, sps_comp, abh_comp)
  
  # Convert to wide format
  if(densities == "angular"){
    dens <- dens %>%
      mutate(sum_angle = atan(dbh_comp / (prox * 100))) %>%
      pivot_wider(id_cols = id,
                  names_from = sps_comp, names_sort = T,
                  values_from = sum_angle, values_fill = 0) %>%
      left_join(neighbors %>% select(id, id_column), by = "id") %>%
      select(-id) %>%
      group_by(get(id_column)) %>%
      select(-id_column) %>%
      summarize(across(.fn = sum))
  } else{
    dens <- dens %>%
      pivot_wider(id_cols = id,
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
  
  # Reorder rows to match input data frame
  output <- output[match(unique(neighbors[, id_column]),
                         output$`get(id_column)`), ]
  
  # Correct name of id column
  names(output)[1] <- id_column
  
  # Return output
  return(output)
  
}
