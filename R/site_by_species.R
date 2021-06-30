

site_by_species <- function(neighbors, id_column, abundance = F){
  
  # Create ID column
  neighbors$id <- 1:nrow(neighbors)
  
  # Convert neighbors to wide format
  wide_nbhds <- neighbors %>%
    pivot_wider(id_cols = c(id, all_of(id_column)),
                names_from = sps_comp, names_sort = T,
                values_from = abh_comp, values_fill = 0) %>%
    select(-id)
  
  # Calculate quantities of each species
  if(abundance == T){    # abundance
    wide_nbhds <- wide_nbhds %>%
      group_by(get(id_column)) %>%
      summarize(across(.fn = ~ sum(. > 0)))
  } else if(abundance == F){    # presence/absence
    wide_nbhds <- wide_nbhds %>%
      group_by(get(id_column)) %>%
      summarize(across(.fn = ~ max(. > 0)))
  }
  
  # Reformat so that site names are row names
  wide_nbhds <- as.data.frame(wide_nbhds)
  output <- wide_nbhds[, 3:ncol(wide_nbhds)]
  row.names(output) <- wide_nbhds[, 1]
  
  # Return output
  return(output)
  
}
