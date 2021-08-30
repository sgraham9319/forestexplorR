#' Create site x species matrix
#' 
#' Takes in a neighborhoods object output by the \code{neighborhoods} function
#' and returns a site x species matrix where each site is a neighborhood. The 
#' user can specify whether the values in the matrix represent abundance or
#' presence/absence.
#' 
#' @param neighbors A neighborhoods object output by the \code{neighborhoods}
#' function.
#' @param id_column Name of column in \code{neighbors} containing site names as
#' a string.
#' @param abundance Boolean specifying whether an abundance or presence/absence
#' site x species matrix is desired (default is presence/absence).
#' @return A site x species matrix where each row represents a distinct site or
#' neighborhood and each column represents a tree species. Row names are the 
#' site names. If \code{abundance = F} presence is indicated by 1 and absence 
#' by 0. If \code{abundance = T} values represent the number of trees of that
#' species in the neighborhood (excluding the focal if neighborhoods are
#' centered on a focal tree).
#' @export
#' @importFrom magrittr %>%
#' @import dplyr

site_by_species <- function(neighbors, id_column, abundance = F){
  
  # Create unique ID for each row
  neighbors$id <- 1:nrow(neighbors)
  
  # Convert neighbors to wide format
  wide_nbhds <- neighbors %>%
    tidyr::pivot_wider(id_cols = c(id, all_of(id_column)),
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
  
  # Reorder rows to match input data frame
  output <- output[match(unique(neighbors[, id_column]), row.names(output)), ]
  
  # Return output
  return(output)
  
}
