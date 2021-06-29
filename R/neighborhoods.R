#' Create neighborhoods
#'
#' Takes a data frame containing mapping information of trees in one or more
#' stands and returns a data frame with a row for each focal-neighbor pair
#'
#' @param mapping Data frame containing tree coordinates.
#' @param stands Vector of names of stands for which neighborhoods are desired.
#' @param radius Numeric vector describing neighborhood radius in meters.
#' @return Neighborhood information for all focal trees in \code{mapping}.

neighborhoods <- function(mapping, stands = "all", radius) {
  
  # Create vector of stand IDs
  if(length(stands) == 1 & stands[1] == "all"){
    stand_list <- unique(mapping$stand_id)
  } else {
    stand_list <- stands
  }
  
  # Convert species column to factor
  mapping$species <- as.factor(mapping$species)
  
  # Check if requested stands appear in mapping data
  for(i in stand_list){
    if(i %in% mapping$stand_id == F){
      print(paste("Warning: skipping stand ", i, " because not found in
                  mapping data", sep = "'"))
    }
  }
  
  # Update stand list
  stand_list <- stand_list[stand_list %in% mapping$stand_id]
  
  # Loop through stands creating neighborhoods for each
  for(i in 1:length(stand_list)){
    
    # Subset to focal stand and calculate abh
    one_stand <- mapping %>%
      filter(stand_id == stand_list[i]) %>%
      mutate(abh = circ_area(dbh / 2)) %>%
      select(tree_id, stand_id, species, dbh, abh, x_coord, y_coord)
    
    # Create distance matrix
    dist_mat <- as.matrix(dist(one_stand[, c("x_coord", "y_coord")]))
    
    # Create competitor species matrix
    sps_mat <- matrix(one_stand$species, nrow = nrow(one_stand), ncol = nrow(one_stand))
    
    # Create competitor dbh matrix
    dbh_mat <- matrix(one_stand$dbh, nrow = nrow(one_stand), ncol = nrow(one_stand))
    
    # Create competitor abh matrix
    abh_mat <- matrix(one_stand$abh, nrow = nrow(one_stand), ncol = nrow(one_stand))
    
    # Create competitor tree_id matrix
    id_mat <- matrix(one_stand$tree_id, nrow = nrow(one_stand), ncol = nrow(one_stand))
    
    # Find values in distance matrix greater than radius or equal to zero - these
    # elements represent pairs of trees that are not in each others neighborhood or
    # are the same tree paired with itself
    non_comp <- which(dist_mat > radius | dist_mat == 0)
    
    # Convert values of these non-competing pairs to NA in all matrices
    dist_mat[non_comp] <- NA
    sps_mat[non_comp] <- NA
    dbh_mat[non_comp] <- NA
    abh_mat[non_comp] <- NA
    id_mat[non_comp] <- NA
    
    # Create distance, species, dbh, abh and size category vectors, excluding NAs
    prox <- dist_mat[!is.na(dist_mat)]
    sps_comp <- sps_mat[!is.na(sps_mat)]
    dbh_comp <- dbh_mat[!is.na(dbh_mat)]
    abh_comp <- abh_mat[!is.na(abh_mat)]
    id_comp <- id_mat[!is.na(id_mat)]
    
    # Combine vectors into competitor information data frame
    comp_dat <- data.frame(id_comp, prox, sps_comp, dbh_comp, abh_comp)
    
    # Get vector of how many rows (competitors) each focal tree should have
    repeats <- unname(apply(dist_mat, 2, non_na_len))
    
    # Repeat rows of one_stand appropriate number of times
    all_rows <- one_stand[rep(1:nrow(one_stand), repeats), ]
    
    # Bind competitor data to one_stand
    one_stand_output <- cbind(all_rows, comp_dat)
    
    # Add to cumulative output
    if(i == 1){
      output <- one_stand_output
    } else {
      output <- bind_rows(output, one_stand_output)
    }
  }
  
  # Reset row names
  rownames(output) <- 1:nrow(output)
  
  # Return output
  output
}
