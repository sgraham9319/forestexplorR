#' Create neighborhoods
#'
#' Takes a data frame containing mapping information of trees in one or more
#' stands and returns neighborhood information for all trees or a set of 
#' user-provided coordinates.
#' 
#' This function requires dbh information for each tree in \code{mapping} in 
#' order to return the dbh of each neighboring tree in each neighborhood. If
#' \code{mapping} contains no dbh data or it is desired to use different dbh
#' data (mapping datasets usually contain dbh at initial measurement, which
#' could be outdated), the argument \code{dbh_data} must be specified.
#' 
#' This function returns a neighborhoods object, which is a data frame where
#' each focal tree or user-provided set of coordinates appears on multiple 
#' lines with each line containing information on one of the trees in its
#' neighborhood. Specifically, each line contains the ID number
#' (\code{id_comp}), species identity (\code{sps_comp}), size (\code{dbh_comp} 
#' and \code{abh_comp}), and distance from the neighborhood center (\code{prox})
#' of the neighbor tree. The neighborhoods object can be passed into a number
#' of other functions in ForestPlotR, including \code{neighborhood_summary} and 
#' \code{site_by_species}.
#' 
#' @param mapping Data frame containing tree coordinates.
#' @param dbh_data Data frame containing a single dbh measurement for each tree
#' id in \code{mapping}. Must have columns \code{tree_id} and \code{dbh} - any 
#' additional columns will be ignored by this function.
#' @param stands Vector of names of stands for which neighborhoods are desired.
#' @param radius Numeric vector describing neighborhood radius in meters.
#' @param coords Data frame containing coordinates for which neighborhoods are
#' desired. Data frame should contain three columns in the following order -
#' location ids, x-coordinates, y-coordinates. Column names are unimportant.
#' @return Neighborhood information for all focal trees in \code{mapping} or,
#' if \code{coords} is provided, for all locations defined by coordinates.
#' @examples
#' # Create neighborhoods for trees in mapping
#' nbhds <- neighborhoods(mapping, stands = c("AB08", "PP17"), radius = 10)
#' 
#' # Create neighborhoods for trees in mapping using more up-to-date dbh data
#' library(dplyr)
#' tree_dbhs <- tree %>%
#'   filter(stand_id == "AB08") %>%
#'   group_by(tree_id) %>%
#'   arrange(year) %>%
#'   summarize(dbh = dbh[n()])
#' nbhds <- neighborhoods(mapping, dbh_data = tree_dbhs, stands = "AB08",
#'                        radius = 10)
#' 
#' # Create neighborhoods for user-provided coordinates
#' locations <- data.frame(
#' loc_id = paste("A", 1:81, sep = ""),
#' x_coord = rep(seq(10, 90, 10), times = 9),
#' y_coord = rep(seq(10, 90, 10), each = 9))
#' nbhds <- neighborhoods(mapping, stands = "AB08", radius = 10,
#'                        coords = locations)
#' @export
#' @importFrom magrittr %>%
#' @import dplyr

neighborhoods <- function(mapping, dbh_data = NULL, stands = "all", radius,
                          coords = NULL) {
  
  # Return error message if no dbh data provided
  if("dbh" %in% names(mapping) == F & "dbh" %in% names(dbh_data) == F){
    stop("No dbh data provided: a column of dbh data must be provided in
    mapping or dbh_data", call. = F)
  }
  
  # Specify the dbh data to use
  if("dbh" %in% names(dbh_data)){
    
    # If dbh in mapping and dbh_data, remove from mapping
    if("dbh" %in% names(mapping)){
      mapping <- mapping %>%
        select(-dbh)
    }
    
    # Bring dbh from dbh_data into mapping
    mapping <- mapping %>%
      left_join(dbh_data %>% 
                  select(tree_id, dbh),
                by = "tree_id")
  }
  
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
    
    # If coordinates provided, add them to one_stand
    if(!is.null(coords)){
      names(coords) <- c("tree_id", "x_coord", "y_coord")
      one_stand <- one_stand %>%
        bind_rows(coords)
      one_stand$stand_id <- rep(stand_list[i], times = nrow(one_stand))
    }
    
    # Create distance matrix
    dist_mat <- as.matrix(dist(one_stand[, c("x_coord", "y_coord")]))
    
    # Create competitor species matrix
    sps_mat <- matrix(one_stand$species, nrow = nrow(one_stand),
                      ncol = nrow(one_stand))
    
    # Create competitor dbh matrix
    dbh_mat <- matrix(one_stand$dbh, nrow = nrow(one_stand),
                      ncol = nrow(one_stand))
    
    # Create competitor abh matrix
    abh_mat <- matrix(one_stand$abh, nrow = nrow(one_stand),
                      ncol = nrow(one_stand))
    
    # Create competitor tree_id matrix
    id_mat <- matrix(one_stand$tree_id, nrow = nrow(one_stand),
                     ncol = nrow(one_stand))
    
    # Find values in distance matrix greater than radius or equal to zero -
    # these elements represent pairs of trees that are not in each others
    # neighborhood or are the same tree paired with itself
    non_comp <- which(dist_mat > radius | dist_mat == 0)
    
    # Convert values of these non-competing pairs to NA in all matrices
    dist_mat[non_comp] <- NA
    sps_mat[non_comp] <- NA
    dbh_mat[non_comp] <- NA
    abh_mat[non_comp] <- NA
    id_mat[non_comp] <- NA
    
    # If coordinates provided, remove focal tree neighborhoods from matrices
    if(!is.null(coords)){
      coord_rows <- (nrow(dist_mat) - (nrow(coords) - 1)):nrow(dist_mat)
      dist_mat <- dist_mat[-coord_rows, coord_rows]
      sps_mat <- sps_mat[-coord_rows, coord_rows]
      dbh_mat <- dbh_mat[-coord_rows, coord_rows]
      abh_mat <- abh_mat[-coord_rows, coord_rows]
      id_mat <- id_mat[-coord_rows, coord_rows]
    }
    
    # Create distance, species, dbh, abh and id vectors, excluding NAs
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
    if(is.null(coords)){
      all_rows <- one_stand[rep(1:nrow(one_stand), repeats), ]
    } else {
      all_rows <- one_stand[rep(coord_rows[1]:nrow(one_stand), repeats),
                            c("tree_id", "stand_id", "x_coord", "y_coord")]
    }
    
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
  
  # Change column names if coordinates provided
  if(!is.null(coords)){
    names(output)[1] <- "loc_id"
  }
  
  # Return output
  output
  
}
