#' Calculate annual growth
#'
#' Requires a data frame where separate measurements of the same tree appear
#' in different rows and returns a data frame with a single row for each tree
#' containing annual growth measurements.
#' 
#' @param data A dataframe containing repeated measurements of growth. Must 
#' have a column named "tree_id" that contains unique tree ID values, a column
#' named "year" containing the year of the measurement, and a column named
#' "dbh" containing the dbh measurements.See built-in dataset \code{tree} for
#' an example.
#' @return A dataframe containing growth rate measurements for each tree:
#' \describe{
#'   \item{tree_id}{unique identification code of the tree}
#'   \item{stand_id}{name of stand in which tree is located - only appears if 
#'   input dataframe contains a \code{stand_id} column}
#'   \item{species}{species identity of tree as a four letter code}
#'   \item{first_record}{year of first size measurement of the tree}
#'   \item{last_record}{year of last size measurement of the tree}
#'   \item{begin_size}{first size measurement as diameter at breast height,
#'   in cm}
#'   \item{final_size}{last size measurement as diameter at breast height,
#'   in cm}
#'   \item{mean_size}{average size across measurements as diameter at breast
#'   height, in cm}
#'   \item{midpoint_size}{average of maximum and minimum size across
#'   measurements as diameter at breast height, in cm}
#'   \item{annual_growth}{annual growth rate as average yearly increase in 
#'   diameter at breast height, in cm/year}
#'   \item{annual_bai}{annual basal area increment as average yearly increase
#'   in area at breast height, in cm^2/year}
#'   \item{size_corr_growth}{square root of \code{annual_growth} divided by
#'   \code{begin_size} to give somewhat normally distributed growth rates}
#'   \item{size_corr_growth_basal}{square root of \code{annual_bai} divided by
#'   \code{begin_size} to give somewhat normally distributed growth rates}
#' }
#' @examples
#' growth_summary(tree)
#' @export
#' @importFrom magrittr %>%
#' @import dplyr

growth_summary <- function(data){
  
  # Check input data format
  if("tree_id" %in% names(data) == F){
    stop("required tree_id column is missing from input")
  }
  if("year" %in% names(data) == F){
    stop("required year column is missing from input")
  }
  if("dbh" %in% names(data) == F){
    stop("required dbh column is missing from input")
  }
  
  # Calculate annual growth
  output <- data %>%
    group_by(tree_id) %>% 
    arrange(year) %>%
    summarize(
      stand_id = ifelse("stand_id" %in% names(data), stand_id[1], NA),
      species = species[1],
      first_record = year[1],
      last_record = year[n()],
      begin_size = dbh[1],
      final_size = dbh[n()],
      mean_size = mean(dbh),
      midpoint_size = (min(dbh) + max(dbh)) / 2
    ) %>% 
    mutate(
      annual_growth = if_else(
        first_record == last_record,
        NA_real_,
        (final_size - begin_size) / (last_record - first_record)),
      annual_bai = if_else(
        first_record == last_record,
        NA_real_,
        (circ_area(final_size / 2) - circ_area(begin_size / 2)) /
          (last_record - first_record)))
  
  # Remove stand_id column if all NA
  if(all(is.na(output$stand_id))){
    output <- output %>%
      select(-stand_id)
  }
  
  # Calculate size corrected growth (radial and basal) for trees where defined
  output$size_corr_growth <- NA
  output$size_corr_growth_basal <- NA
  valid_rows <- which(output$annual_growth >= 0)
  output$size_corr_growth[valid_rows] <- 
    sqrt(output$annual_growth[valid_rows] /
           output$begin_size[valid_rows])
  output$size_corr_growth_basal[valid_rows] <-
    sqrt(output$annual_bai[valid_rows] /
           output$begin_size[valid_rows])
  
  # Give warning if some trees showed negative growth
  neg_grow <- length(which(output$annual_growth < 0))
  if(neg_grow > 0){
    print(paste("Warning:", neg_grow, "trees exhibited negative annual growth"))
  }
  
  return(output)
  
}