#' Calculate between-census annual growth
#' 
#' Requires a data frame where separate measurements of the same tree appear
#' in different rows and returns a data frame containing annual growth between
#' each pair of consecutive measurements for each tree.
#' 
#' @param data A data frame containing repeated dbh measurements of trees. Must 
#' have a column named "tree_id" that contains unique tree ID values, a column
#' named "year" containing the year of the measurement, and a column named
#' "dbh" containing the dbh measurements. See built-in dataset \code{tree} for
#' an example.
#' @return A data frame containing between-census annual growth of each tree:
#' \describe{
#'   \item{tree_id}{unique identification code of the tree}
#'   \item{start_year}{first year of between-census period}
#'   \item{start_dbh}{diameter at breast height at beginning of between-census
#'   period, in cm}
#'   \item{end_year}{last year of between-census period}
#'   \item{end_dbh}{diameter at breast height at end of between-census
#'   period, in cm}
#'   \item{annual_growth}{annual growth rate as average yearly increase in 
#'   diameter at breast height, in cm/year}
#'   \item{...}{any other columns that appeared in \code{data}}
#' }
#' @examples
#' detailed_growth(tree)

detailed_growth <- function(data){
  
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
    filter(!is.na(year) & !is.na(dbh)) %>%
    group_by(tree_id) %>% 
    arrange(tree_id, year) %>%
    rename(start_dbh = dbh,
           start_year = year) %>%
    mutate(year_diff = c(diff(start_year), NA),
           dbh_diff = c(diff(start_dbh), NA),
           end_year = start_year + year_diff,
           end_dbh = start_dbh + dbh_diff,
           annual_growth = dbh_diff / year_diff) %>%
    select(-year_diff, -dbh_diff) %>%
    filter(!is.na(end_year))
  
  return(output)
  
}