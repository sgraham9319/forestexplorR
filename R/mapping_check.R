#' Check tree mapping data
#' 
#' Checks a tree mapping dataset for missing data and formatting that could 
#' lead to errors when applying other functions in this package. Outputs a
#' table of the trees with data issues and a table summarizing the number of
#' trees in the mapping dataset that have data issues.
#' 
#' The data issues checked for are: presence of required columns, duplicated
#' tree ids, missing x or y coordinates, x or y coordinates outside the
#' expected range, missing stand id or species information. This function does
#' not check for misspelled stand ids or species, which should be checked 
#' independently.
#' 
#' @param map_data Data frame containing tree mapping data. Should contain the
#' columns \code{tree_id}, \code{stand_id}, \code{species}, \code{x_coord}, and
#' \code{y_coord}. Any additional columns will be ignored by this function.
#' @param max_x Maximum expected x coordinate (i.e. should be 100 if the stands
#' are 100 x 100 m).
#' @param max_y Maximum expected y coordinate.
#' @return A list containing two elements:
#' \itemize{
#'   \item{\code{problem_trees} is a data frame containing the rows of
#'   \code{map_data} that contain data issues. An additional column describes
#'   the identified issue}
#'   \item{\code{issue_summary} is a data frame that shows the number and 
#'   percentage of trees with at least one issue and with each of the specific
#'   issues}
#' }
#' @examples
#' map_check_test <- mapping_check(messy_mapping, 100, 100)
#' @export
#' @importFrom magrittr %>%
#' @import dplyr

mapping_check <- function(map_data, max_x, max_y){
  
  # Check required columns are present
  if(sum(c("tree_id", "stand_id", "species", "x_coord", "y_coord") %in%
         names(map_data)) != 5){
    print("Warning: one or more of the required columns are missing")
    print("Required columns: tree_id, stand_id, species, x_coord, y_coord")
  } else {
    
    # Identify trees with missing or nonsensical data
    problem_trees <- bind_rows(
      map_data %>%
        filter(duplicated(tree_id)) %>%
        mutate(issue = "duplicated tree_id"),
      map_data %>%
        filter(is.na(map_data$y_coord) | is.na(map_data$x_coord)) %>%
        mutate(issue = "missing coordinates"),
      map_data %>%
        filter(!is.na(y_coord) & !is.na(x_coord)) %>%
        filter(y_coord < 0 | y_coord > max_y | x_coord < 0 | x_coord > max_x) %>%
        mutate(issue = "coordinates out of range"),
      map_data %>%
        filter(is.na(stand_id)) %>%
        mutate(issue = "missing stand_id"),
      map_data %>%
        filter(is.na(species)) %>%
        mutate(issue = "missing species"))
    
    # Get number and percent of trees in each issue category
    issue_summary <- problem_trees %>%
      group_by(issue) %>%
      summarize(count = n()) %>%
      mutate(pct = 100 * (count / nrow(map_data)))
    
    # Get number and percent of trees with at least one issue
    unique_issue_trees <- data.frame(
      issue = "at least one issue",
      count = length(unique(problem_trees$tree_id)),
      pct = 100 * (length(unique(problem_trees$tree_id)) / nrow(map_data)))
    
    # Combine into single output table
    issue_output <- bind_rows(unique_issue_trees, issue_summary)
    
    # Return both output tables as a list
    return(list(problem_trees = problem_trees, issue_summary = issue_output))
    
  }
}
