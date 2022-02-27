#' Check tree measurement data
#' 
#' Checks a tree measurements dataset for missing data and formatting that could 
#' lead to errors when applying other functions in this package. Outputs a
#' table of the trees with data issues and a table summarizing the number of
#' trees in the mapping dataset that have data issues.
#' 
#' The data issues checked for are: presence of required columns, single tree
#' ids referring to multiple trees, trees having no associated mapping data,
#' missing dbh, stand id, species, or measurement year. The provided
#' \code{tree_data} is also checked for the presence of a \code{mort} column
#' containing mortality data; if this column is found, a check for missing
#' mortality data is also performed. This function does not check for misspelled
#' stand ids or species, which should be checked independently. The common issue
#' of negative growth rates resulting from measurement error are not checked
#' here, but are checked by \code{growth_summary}.
#' 
#' Tree ids indicated to have data issues according to this function are not 
#' necessarily unusable. For instance, missing year data could be inferred from
#' knowledge on when certain stands were measured. The tree ids with data issues
#' should therefore be investigated further rather than being excluded from 
#' further analyses right away.
#' 
#' @param tree_data Data frame containing tree measurement data where each row
#' represents a single measurement of a single tree. Should contain the columns
#' \code{tree_id}, \code{stand_id}, \code{species}, \code{year} and \code{dbh}.
#' If a \code{mort} column (i.e. mortality data used for the 
#' \code{mortality_model} function) it will be checked for missing values but no
#' warning will be produced if this column is absent. Any additional columns
#' will be ignored by this function.
#' @param map_data Data frame containing tree mapping data. Should contain the
#' columns \code{tree_id}, \code{stand_id}, \code{species}, \code{x_coord},
#' and \code{y_coord}. Any additional columns will be ignored by this function.
#' Note that this function does not check the mapping data, which should first 
#' be checked with the function \code{mapping_check}.
#' @return A list containing two elements:
#' \itemize{
#'   \item{\code{problem_trees} is a data frame containing the tree ids found
#'   to have data issues and a description of the issue}
#'   \item{\code{issue_summary} is a data frame that shows the number and 
#'   percentage of trees with at least one issue and with each of the specific
#'   issues}
#' }
#' @examples
#' tree_check_test <- tree_check(messy_tree, mapping)
#' @export
#' @importFrom magrittr %>%
#' @import dplyr

tree_check <- function(tree_data, map_data){
  
  # Check required columns are present
  if(sum(c("tree_id", "stand_id", "species", "year", "dbh") %in%
         names(tree_data)) != 5){
    print("Warning: one or more of the required columns are missing")
    print("Required columns: tree_id, stand_id, species, year, dbh")
  } else {
    
    # Identify trees with missing data
    problem_ids <- bind_rows(
      tree_data %>%
        group_by(tree_id) %>%
        summarize(num_stands = length(unique(stand_id)),
                  num_species = length(unique(species))) %>%
        filter(num_stands > 1 | num_species > 1) %>%
        mutate(issue = "refers to multiple trees") %>%
        select(tree_id, issue),
      tree_data %>%
        filter(tree_id %in% map_data$tree_id == F) %>%
        mutate(issue = "no mapping data") %>%
        select(tree_id, issue),
      tree_data %>%
        filter(is.na(dbh)) %>%
        mutate(issue = "missing dbh") %>%
        select(tree_id, issue),
      tree_data %>%
        filter(is.na(stand_id)) %>%
        mutate(issue = "missing stand id") %>%
        select(tree_id, issue),
      tree_data %>%
        filter(is.na(species)) %>%
        mutate(issue = "missing species id") %>%
        select(tree_id, issue),
      tree_data %>%
        filter(is.na(year)) %>%
        mutate(issue = "missing measurement year") %>%
        select(tree_id, issue))
    
    # If mort column is present, check for missing values
    if("mort" %in% names(tree_data)){
      problem_ids <- bind_rows(
        problem_ids,
        tree_data %>%
          filter(is.na(mort)) %>%
          mutate(issue = "missing mortality record") %>%
          select(tree_id, issue)
      )
    }
    
    # Calculate number of tree ids in tree data
    inds <- length(unique(tree_data$tree_id))
    
    # Get number and percent of trees in each issue category
    issue_summary <- problem_ids %>%
      group_by(issue) %>%
      summarize(count = n()) %>%
      mutate(pct = 100 * (count / inds))
    
    # Get number and percent of trees with at least one issue
    unique_issue_trees <- data.frame(
      issue = "at least one issue",
      count = length(unique(problem_ids$tree_id)),
      pct = 100 * (length(unique(problem_ids$tree_id)) / inds))
    
    # Combine into single output table
    issue_output <- bind_rows(unique_issue_trees, issue_summary)
    
    # Print summary message
    if(nrow(problem_ids) == 0){
      print("No potential formatting problems detected")
    } else {
      print("Potential formatting problems detected: please review output and correct errors or remove problem trees if necessary before continuing analysis")
    }
    
    # Return both output tables as a list
    return(list(problem_trees = problem_ids, issue_summary = issue_output))
  }
}
