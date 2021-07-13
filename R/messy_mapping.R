#' Uncleaned mapping data for 9166 trees.
#' 
#' A dataset containing x and y coordinates for trees in 15 mapped forest
#' stands located in Mount Rainier National Park, WA, USA. This is the same
#' dataset as \code{mapping} with a variety of errors intentionally added in
#' order to test and demonstrate the \code{mapping_check} function.
#' 
#' @format A data frame with 9166 rows and 8 variables:
#' \describe{
#'   \item{tree_id}{unique identification code of the tree}
#'   \item{stand_id}{name of plot in which tree is located}
#'   \item{tag}{unique tag number of the tree}
#'   \item{species}{species identity of tree as a four letter code}
#'   \item{year}{year when mapping data was collected}
#'   \item{dbh}{diameter at breast height of the tree at the time of
#'   measurement, in cm}
#'   \item{x_coord}{within-plot x coordinate of tree, in meters from origin}
#'   \item{y_coord}{within-plot y coordinate of tree, in meters from origin}
#' }
#' @source Raw data (without the induced errors) are available from the Pacific
#' Northwest Permanent Sample Plot Network on request
#' \url{http://pnwpsp.forestry.oregonstate.edu/data}
"messy_mapping"