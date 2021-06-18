#' Repeated size measurements for 9155 trees.
#' 
#' A dataset containing repeated size and health condition measurements for
#' 9155 trees in 15 mapped forest stands located in Mount Rainier National
#' Park, WA, USA. Each row represents a single measurement of a single tree.
#' These plots are part of the Pacific Northwest Permanent Sample Plot Network 
#' (\url{http://pnwpsp.forestry.oregonstate.edu/})
#' 
#' @format A data frame with 60114 rows and 25 variables:
#' \describe{
#'   \item{tree_id}{unique identification code of the tree}
#'   \item{stand_id}{name of plot in which tree is located}
#'   \item{plot}{subplot in which tree is located}
#'   \item{tag}{unique tag number of the tree}
#'   \item{species}{species identity of tree as a four letter code}
#'   \item{year}{year when mapping data was collected}
#'   \item{tree_status}{code defining tree health condition}
#'   \item{dbh}{diameter at breast height of the tree at the time of
#'   measurement, in cm}
#'   \item{x_coord}{within-plot x coordinate of tree, in meters from origin}
#'   \item{y_coord}{within-plot y coordinate of tree, in meters from origin}
#'   \item{size_cat}{size category of tree, regular if dbh > 15 but small if
#'   dbh between 5 and 15}
#' }
#' @source Data are available from the Pacific Northwest Permanent Sample Plot
#' Network on request \url{http://pnwpsp.forestry.oregonstate.edu/data}
"tree"