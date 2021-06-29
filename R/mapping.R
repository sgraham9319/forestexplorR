#' Mapping information for 9155 trees.
#' 
#' A dataset containing x and y coordinates for all trees in 15 mapped forest
#' stands located in Mount Rainier National Park, WA, USA. These plots are 
#' part of the Pacific Northwest Permanent Sample Plot Network 
#' (\url{http://pnwpsp.forestry.oregonstate.edu/})
#' 
#' @format A data frame with 9155 rows and 8 variables:
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
#' @source Data are available from the Pacific Northwest Permanent Sample Plot
#' Network on request \url{http://pnwpsp.forestry.oregonstate.edu/data}
"mapping"