#' Repeated size measurements for 9155 trees.
#' 
#' A dataset containing repeated size measurements for 9155 trees in 15 mapped
#' forest stands located in Mount Rainier National Park, WA, USA. Each row
#' represents a single measurement of a single tree. These plots are part of
#' the Pacific Northwest Permanent Sample Plot Network 
#' (\url{http://pnwpsp.forestry.oregonstate.edu/})
#' 
#' @format A data frame with 60114 rows and 6 variables:
#' \describe{
#'   \item{tree_id}{unique identification code of the tree}
#'   \item{stand_id}{name of plot in which tree is located}
#'   \item{tag}{unique tag number of the tree}
#'   \item{species}{species identity of tree as a four letter code}
#'   \item{year}{year when measurement was taken}
#'   \item{dbh}{diameter at breast height of the tree at the time of
#'   measurement, in cm}
#' }
#' @source Data are available from the Pacific Northwest Permanent Sample Plot
#' Network on request \url{http://pnwpsp.forestry.oregonstate.edu/data}
"tree"