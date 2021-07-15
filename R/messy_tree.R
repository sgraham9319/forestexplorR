#' Uncleaned repeated size measurement data.
#' 
#' A dataset containing repeated size measurements for trees in 2 mapped
#' forest stands located in Mount Rainier National Park, WA, USA. Each row
#' represents a single measurement of a single tree. This is a subset of the
#' \code{tree} dataset with a variety of errors intentionally added in order to
#' test and demonstrate the \code{tree_check} function.
#' 
#' @format A data frame with 6501 rows and 6 variables:
#' \describe{
#'   \item{tree_id}{unique identification code of the tree}
#'   \item{stand_id}{name of plot in which tree is located}
#'   \item{tag}{unique tag number of the tree}
#'   \item{species}{species identity of tree as a four letter code}
#'   \item{year}{year when measurement was taken}
#'   \item{dbh}{diameter at breast height of the tree at the time of
#'   measurement, in cm}
#' }
#' @source Raw data (without the induced errors) are available from the Pacific
#' Northwest Permanent Sample Plot Network on request
#' \url{http://pnwpsp.forestry.oregonstate.edu/data}
"messy_tree"