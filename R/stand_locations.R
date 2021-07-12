#' Locations of 15 forest stands.
#' 
#' A data frame containing latitudes, longitudes, and y-azimuths for 15 mapped
#' forest stands located in Mount Rainier National Park, WA, USA. These forest
#' stands are part of the Pacific Northwest Permanent Sample Plot Network 
#' (\url{http://pnwpsp.forestry.oregonstate.edu/})
#' 
#' @format A data frame with 15 rows and 4 variables:
#' \describe{
#'   \item{stand_id}{unique identification code of each stand}
#'   \item{y_azim}{y-azimuth of the local coordinate system of each stand}
#'   \item{latitude}{latitude of stand origin (i.e. local coordinates 0, 0) in
#'   decimal degrees}
#'   \item{longitude}{longitude of stand origin (i.e. local coordinates 0, 0) in
#'   decimal degrees}
#' }
#' @source Data are available from the Pacific Northwest Permanent Sample Plot
#' Network on request \url{http://pnwpsp.forestry.oregonstate.edu/data}
"stand_locations"