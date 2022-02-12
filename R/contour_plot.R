#' Create contour plot.
#' 
#' Creates a contour plot using values of a specified variable associated with
#' coordinates in the forest stand.
#' 
#' This function takes point estimates of some variable and interpolates a 
#' continuous contour plot across the forest stand. For some variables, such as
#' soil nutrient levels, the contour will remain accurate up to the stand
#' boundary. However, if the variable represents a summary of the neighborhood
#' around the points such as local tree density, any point whose neighborhood
#' overlaps the stand boundary will be an underestimate. In the latter case,
#' we recommend users correct for edge effects in their implementation of
#' \code{neighborhood_summary()}
#' 
#' @param grid_vals Data frame containing coordinates and values of the contour
#' variable for a series of points. Columns containing x and y coordinates must
#' be named \code{x_coord} and \code{y_coord} respectively but the name of the
#' contour variable column can be anything.
#' @param value Name of contour variable column as a string.
#' @return A contour plot will appear in the plotting window.
#' @examples
#' library(dplyr, warn.conflicts = F)
#' # Create coordinate grid
#' locations <- data.frame(
#'   loc_id = paste("A", 1:441, sep = ""),
#'   x_coord = rep(seq(0, 100, 5), times = 21),
#'   y_coord = rep(seq(0, 100, 5), each = 21))
#'   
#' # Calculate density for each point
#' nbhds <- neighborhoods(mapping, stands = "AB08", radius = 10,
#'                        coords = locations)
#' nbhd_summ <- neighborhood_summary(nbhds, "loc_id", 10, "angular",
#'                                   edge_correction = T, x_limit = 100,
#'                                   y_limit = 100)
#' 
#' # Format data for contour_plot
#' grid_data <- locations %>%
#'   left_join(nbhd_summ %>%
#'     select(loc_id, all_angle_sum),
#'     by = "loc_id") %>%
#'   select(-loc_id)
#' 
#' # Create contour plot
#' contour_plot(grid_data, "all_angle_sum")
#' @export
#' @importFrom magrittr %>%
#' @import dplyr

contour_plot <- function(grid_vals, value, edge_handling = "none", rad = NULL,
                         max_x = NULL, max_y = NULL){
  
  # Create contour plot
  plotly::plot_ly(x = grid_vals$x_coord,
                  y = grid_vals$y_coord,
                  z = grid_vals[, value],
                  type = "contour")
    
}
