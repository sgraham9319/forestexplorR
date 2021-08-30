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
#' overlaps the stand boundary will be an underestimate.
#' 
#' This function can handle these edge effects by calculating the proportion of
#' each point's neighborhood included in the data and multiplying their value
#' to estimate a complete neighborhood. This is achieved by specifying
#' \code{edge_handling = "multiply"} - see Examples.
#' 
#' @param grid_vals Data frame containing coordinates and values of the contour
#' variable for a series of points. Columns containing x and y coordinates must
#' be named \code{x_coord} and \code{y_coord} respectively but the name of the
#' contour variable column can be anything.
#' @param value Name of contour variable column as a string.
#' @param edge_handling String specifying whether values close to the plot edge
#' should be adjusted. Possible values are \code{"none"} (default) and
#' \code{"multiply"}. See Details.
#' @param rad Distance from edge of plot below which values should be adjusted.
#' Must be in same units as the coordinates in \code{grid_vals}
#' @param max_x Maximum possible x coordinate of forest plot.
#' @param max_y Maximum possible y coordinate of forest plot.
#' @return A contour plot will appear in the plotting window.
#' @examples
#' # Create coordinate grid
#' locations <- data.frame(
#'   loc_id = paste("A", 1:441, sep = ""),
#'   x_coord = rep(seq(0, 100, 5), times = 21),
#'   y_coord = rep(seq(0, 100, 5), each = 21))
#'   
#' # Calculate density for each point
#' nbhds <- neighborhoods(mapping, "AB08", 10, coords = locations)
#' nbhd_summ <- neighborhood_summary(nbhds, "loc_id", 10, "angular")
#' 
#' # Format data for contour_plot
#' grid_data <- locations %>%
#'   left_join(nbhd_summ %>%
#'     select(loc_id, all_angle_sum),
#'     by = "loc_id") %>%
#'   select(-loc_id)
#' 
#' # Create contour plot
#' contour_plot(grid_data, "all_angle_sum", edge_handling = "multiply",
#'   rad = 10, max_x = 100, max_y = 100)

contour_plot <- function(grid_vals, value, edge_handling = "none", rad = NULL,
                         max_x = NULL, max_y = NULL){
  
  # Check if edge handling requested
  if(edge_handling == "none"){
    
    # Create contour plot
    plotly::plot_ly(x = grid_vals$x_coord,
                    y = grid_vals$y_coord,
                    z = grid_vals[, value],
                    type = "contour")
    
  } else if(edge_handling == "multiply"){
    
    # Return error if necessary parameters not provided
    if(any(is.null(rad), is.null(max_x), is.null(max_y))){
      print("Error: arguments rad, max_x and max_y must be provided")
    } else {
      
      # Update edge values
      plot_data <- grid_vals %>%
        mutate(left_dist = x_coord,
               right_dist = max_x - x_coord,
               bot_dist = y_coord,
               top_dist = max_y - y_coord,
               tb_dist = if_else(bot_dist <= top_dist, bot_dist, top_dist),
               lr_dist = if_else(left_dist <= right_dist, left_dist, right_dist),
               tb_sag = if_else(tb_dist < rad, rad - tb_dist, 0),
               lr_sag = if_else(lr_dist < rad, rad - lr_dist, 0),
               tb_seg = seg_area(rad, tb_sag),
               lr_seg = seg_area(rad, lr_sag),
               tb_chord = chord_len(rad, tb_sag),
               lr_chord = chord_len(rad, lr_sag),
               tb_tri = if_else(tb_chord / 2 > lr_dist, tb_chord / 2 - lr_dist, 0),
               lr_tri = if_else(lr_chord / 2 > tb_dist, lr_chord / 2 - tb_dist, 0),
               missing_area = tb_seg + lr_seg - 0.5 * (tb_tri * lr_tri),
               prop_measured = 1 - (missing_area / (pi * rad ^ 2))) %>%
        select(x_coord, y_coord, all_of(value), prop_measured) %>%
        mutate(adj_value = get(value) / prop_measured)
      
      # Create plot
      plotly::plot_ly(plot_data,
                      x = plot_data$x_coord,
                      y = plot_data$y_coord,
                      z = plot_data$adj_value,
                      type = "contour")
    }
  }
}
