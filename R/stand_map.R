#' Create a stand map.
#' 
#' Creates a map of the forest stand where each tree is a point, sized relative
#' to its most recent diameter at breast height (dbh) measurement. Trees for
#' which dbh data is missing from \code{dbh_data} will be plotted in blue
#' instead of the standard red and assigned a size equal to mean dbh in the
#' plot. Points can optinally be labeled using a column included in
#' \code{map_data}, such as tag numbers of the trees.
#' 
#' The x and y limits of the stand map can be specified because if a large
#' stand (>50m in width or length) is mapped all at once the labels probably
#' will not be readable.
#' 
#' @param map_data Data frame of tree coordinates. Must include the columns:
#' \code{tree_id}, \code{x_coord}, and \code{y_coord}. If points (trees) are to
#' be labeled in the plot, \code{map_data} must also contain a column containing
#' the labels (e.g. tag numbers).
#' @param dbh_data Data frame of tree measurement data. Must include the
#' columns: \code{tree_id}, \code{dbh}, and \code{year}.
#' @param x_limit Numeric vector of length two indicating the x coordinate
#' range to be included.
#' @param y_limit Numeric vector of length two indicating the y coordinate
#' range to be included.
#' @param labels Name of column in \code{map_data} containing labels for trees
#' in the stem map, provided as a string. If \code{labels} argument is not 
#' provided, trees will not be labeled in the stem map.
#' @return Stand map will be plotted in the plotting window.
#' @examples
#' library(dplyr)
#' # Isolate mapping and tree data for one stand
#' one_stand_map <- mapping %>%
#'   filter(stand_id == "AB08")
#' one_stand_tree <- tree %>%
#'   filter(stand_id == "AB08")
#'
#' # Create map of one quarter of the stand using tag numbers as labels
#' stand_map(one_stand_map, one_stand_tree, c(50, 100), c(50, 100), "tag")
#' 
#' @export
#' @importFrom magrittr %>%
#' @import dplyr

stand_map <- function(map_data, dbh_data, x_limit, y_limit, labels = NULL){
  
  # Remove dbh column from mapping if present
  if("dbh" %in% names(map_data)){
    map_data <- map_data %>%
      select(-dbh)
  }
  
  # Subset to trees within x and y limits
  plot_data <- map_data %>%
    filter(x_coord >= x_limit[1] &
             x_coord <= x_limit[2] &
             y_coord >= y_limit[1] &
             y_coord <= y_limit[2])
  
  # Extract most recent dbh value for each tree
  most_recent_dbh <- dbh_data %>%
    group_by(tree_id) %>%
    arrange(year) %>%
    summarize(dbh = dbh[n()])
  
  # Attach dbh data
  plot_data <- plot_data %>%
    left_join(most_recent_dbh, by = "tree_id")
  
  # Add color column to plot_data, making NA dbh trees blue
  plot_data <- plot_data %>%
    mutate(pt_color = (if_else(is.na(dbh), "2", "1")))
  
  # Change any NA dbh values to mean dbh
  plot_data <- plot_data %>%
    mutate(dbh = if_else(is.na(dbh), mean(dbh, na.rm = T), dbh))
  
  # Create plot
  p <- ggplot2::ggplot(plot_data) +
    ggplot2::geom_point(ggplot2::aes(x = x_coord, y = y_coord, size = dbh,
                                     fill = pt_color),
                        shape = 21, show.legend = F) + 
    ggplot2::scale_size_area(max_size = 5) +
    ggplot2::labs(x = "X coordinate", y = "Y coordinate") +
    ggplot2::theme_bw()
  
  # Add point labels if requested
  if(!is.null(labels)){
    p <- p +
      ggplot2::geom_text(ggplot2::aes(x = x_coord, y = y_coord,
                                      label = get(labels)),
                         size = 3, vjust = -1,
                         position = ggplot2::position_jitter(width = 0.5,
                                                             height = 0.2))
  }
  
  # Output plot
  p
}