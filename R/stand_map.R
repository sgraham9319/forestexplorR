#' Create a stand map.
#' 
#' Creates a map of the forest stand where each tree is a point, sized relative
#' to its most recent diameter at breast height (dbh) measurement and labeled
#' with its tag number. Trees for which dbh data is missing from \code{dbh_data}
#' will be plotted in blue instead of the standard red and assigned a size
#' equal to mean dbh in the plot.
#' 
#' The x and y limits of the stand map can be specified because if a large
#' stand (>50m in width or length) is mapped all at once the labels probably
#' will not be readable.
#' 
#' @param map_data Data frame of tree coordinates. Must include the columns:
#' \code{tree_id}, \code{x_coord}, \code{y_coord}, and \code{tag} (tag numbers
#' will be used to label trees on stand map).
#' @param dbh_data Data frame of tree measurement data. Must include the
#' columns: \code{tree_id}, \code{dbh}, and \code{year}
#' @param x_limit Numeric vector of length two indicating the x coordinate
#' range to be included.
#' @param y_limit Numeric vector of length two indicating the y coordinate
#' range to be included.
#' @return Stand map will be plotted in the plotting window.
#' @examples
#' library(dplyr)
#' # Isolate mapping and tree data for one stand
#' one_stand_map <- mapping %>%
#'   filter(stand_id == "AB08")
#' one_stand_tree <- tree %>%
#'   filter(stand_id == "AB08")
#'
#' # Create map of one quarter of the stand 
#' stand_map(one_stand_map, one_stand_tree, c(50, 100), c(50, 100))
#' @export
#' @importFrom magrittr %>%
#' @import dplyr

stand_map <- function(map_data, dbh_data, x_limit, y_limit){
  
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
  ggplot2::ggplot(plot_data) +
    ggplot2::geom_point(ggplot2::aes(x = x_coord, y = y_coord, size = dbh,
                                     fill = pt_color),
                        shape = 21, show.legend = F) + 
    ggplot2::scale_size_area(max_size = 5) +
    ggplot2::geom_text(ggplot2::aes(x = x_coord, y = y_coord, label = tag),
                       size = 3, vjust = -1,
                       position = ggplot2::position_jitter(width = 0.5,
                                                           height = 0.2)) +
    ggplot2::labs(x = "X coordinate", y = "Y coordinate") +
    ggplot2::theme_bw()
}