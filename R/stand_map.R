#' Create a stand map.
#' 
#' Creates a map of the forest stand where each tree is a point, sized relative
#' to its diameter at breast height (dbh) and labeled with its tag number.
#' 
#' The x and y limits of the stand map can be specified because if a large
#' stand (>50m in width or length) is mapped all at once the labels probably
#' will not be readable.
#' 
#' @param map_data Data frame of tree coordinates. Must include the columns:
#' \code{x_coord}, \code{y_coord}, \code{dbh}, and \code{tag} (tag numbers will
#' be used to label trees on stand map).
#' @param x_limit Numeric vector of length two indicating the x coordinate
#' range to be included.
#' @param y_limit Numeric vector of length two indicating the y coordinate
#' range to be included.
#' @return Stand map will be plotted in the plotting window.
#' @examples 
#' # Isolate mapping data for one stand
#' one_stand <- mapping %>%
#'   filter(stand_id == "AB08")
#' 
#' # Create stand map
#' stand_map(one_stand, c(50, 100), c(50, 100))
#' @export
#' @importFrom magrittr %>%
#' @import dplyr
#' 

stand_map <- function(map_data, x_limit, y_limit){
  
  # Subset to trees within x and y limits
  plot_data <- map_data %>%
    filter(x_coord >= x_limit[1] &
             x_coord <= x_limit[2] &
             y_coord >= y_limit[1] &
             y_coord <= y_limit[2])
  
  # Create plot
  ggplot2::ggplot(plot_data) +
    ggplot2::geom_point(ggplot2::aes(x = x_coord, y = y_coord, size = dbh),
                        shape = 21, fill = "red", show.legend = F) + 
    ggplot2::scale_size_area(max_size = 5) +
    ggplot2::geom_text(ggplot2::aes(x = x_coord, y = y_coord, label = tag),
                       size = 3, vjust = -1,
                       position = ggplot2::position_jitter(width = 0.5,
                                                           height = 0.2)) +
    ggplot2::labs(x = "X coordinate", y = "Y coordinate") +
    ggplot2::theme_bw()
}
