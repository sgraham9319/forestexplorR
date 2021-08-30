#' Plot size class distributions
#' 
#' Creates tree size class distribution plots for one or more stands.
#' 
#' The user-provided data frame \code{size_data} can contain trees from an 
#' unlimited number of stands but it is recommended that no more that 10 stands
#' be specified in the \code{stands} argument to maintain resolution in the 
#' resulting multi-panel plot. If \code{size_data} contains multiple size
#' measurements for a single tree, the \code{tree_id} and \code{year} columns
#' will be used to select the most recent size measurement for plotting.
#' 
#' @param size_data Data frame containing size information of trees. Must 
#' include the columns: \code{tree_id}, \code{stand_id}, \code{dbh}, and
#' \code{year} (i.e. year in which size measurement was taken).
#' @param stands Character vector of names of stands for which size class
#' distributions are desired. To maintain plot clarity it is recommended that
#' 10 stands or fewer are specified.
#' @param bin_size Integer representing desired width of size classes. Default
#' is 10 dbh units.
#' @return A single plot containing a panel for each of the specified
#' \code{stands}.
#' @examples
#' # Single stand
#' size_dist(size_data = tree, stands = "AB08", bin_size = 5)
#' 
#' # Multiple stands
#' size_dist(size_data = tree, stands = c("AM16", "AO03", "AR07"))
#' @export
#' @importFrom magrittr %>%
#' @import dplyr

size_dist <- function(size_data, stands, bin_size = 10){
  
  # Subset tree size data to focal stands
  size_data <- size_data %>%
    filter(stand_id %in% stands)
  
  # Extract most recent size measurement per tree
  size_data <- size_data %>%
    group_by(tree_id) %>%
    arrange(desc(year)) %>%
    filter(row_number() == 1)
  
  # Make plot
  if(length(stands) > 1){
    ggplot2::ggplot(data = size_data, aes(dbh)) +
      ggplot2::geom_histogram(aes(fill=stand_id),
                              breaks = seq(0, max(size_data$dbh) +
                                             bin_size, bin_size)) +
      ggplot2::labs(y = "# trees", x = "Diameter at breast height (cm)") +
      ggplot2::theme_classic() +
      ggplot2::facet_grid(rows = vars(size_data$stand_id)) +
      ggplot2::theme(legend.position = "none")
  } else {
    ggplot2::ggplot(data = size_data, aes(dbh)) +
      ggplot2::geom_histogram(breaks = seq(0, max(size_data$dbh) +
                                             bin_size, bin_size)) +
      ggplot2::labs(title = stands, y = "# trees",
                    x = "Diameter at breast height (cm)") +
      ggplot2::theme_classic()
  }
}
