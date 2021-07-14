library(plotly)

# Create coordinate grid
locations <- data.frame(
  loc_id = paste("A", 1:441, sep = ""),
  x_coord = rep(seq(0, 100, 5), times = 21),
  y_coord = rep(seq(0, 100, 5), each = 21))

# Create neighborhoods
nbhds <- neighborhoods(mapping, stands = "AB08", radius = 10,
                       coords = locations)

# Create neighborhood summaries
nbhd_summ <- neighborhood_summary(nbhds, "loc_id", 10, "angular")

# Define gridded values dataset
grid_data <- locations %>%
  left_join(nbhd_summ %>%
              select(loc_id, all_angle_sum),
            by = "loc_id") %>%
  select(-loc_id)

# Define function
contour_plot <- function(grid_vals, value, edge_handling = "none", rad = NULL,
                         max_x = NULL, max_y = NULL){
  
  # Check if edge handling requested
  if(edge_handling == "none"){
    
    # Create contour plot
    plot_ly(x = grid_vals$x_coord,
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
      plot_ly(plot_data,
              x = plot_data$x_coord,
              y = plot_data$y_coord,
              z = plot_data$adj_value,
              type = "contour")
    }
  }
}

# Examples
contour_plot(grid_data, "all_angle_sum")
contour_plot(grid_data, "all_angle_sum", "multiply", 10, 100, 100)
