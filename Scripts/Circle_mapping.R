library(dplyr)

map_circ <- mapping %>%
  mutate(y_dist = y_coord - 50,
         x_dist = x_coord - 50,
         dist_from_center = sqrt((x_dist ^ 2) + (y_dist ^ 2))) %>%
  filter(dist_from_center <= 50) %>%
  mutate(bearing_rad = case_when(
    x_dist == 0 & y_dist == 0 ~ 0,
    x_dist >= 0 & y_dist >= 0 ~ atan(x_dist / y_dist),
    x_dist >= 0 & y_dist < 0 ~ pi - atan(x_dist / abs(y_dist)),
    x_dist < 0 & y_dist < 0 ~ pi + atan(abs(x_dist) / abs(y_dist)),
    x_dist < 0 & y_dist >= 0 ~ (2 * pi) - atan(abs(x_dist) / y_dist)),
    bearing_deg = bearing_rad / (pi / 180)
  )

# Need to extract neighborhoods using dist_from_center and bearing_deg

# Check it did what I expected
one_stand_map <- map_circ %>%
  filter(stand_id == "AB08")
one_stand_tree <- tree %>%
  filter(stand_id == "AB08")

# Create map of one quarter of the stand using tag numbers as labels
stand_map(one_stand_map, one_stand_tree, c(0, 100), c(0, 100), "tag")
