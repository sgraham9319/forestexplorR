devtools::load_all()
library(dplyr)

# Extract mapping data for one stand to practice on
one_stand <- mapping %>%
  filter(stand_id == "AB08")

# User provides limits of x and y coordinates (don't assume lower limit is 0
# because they may want to look at only a portion of the plot)
x_limit <- c(50, 100)
y_limit <- c(50, 100)

# Exclude any points beyond x or y limits
one_stand <- one_stand %>%
  filter(x_coord >= x_limit[1] &
           x_coord <= x_limit[2] &
           y_coord >= y_limit[1] &
           y_coord <= y_limit[2])

# Define point size for each tree
one_stand <- one_stand %>%
  mutate(point_size = 0.5 + (dbh / max(dbh)))

# Adjust plotting window
par(mar = c(2, 2, 1, 1))

# Create plot
plot(NULL, ylim = y_limit, xlim = x_limit, yaxt = "n", xaxt = "n",
     ylab = "", xlab = "", yaxs = "i", xaxs = "i")
points(x = one_stand$x_coord, y = one_stand$y_coord,
       pch = 21, bg = "blue", cex = one_stand$point_size)

# Label trees with tag numbers
text(x = one_stand$x_coord, y = one_stand$y_coord + 1.2,
     labels = one_stand$tag, cex = 0.6)

# Label axes
mtext(text = x_limit[1], side = 1, line = 0.5, at = x_limit[1])
mtext(text = x_limit[2], side = 1, line = 0.5, at = x_limit[2])
mtext(text = y_limit[1], side = 2, line = 0.5, at = y_limit[1], las = 1)
mtext(text = y_limit[2], side = 2, line = 0.5, at = y_limit[2], las = 1)
mtext(text = "x-coordinate", side = 1, line = 0.6, cex = 1.2)
mtext(text = "y-coordinate", side = 2, line = 0.6, cex = 1.2)

# Reset plotting window
#dev.off()
