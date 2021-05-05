# Testing out color palette options.

color_scale <- 255

red_start <- c(240, 7, 7)/color_scale
red_end <- c(145, 54, 54)/color_scale
orange_start <- c(227, 125, 16)/color_scale
orange_end <- c(145, 116, 54)/color_scale
green_start <- c(68, 194, 10)/color_scale
green_end <- c(67, 122, 42)/color_scale
blue_start <- c(90, 141, 242)/color_scale
blue_end <- c(54, 83, 145)/color_scale
purple_start <- c(200, 90, 242)/color_scale
purple_end <- c(116, 58, 138)/color_scale

palette_length <- 6
palette_no <- 5
red_palette <- colorRampPalette(c(rgb(red_start[1],
                                      red_start[2],
                                      red_start[3], 1),
                                  rgb(red_end[1],
                                      red_end[2],
                                      red_end[3], 1)))(palette_length)
orange_palette <- colorRampPalette(c(rgb(orange_start[1],
                                         orange_start[2],
                                         orange_start[3], 1),
                                     rgb(orange_end[1],
                                         orange_end[2],
                                         orange_end[3], 1)))(palette_length)
green_palette <- colorRampPalette(c(rgb(green_start[1],
                                        green_start[2],
                                        green_start[3], 1),
                                    rgb(green_end[1],
                                        green_end[2],
                                        green_end[3], 1)))(palette_length)
blue_palette <- colorRampPalette(c(rgb(blue_start[1],
                                       blue_start[2],
                                       blue_start[3], 1),
                                   rgb(blue_end[1],
                                       blue_end[2],
                                       blue_end[3], 1)))(palette_length)
purple_palette <- colorRampPalette(c(rgb(purple_start[1],
                                         purple_start[2],
                                         purple_start[3], 1),
                                     rgb(purple_end[1],
                                         purple_end[2],
                                         purple_end[3], 1)))(palette_length)
palette <- c(red_palette, orange_palette, green_palette, blue_palette, purple_palette)

plot_df <- data.frame(x = rep(1:palette_length, palette_no),
                      y = rep(1:palette_no, each = palette_length),
                      type = factor(1:length(palette)))
ggplot(plot_df, aes(x = x, y = y, color = type)) +
  geom_point(size = 10) +
  scale_color_manual(values = palette)
