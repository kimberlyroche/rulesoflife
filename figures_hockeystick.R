source("path_fix.R")

library(rulesoflife)
library(tidyverse)
library(ggraph)
library(igraph)

source("ggplot_fix.R")

# ------------------------------------------------------------------------------
#   Pull high confidence taxon pairs
# ------------------------------------------------------------------------------

data <- load_data(tax_level = "ASV")
plot_dir <- check_dir(c("output", "figures"))
rug_obj <- summarize_Sigmas(output_dir = "asv_days90_diet25_scale1")
scores_pieces <- apply(rug_obj$rug, 2, function(x) {
  calc_universality_score(x, return_pieces = TRUE)
})

p1_factor <- cut(scores_pieces[1,], breaks = c(0, 0.6, 0.7, 0.8, 0.9, Inf))
levels(p1_factor) <- c("<60%", "60-70%", "70-80%", "80-90%", ">90%")
plot_df <- data.frame(p1 = p1_factor,
                      p2 = scores_pieces[2,])
p <- ggplot(plot_df, aes(x = p1, y = p2)) +
  geom_violin() +
  labs(x = "proportion agreement in sign",
       y = "average absolute correlation")
ggsave(file.path(plot_dir, "piecewise_universality_scores.png"),
       p,
       units = "in",
       dpi = 100,
       height = 5,
       width = 6)
show(p)
