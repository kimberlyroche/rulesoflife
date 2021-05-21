source("path_fix.R")

library(rulesoflife)
library(tidyverse)

source("ggplot_fix.R")

data <- load_data(tax_level = "ASV")
output_dir <- "asv_days90_diet25_scale1"
rug_obj <- summarize_Sigmas(output_dir = output_dir)

universalities <- apply(rug_obj$rug, 2, calc_universality_score)
consensus_sign <- apply(rug_obj$rug, 2, calc_consensus_sign)

plot_df <- data.frame(universality = universalities,
                      universality_sign = consensus_sign,
                      tax_distance = rug_phylogenetic_distances(rug_obj,
                                                                data$taxonomy))
plot_df$universality_sign <- factor(plot_df$universality_sign)
levels(plot_df$universality_sign) <- c("negative", "neutral", "positive")

ggplot(plot_df, aes(x = tax_distance, y = universality, fill = universality_sign)) +
  geom_point(shape = 21, size = 2, alpha = 1) +
  scale_fill_manual(values = c("#5680e3", "#bbbbbb", "#e35656")) +
  facet_wrap(. ~ universality_sign) +
  labs(x = "phylogenetic distance",
       y = "universality score",
       fill = "Consensus sign")
ggsave(file.path(check_dir(c("output", "figures")), "phylogenetic_distance_vs_universality.png"),
       units = "in",
       dpi = 100,
       height = 4,
       width = 10)
