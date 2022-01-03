library(rulesoflife)
library(tidyverse)
library(cowplot)

# Lowest universality scores -- how do these look in terms of "pieces" of the
# score?

rug <- readRDS("output/rug_asv.rds")
scores <- apply(rug$rug, 2, calc_universality_score)
scores2 <- apply(rug$rug, 2, function(x) {
  calc_universality_score(x, return_pieces = TRUE)
})

lowest <- which(scores < 0.05)

p1 <- ggplot(data.frame(x = scores2[1,lowest]), aes(x = x)) +
  geom_histogram(color = "white") +
  theme_bw() +
  labs(x = "proportion agreement in sign")

p2 <- ggplot(data.frame(x = scores2[2,lowest]), aes(x = x)) +
  geom_histogram(color = "white") +
  theme_bw() +
  labs(x = "median correlation strength")

p <- plot_grid(p1, p2)

ggsave("output/figures/lowest_universality_distributions.png",
       p,
       dpi = 100,
       units = "in",
       height = 4,
       width = 7)

plot_df <- NULL
for(i in 1:length(lowest)) {
  plot_df <- rbind(plot_df,
                   data.frame(correlation = rug$rug[,lowest[i]],
                              idx = lowest[i]))
}

p <- ggplot(plot_df %>% filter(idx %in% lowest[sample(1:length(lowest), size = 50)]),
            aes(x = factor(idx), y = correlation)) +
  geom_violin() +
  theme_bw() +
  labs(x = "pair index")

ggsave("output/figures/lowest_universality_distributions2.png",
       p,
       dpi = 100,
       units = "in",
       height = 4,
       width = 12)
