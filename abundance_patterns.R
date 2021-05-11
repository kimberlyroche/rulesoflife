library(rulesoflife)
library(driver)
library(ggplot2)
library(dplyr)

source("ggplot_fix.R")

source("plot_abundance_distros.R")

# Pull data
rug_obj <- summarize_Sigmas(output_dir = "asv_days90_diet25_scale1")
scores <- apply(rug_obj$rug, 2, calc_universality_score)
consensus_signs <- apply(rug_obj$rug, 2, calc_consensus_sign)

# Get mean CLR abundance of all taxa
data <- load_data(tax_level = "ASV")
clr.counts <- clr_array(data$counts + 0.5, parts = 1)
clr.means <- rowMeans(clr.counts)

# col_order <- order(colMeans(rug_obj$rug))
# plot_rug(rug_obj$rug, canonical_col_order = col_order)

k <- 100
# pairs <- list(bottom = col_order[1:k],
#               middle = col_order[ceiling(ncol(rug_obj$rug)/2-k/2):floor(ncol(rug_obj$rug)/2+k/2)],
#               top = col_order[(length(col_order)-k+1):length(col_order)])

temp_df <- data.frame(index = 1:length(scores),
                      score = scores,
                      sign = consensus_signs)
pairs <- list(bottom = temp_df %>% filter(sign < 0) %>% arrange(desc(score)) %>% slice(1:k) %>% pull(index),
              top = temp_df %>% filter(sign > 0) %>% arrange(desc(score)) %>% slice(1:k) %>% pull(index),
              middle = temp_df %>% arrange(score) %>% slice(1:k) %>% pull(index))
# There's a sort of keystone-ness going on in the top and bottom pairs (in terms
# of correlation): the same taxa occur frequently in these pairs.
# Optionally downsample the number of unique taxa in the "middle" samples; since
# there are more unique taxa here than in the extreme negative or positive pairs
# the plots look very dense without downsampling.
# pairs[["middle"]] <- sample(pairs[["middle"]], size = 50)

for(use in names(pairs)) {
  plot_abundance_distros(rug_obj, clr.means, pairs[[use]], paste0(use, "_", k))
}

# Quick bipartite plots
get_paired_abundance <- function(pairs, clr.means) {
  pair1 <- rug_obj$tax_idx1[pairs]
  pair2 <- rug_obj$tax_idx2[pairs]

  clr1 <- clr.means[pair1]
  clr2 <- clr.means[pair2]

  k <- length(pairs)
  max_pair <- sapply(1:k, function(x) {
    max(clr1[x], clr2[x])
  })
  min_pair <- sapply(1:k, function(x) {
    min(clr1[x], clr2[x])
  })
  return(list(min = min_pair, max = max_pair))
}

neutral_pairs <- get_paired_abundance(pairs$middle, clr.means)
negative_pairs <- get_paired_abundance(pairs$bottom, clr.means)
positive_pairs <- get_paired_abundance(pairs$top, clr.means)

line_df <- data.frame(x = rep(1, k), xend = rep(2, k),
                      y = neutral_pairs$min, yend = neutral_pairs$max,
                      type = rep("neutral pairs", k))
line_df <- rbind(line_df,
                 data.frame(x = rep(1, k), xend = rep(2, k),
                            y = negative_pairs$min, yend = negative_pairs$max,
                            type = rep("top 100 positive pairs", k)))
line_df <- rbind(line_df,
                 data.frame(x = rep(1, k), xend = rep(2, k),
                            y = positive_pairs$min, yend = positive_pairs$max,
                            type = rep("top 100 negative pairs", k)))
line_df$type <- factor(line_df$type, levels = c("neutral pairs", "top 100 positive pairs", "top 100 negative pairs"))
point_df <- data.frame(x = c(line_df$x, line_df$xend),
                       y = c(line_df$y, line_df$yend),
                       type = line_df$type)

ggplot(line_df) +
  geom_segment(aes(x = factor(x), xend = factor(xend), y = y, yend = yend, color = type), alpha = 0.5) +
  geom_point(data = point_df, aes(x = factor(x), y = y, color = type), alpha = 0.5) +
  facet_wrap(. ~ type) +
  scale_x_discrete(breaks = c(1, 2),
                   labels = c("less abundant", "more abundant")) +
  scale_color_manual(values = c("#888888", "red", "navy")) +
  theme_bw() +
  theme(axis.title.x=element_blank()) +
  labs(x = NULL,
       y = "mean CLR abundance",
       color = "Pair type")
plot_dir <- check_dir(c("output", "figures"))
ggsave(file.path(plot_dir, "abundance_pairs_bipartite.png"),
       units = "in",
       dpi = 100,
       height = 5,
       width = 9)
