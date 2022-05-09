source("path_fix.R")

library(tidyverse)
library(rulesoflife)
library(driver)
library(RColorBrewer)

# ------------------------------------------------------------------------------
#
#   Supplemental Figure 7 - abundance is not a driver of correlation;
#               a) distributions of differences in mean CLR in top 2.5% most
#               universal vs. all other pairs and b) distributions of mean CLR
#               for each partner in same pairs
#
# ------------------------------------------------------------------------------

# Pull top 2.5% most universal pairs
rug_asv <- summarize_Sigmas(output_dir = "asv_days90_diet25_scale1")
scores <- apply(rug_asv$rug, 2, calc_universality_score)
scores_piecewise <- apply(rug_asv$rug, 2, function(x) calc_universality_score(x = x, return_pieces = TRUE))

# Get average CLR abundance for all taxa
data <- load_data(tax_level = "ASV")
clr.counts <- clr_array(data$counts + 0.5, parts = 1)
clr.means <- rowMeans(clr.counts)

plot_df <- data.frame(pair_idx = 1:ncol(rug_asv$rug),
                      tax_idx1 = rug_asv$tax_idx1,
                      tax_idx2 = rug_asv$tax_idx2,
                      score = scores,
                      mad = scores_piecewise[1,],
                      prop_agree = scores_piecewise[2,])
plot_df$mean1 <- clr.means[plot_df$tax_idx1]
plot_df$mean2 <- clr.means[plot_df$tax_idx2]

threshold <- quantile(scores, probs = c(0.975))
plot_df$top <- sapply(scores, function(x) x > threshold)

# Rearrange taxa 1/2 to be lower-higher
for(i in 1:nrow(plot_df)) {
  if(plot_df$mean2[i] < plot_df$mean1[i]) {
    t1 <- plot_df$tax_idx1[i]
    plot_df$tax_idx1[i] <- plot_df$tax_idx2[i]
    plot_df$tax_idx2[i] <- t1
    m1 <- plot_df$mean1[i]
    plot_df$mean1[i] <- plot_df$mean2[i]
    plot_df$mean2[i] <- m1
  }
}

# ------------------------------------------------------------------------------
#
#   Statistical test(s): Does average abundance of either partner microbe
#                        predict universality score?
#
# ------------------------------------------------------------------------------

fit <- lm(log(score) ~ mean1*mean2, data = plot_df)
coef(summary(fit))

p0 <- ggplot(plot_df, aes(x = mean1, y = log(score))) +
  geom_point() +
  theme_bw() +
  labs(x = "CLR ASV mean",
       y = "log universality score") #+
  # xlim(c(-2.5, 3.75))

# ------------------------------------------------------------------------------
#
#   Visual confirmation: plot a) distributions of differences in mean CLR and b)
#              distributions of mean CLR for each partner in pair
#
# ------------------------------------------------------------------------------

plot_df <- plot_df %>%
  mutate(delta = mean2 - mean1)

plot_df$top_factor <- factor(plot_df$top)
levels(plot_df$top_factor) <- c("Bottom 97.5%", "Top 2.5%")

p2 <- ggplot(plot_df, aes(x = delta, y = top_factor, fill = top_factor)) +
  geom_density_ridges() +
  theme_bw() +
  labs(x = "difference in mean CLR abundance") +
  scale_fill_manual(values = brewer.pal(n = 4, name = "RdPu")[1:2]) +
  coord_cartesian(clip = "off") +
  scale_y_discrete(expand = expansion(add = c(0.2, 2))) +
  theme(axis.title.y = element_blank(),
        legend.position = "none",
        legend.background = element_rect(fill='transparent')) +
  guides(fill = guide_legend(reverse = TRUE))

top_df <- plot_df %>%
  filter(top)
bottom_df <- plot_df %>%
  filter(!top)

plot_df2 <- rbind(data.frame(mean = top_df$mean1,
                             partner = "Lesser",
                             top = TRUE),
                  data.frame(mean = top_df$mean2,
                             partner = "Greater",
                             top = TRUE))
plot_df2 <- rbind(plot_df2,
                  rbind(data.frame(mean = bottom_df$mean1,
                                   partner = "Lesser",
                                   top = FALSE),
                        data.frame(mean = bottom_df$mean2,
                                   partner = "Greater",
                                   top = FALSE)))

plot_df2$top <- factor(plot_df2$top)
levels(plot_df2$top) <- c("Bottom 97.5%", "Top 2.5%")

p1 <- ggplot(plot_df2, aes(x = mean, y = top, fill = partner)) +
  geom_density_ridges(alpha = 0.5, scale = 1.2) +
  theme_bw() +
  labs(x = "mean CLR abundance",
       fill = "Partner\nabundance") +
  scale_fill_manual(values = brewer.pal(n = 4, name = "RdPu")[c(4,1)]) +
  coord_cartesian(clip = "off") +
  scale_y_discrete(expand = expansion(add = c(0.2, 2))) +
  theme(axis.title.y = element_blank(),
        legend.position = c(0.8, 0.8),
        legend.background = element_rect(fill='transparent')) +
  guides(fill = guide_legend(reverse = TRUE))

# p <- plot_grid(p1, p2, ncol = 2,
#                scale = 0.95,
#                rel_widths = c(1, 0.9),
#                labels = c("A", "B"),
#                label_size = 22,
#                label_x = 0.02,
#                label_y = 1.02) +
#   theme(plot.background = element_rect(fill = "white", color = "white"))

prow1 <- plot_grid(NULL, NULL, p0, NULL, ncol = 4,
                   scale = 0.95,
                   rel_widths = c(0.1, 0.1, 1, 0.2),
                   labels = c("", "A", "", ""),
                   label_size = 22,
                   label_x = 0.02,
                   label_y = 1.02)
prow2 <- plot_grid(p1, p2, ncol = 2,
                   scale = 0.95,
                   rel_widths = c(1, 0.9),
                   labels = c("B", "C"),
                   label_size = 22,
                   label_x = 0.02,
                   label_y = 1.02)
p <- plot_grid(prow1, NULL, prow2, ncol = 1,
               rel_heights = c(1, 0.05, 0.9)) +
  theme(plot.background = element_rect(fill = "white", color = "white"))

ggsave(file.path("output", "figures", "S7.svg"),
       p,
       dpi = 100,
       units = "in",
       height = 9,
       width = 8)
