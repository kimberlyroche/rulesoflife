source("path_fix.R")

library(tidyverse)
library(rulesoflife)
library(driver)

# ------------------------------------------------------------------------------
#
#   Figure S6 - "bigraph": show there are similar trends in abundance of pairs
#               in the top 2.5% and bottom 97.5% of universality scores
#
# ------------------------------------------------------------------------------

# Pull top 2.5% most universal pairs
rug_asv <- summarize_Sigmas(output_dir = "asv_days90_diet25_scale1")
scores <- apply(rug_asv$rug, 2, calc_universality_score)

# Get average CLR abundance for all taxa
data <- load_data(tax_level = "ASV")
clr.counts <- clr_array(data$counts + 0.5, parts = 1)
clr.means <- rowMeans(clr.counts)

plot_df <- data.frame(pair_idx = 1:ncol(rug_asv$rug),
                      tax_idx1 = rug_asv$tax_idx1,
                      tax_idx2 = rug_asv$tax_idx2,
                      score = scores)
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
#   Version 1: plot a) distributions of differences in mean CLR and b)
#              distributions of mean CLR for each partner in pair
#
# ------------------------------------------------------------------------------

plot_df <- plot_df %>%
  mutate(delta = mean2 - mean1)

plot_df$top_factor <- factor(plot_df$top)
levels(plot_df$top_factor) <- c("Bottom 97.5%", "Top 2.5%")

p1 <- ggplot(plot_df, aes(x = delta, y = top_factor, fill = top_factor)) +
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

p2 <- ggplot(plot_df2, aes(x = mean, y = top, fill = partner)) +
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

p <- plot_grid(p2, p1, ncol = 2,
               scale = 0.95,
               rel_widths = c(1, 0.9),
               labels = c("a", "b"),
               label_size = 22,
               label_x = 0.02,
               label_y = 1.02) +
  theme(plot.background = element_rect(fill = "white", color = "white"))

ggsave(file.path("output", "figures", "S6.png"),
       p,
       dpi = 100,
       units = "in",
       height = 4,
       width = 8)

# ------------------------------------------------------------------------------
#
#   Version 2: bigraph
#
# ------------------------------------------------------------------------------

# Subset huge number of pairs in the bottom 97.5%
bottom_df_subset <- bottom_df %>%
  mutate(key = sample(1:nrow(bottom_df))) %>%
  arrange(key) %>%
  slice(1:1000)

point_df <- data.frame(x = "lower",
                       y = top_df %>% pull(mean1),
                       top = TRUE)
point_df <- rbind(point_df,
                  data.frame(x = "upper",
                             y = top_df %>% pull(mean2),
                             top = TRUE))
point_df <- rbind(point_df,
                  data.frame(x = "lower",
                             y = bottom_df_subset %>% pull(mean1),
                             top = FALSE))
point_df <- rbind(point_df,
                  data.frame(x = "upper",
                             y = bottom_df_subset %>% pull(mean2),
                             top = FALSE))

line_df <- data.frame(x = "lower",
                      xend = "upper",
                      y = top_df %>% pull(mean1),
                      yend = top_df %>% pull(mean2),
                      top = TRUE)
line_df <- rbind(line_df,
                 data.frame(x = "lower",
                            xend = "upper",
                            y = bottom_df_subset %>% pull(mean1),
                            yend = bottom_df_subset %>% pull(mean2),
                            top = FALSE))

line_alpha <- 0.25
line_size <- 0.75

p1 <- ggplot() +
  geom_segment(data = line_df %>% filter(!top),
               aes(x = factor(x), xend = factor(xend), y = y, yend = yend),
               alpha = line_alpha/2,
               size = line_size) +
  geom_point(data = point_df %>% filter(!top),
             aes(x = x, y = y)) +
  theme_bw() +
  labs(x = "\nmember of taxon pair",
       y = "mean empirical CLR abundance\n") +
  coord_cartesian(clip = "off") +
  scale_x_discrete(expand = expansion(add = c(0.2, 0.2))) +
  ylim(c(-2.5, 6.1))

p2 <- ggplot() +
  geom_segment(data = line_df %>% filter(top),
               aes(x = factor(x), xend = factor(xend), y = y, yend = yend),
               alpha = line_alpha,
               size = line_size) +
  geom_point(data = point_df %>% filter(top),
             aes(x = x, y = y)) +
  theme_bw() +
  labs(x = "\nmember of taxon pair",
       y = "mean empirical CLR abundance\n") +
  coord_cartesian(clip = "off") +
  scale_x_discrete(expand = expansion(add = c(0.2, 0.2))) +
  ylim(c(-2.5, 6.1))

p <- plot_grid(p1, p2, ncol = 2,
          scale = 0.9,
          rel_widths = c(1, 0.9),
          labels = c("a", "b"),
          label_size = 22,
          label_x = 0.0,
          label_y = 1.02) +
  theme(plot.background = element_rect(fill = "white", color = "white"))

ggsave(file.path("output", "figures", "S6_alt.png"),
       p,
       dpi = 100,
       units = "in",
       height = 5,
       width = 7)
