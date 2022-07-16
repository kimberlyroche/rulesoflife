source("path_fix.R")

library(tidyverse)
library(rulesoflife)
library(driver)
library(ggridges)
library(cowplot)

# ------------------------------------------------------------------------------
#
#   Supplemental Figure S5 - "hockey stick" plots for family and phylum levels
#                            and stacked density plots of universality scores
#                            for phyla and families
#
# ------------------------------------------------------------------------------

source("thresholds.R")

rug_fam <- summarize_Sigmas(output_dir = "fam_days90_diet25_scale1")

fam_scores_pieces <- apply(rug_fam$rug, 2, function(x) calc_universality_score(x, return_pieces = TRUE))
fam_score_df <- data.frame(x = fam_scores_pieces[1,], y = fam_scores_pieces[2,])

# Code below copied from `figures_main.R`
fam_score_df$sign <- apply(rug_fam$rug, 2, calc_consensus_sign)
fam_score_df$sign <- factor(fam_score_df$sign, levels = c(1, -1))
levels(fam_score_df$sign) <- c("positive", "negative")
# 2) Percent significant observations for this taxon pair
fam_score_df$signif <- apply(rug_fam$rug, 2, function(x) {
  sum(x < thresholds %>% filter(type == "ASV") %>% pull(lower) | x > thresholds %>% filter(type == "ASV") %>% pull(upper))/length(x)
})

rug_phy <- summarize_Sigmas(output_dir = "phy_days90_diet25_scale1")

phy_scores_pieces <- apply(rug_phy$rug, 2, function(x) calc_universality_score(x, return_pieces = TRUE))
phy_score_df <- data.frame(x = phy_scores_pieces[1,], y = phy_scores_pieces[2,])

# Code below copied from `figures_main.R`
phy_score_df$sign <- apply(rug_phy$rug, 2, calc_consensus_sign)
phy_score_df$sign <- factor(phy_score_df$sign, levels = c(1, -1))
levels(phy_score_df$sign) <- c("positive", "negative")
# 2) Percent significant observations for this taxon pair
phy_score_df$signif <- apply(rug_phy$rug, 2, function(x) {
  sum(x < thresholds %>% filter(type == "ASV") %>% pull(lower) | x > thresholds %>% filter(type == "ASV") %>% pull(upper))/length(x)
})

score_df <- rbind(cbind(fam_score_df, type = "Family/order/class"),
                  cbind(phy_score_df, type = "Phylum"))

p1 <- ggplot(score_df %>% filter(sign != 0)) +
  geom_point(mapping = aes(x = x, y = y, fill = signif, color = sign),
             size = 2,
             shape = 21,
             stroke = 1) +
  scale_fill_distiller(palette = "PuRd", direction = 1,
                       guide = guide_colorbar(frame.colour = "black",
                                              ticks.colour = "black")) +
  scale_color_manual(values = c("#000000", "#888888")) +
  facet_wrap(. ~ type) +
  theme_bw() +
  ylim(c(0,1)) +
  theme(legend.title = element_text(margin = margin(b = 5))) +
  labs(x = "proportion shared sign",
       y = "median association strength",
       fill = "Proportion\nsignificant\nobservations",
       color = "Consensus\ncorrelation sign")

# ------------------------------------------------------------------------------
#   ggridges stacked histograms of universality scores at the levels higher
#   than ASV
# ------------------------------------------------------------------------------

plot_df <- data.frame(score = apply(rug_phy$rug, 2, calc_universality_score),
                      type = "Phylum")
plot_df <- rbind(plot_df,
                 data.frame(score = apply(rug_fam$rug, 2, calc_universality_score),
                            type = "Family/order/class"))

plot_df2 <- thresholds_scores %>%
  filter(type != "ASV")
plot_df2$type <- factor(plot_df2$type)

p2 <- ggplot(plot_df, aes(x = score, y = type, fill = type)) +
  geom_density_ridges(stat = "binline", bins = 20, scale = 0.9, draw_baseline = TRUE) +
  geom_segment(data = plot_df2,
               aes(x = x0, xend = x0, y = as.numeric(type), yend = as.numeric(type) + .9),
               color = "black",
               size = 1,
               linetype = "dashed") +
  theme_bw() +
  labs(x = "universality score",
       y = "") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme_ridges(grid = FALSE, center_axis_labels = TRUE) +
  theme(legend.position = "none") +
  scale_fill_brewer(palette = "Blues")

p2_padded <- plot_grid(NULL, p2, ncol = 2,
                       rel_widths = c(0.03, 1),
                       labels = c("", "B"),
                       label_size = 18,
                       label_x = 0.05)

p <- plot_grid(p1, NULL, p2_padded, ncol = 3,
               rel_widths = c(2, -0.15, 1.25),
               labels = c("A", "", ""),
               label_size = 18,
               # label_x = -0.02,
               scale = 0.95)

ggsave(file.path("output", "figures", "S5.svg"),
       p,
       dpi = 100,
       units = "in",
       height = 4,
       width = 12)
