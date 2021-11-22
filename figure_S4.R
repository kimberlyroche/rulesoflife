source("path_fix.R")

library(tidyverse)
library(rulesoflife)
library(driver)
library(cowplot)

# ------------------------------------------------------------------------------
#
#   Supplemental Figure S4 - "hockey stick" plots for family and phylum levels
#
# ------------------------------------------------------------------------------

# These are hard-coded. See `figures_supplemental.R` for the calculation of
# "important"/significant correlations from permuted data.
thresholds <- data.frame(type = factor(c("Phylum", "Family", "ASV")),
                         lower = c(-0.303, -0.256, -0.263),
                         upper = c(0.149, 0.207, 0.254))

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

score_df <- rbind(cbind(fam_score_df, type = "Family"),
                  cbind(phy_score_df, type = "Phylum"))

p <- ggplot(score_df %>% filter(sign != 0)) +
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

ggsave("output/figures/S4.svg",
       p,
       dpi = 100,
       units = "in",
       height = 4,
       width = 8)
