source("path_fix.R")

library(tidyverse)
library(rulesoflife)
library(RColorBrewer)
library(cowplot)
library(grid)
library(ggridges)
library(ggraph)
library(igraph)
library(scales)

# source("thresholds.R")

rug_fam <- summarize_Sigmas(output_dir = "fam_days90_diet25_scale1")$rug
rug_phy <- summarize_Sigmas(output_dir = "phy_days90_diet25_scale1")$rug

# ------------------------------------------------------------------------------
#   "Hockey stick" plot
# ------------------------------------------------------------------------------

plot_hockeystick <- function(rug) {
  scores_pieces <- apply(rug, 2, function(x) calc_universality_score(x, return_pieces = TRUE))
  scores <- apply(scores_pieces, 2, function(x) x[1]*x[2])
  score_df <- data.frame(x = scores_pieces[1,], y = scores_pieces[2,])

  # Check number of pairs with 100% sign agreement
  sign_agree <- which(scores_pieces[1,] == 1)
  sign_disagree <- which(scores_pieces[1,] == 0.5)
  cat(paste0(length(sign_agree), " pairs have 100% sign agreement\n"))
  cat(paste0(round(median(abs(rug[,sign_agree])), 3), " median assocation strength (100% agree)\n"))
  cat(paste0(round(median(abs(rug[,sign_disagree])), 3), " median assocation strength (<100% agree)\n"))

  # Additional labelings
  # 1) Consensus sign
  score_df$sign <- apply(rug, 2, calc_consensus_sign)
  score_df$sign <- factor(score_df$sign, levels = c(1, -1))
  levels(score_df$sign) <- c("positive", "negative")
  # 2) Percent significant observations for this taxon pair
  score_df$signif <- apply(rug, 2, function(x) {
    sum(x < thresholds %>% filter(type == "ASV") %>% pull(lower) | x > thresholds %>% filter(type == "ASV") %>% pull(upper))/length(x)
  })

  sign_palette <- c("red", "#0047AB")
  names(sign_palette) <- c(1, -1)

  dummy <- ggplot(data.frame(x = 1:2, y = 1:2, sign = c(1, -1)),
                  aes(x = x, y = y, fill = factor(sign))) +
    geom_point(size = 3, shape = 21) +
    scale_fill_manual(values = sign_palette, labels = c("positive", "negative")) +
    theme_bw() +
    labs(fill = "Consensus sign") +
    guides(alpha = "none",
           color = "none")

  legend <- get_legend(dummy)

  p <- ggplot() +
    geom_point(data = score_df %>% filter(sign == "negative"),
               mapping = aes(x = x, y = y),
               size = 2,
               shape = 21,
               fill = muted("navy")) +
    geom_point(data = score_df %>% filter(sign == "positive"),
               mapping = aes(x = x, y = y),
               size = 2,
               shape = 21,
               fill = "red") +
    theme_bw() +
    ylim(c(0,1)) +
    theme(axis.text.x = element_text(size = 11),
          axis.title.x = element_text(size = 13),
          axis.text.y = element_text(size = 11),
          axis.title.y = element_text(size = 13)) +
    labs(x = "proportion shared sign",
         y = "median correlation strength",
         color = "Consensus\ncorrelation sign")

  list(p = p + theme(legend.position = "none"), legend = legend)
}

p_fam_pieces <- plot_hockeystick(rug_fam)
p_phy_pieces <- plot_hockeystick(rug_phy)

row1 <- plot_grid(p_fam_pieces$p, p_phy_pieces$p, legend,
                  ncol = 3,
                  rel_widths = c(1, 1, 0.3),
                  labels = c("A", "B"),
                  label_y = 1.02,
                  label_size = 16,
                  scale = 0.95)

# ------------------------------------------------------------------------------
#   ggridges stacked histograms of universality scores at the ASV level
# ------------------------------------------------------------------------------

plot_univ_hist <- function(rug) {
  ggplot(data.frame(score = apply(rug, 2, calc_universality_score)),
         aes(x = score, y = 1)) +
    geom_density_ridges(stat = "binline", bins = 20, scale = 0.9, draw_baseline = TRUE) +
    theme_bw() +
    labs(x = "universality score",
         y = "") +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    theme_ridges(grid = FALSE, center_axis_labels = TRUE) +
    theme(legend.position = "none") +
    xlim(c(-0.05, 0.7)) +
    scale_fill_brewer(palette = "Blues") +
    theme(axis.text = element_text(size = 9),
          axis.title = element_text(size = 11))
}

row2 <- plot_grid(plot_univ_hist(rug_fam), plot_univ_hist(rug_phy), NULL,
                  ncol = 3,
                  rel_widths = c(1, 1, 0.2),
                  labels = c("C", "D"),
                  label_y = 1.02,
                  label_size = 16,
                  scale = 0.95)

p <- plot_grid(row1, row2,
               ncol = 1,
               rel_heights = c(1, 0.8))

ggsave(file.path("output", "figures", "other-hockeysticks.svg"),
       p,
       dpi = 100,
       units = "in",
       height = 7,
       width = 10)
