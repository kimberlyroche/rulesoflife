source("path_fix.R")

library(tidyverse)
library(rulesoflife)
library(driver)
library(RColorBrewer)
library(ggridges)
library(cowplot)

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
                      prop_agree = scores_piecewise[1,],
                      mcs = scores_piecewise[2,])
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

fit <- lm(mcs ~ mean1*mean2, data = plot_df)
cat(paste0("Mean CLR abundance of partner 1 is positively associated with median correlation strength: ",
           "beta = ", round(coef(summary(fit))[2,1], 3), ", p-value = ", round(coef(summary(fit))[2,4], 3), "\n"))
var_expl <- (var(plot_df$mcs) - var(fit$residuals))/var(plot_df$mcs)
cat(paste0("\tPercent variance explained: ", round(var_expl*100,3), "%\n"))

p1 <- ggplot(plot_df, aes(x = mean1, y = mcs)) +
  geom_point(size = 3, shape = 21, fill = "#999999") +
  geom_smooth(color = "black", alpha = 0.66, method = "lm") +
  theme_bw() +
  labs(x = "CLR ASV mean",
       y = "median correlation strength") +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))

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
  geom_density_ridges(quantile_lines = TRUE) +
  theme_bw() +
  labs(x = "difference in mean CLR abundance") +
  scale_fill_manual(values = brewer.pal(n = 4, name = "RdPu")[1:2]) +
  coord_cartesian(clip = "off") +
  scale_y_discrete(expand = expansion(add = c(0.2, 2))) +
  theme(axis.title.y = element_blank(),
        legend.position = "none",
        legend.background = element_rect(fill='transparent')) +
  guides(fill = guide_legend(reverse = TRUE)) +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14))

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

p3 <- ggplot(plot_df2, aes(x = mean, y = top, fill = partner)) +
  geom_density_ridges(alpha = 0.65, scale = 1.2, quantile_lines = TRUE) +
  theme_bw() +
  labs(x = "mean CLR abundance",
       fill = "Partner\nabundance") +
  scale_fill_manual(values = brewer.pal(n = 4, name = "RdPu")[c(4,1)]) +
  coord_cartesian(clip = "off") +
  scale_y_discrete(expand = expansion(add = c(0.2, 2))) +
  theme(axis.title.y = element_blank(),
        legend.position = c(0.8, 0.8),
        legend.background = element_rect(fill = 'transparent')) +
  guides(fill = guide_legend(reverse = TRUE)) +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14))

prow1 <- plot_grid(NULL, NULL, p1, NULL, ncol = 4,
                   scale = 0.95,
                   rel_widths = c(0.2, 0.1, 1, 0.2),
                   labels = c("", "A", "", ""),
                   label_size = 18,
                   label_x = 0.02,
                   label_y = 1.02)
prow2 <- plot_grid(p3, p2, ncol = 2,
                   scale = 0.95,
                   rel_widths = c(1, 0.9),
                   labels = c("B", "C"),
                   label_size = 19,
                   label_x = 0.02,
                   label_y = 1.02)
p <- plot_grid(prow1, NULL, prow2, ncol = 1,
               rel_heights = c(1, 0.05, 0.8)) +
  theme(plot.background = element_rect(fill = "white", color = "white"))

ggsave(file.path("output", "figures", "abundance.svg"),
       p,
       dpi = 100,
       units = "in",
       height = 9,
       width = 9)

# ------------------------------------------------------------------------------
#   Stats
# ------------------------------------------------------------------------------

fit <- wilcox.test(plot_df2 %>% filter(top == "Top 2.5%" & partner == "Lesser") %>% pull(mean),
                   plot_df2 %>% filter(top != "Top 2.5%" & partner == "Lesser") %>% pull(mean),
                   alternative = "two.sided",
                   conf.int = TRUE)
cat(paste0("Difference in location (panel B, lesser): ", round(fit$estimate, 3), " (p-value = ",
           round(fit$p.value, 5), ")\n"))

fit <- wilcox.test(plot_df2 %>% filter(top == "Top 2.5%" & partner == "Greater") %>% pull(mean),
                   plot_df2 %>% filter(top != "Top 2.5%" & partner == "Greater") %>% pull(mean),
                   alternative = "two.sided",
                   conf.int = TRUE)
cat(paste0("Difference in location (panel B, greater): ", round(fit$estimate, 3), " (p-value = ",
           round(fit$p.value, 5), ")\n"))

fit <- wilcox.test(plot_df %>% filter(top_factor == "Top 2.5%") %>% pull(delta),
                   plot_df %>% filter(top_factor != "Top 2.5%") %>% pull(delta),
                   alternative = "two.sided",
                   conf.int = TRUE)
cat(paste0("Difference in location (panel C): ", round(fit$estimate, 3), " (p-value = ",
           round(fit$p.value, 5), ")\n"))
