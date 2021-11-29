source("path_fix.R")

library(tidyverse)
library(rulesoflife)
library(RColorBrewer)
library(cowplot)

# ------------------------------------------------------------------------------
#
#   Figure 3 - three "rug" plots at different taxonomic levels plus breakout
#              plots of negatively and positively correlated highly "universal"
#              taxa
#
# ------------------------------------------------------------------------------

rugs <- list(a = summarize_Sigmas(output_dir = "phy_days90_diet25_scale1"),
             b = summarize_Sigmas(output_dir = "fam_days90_diet25_scale1"),
             c = summarize_Sigmas(output_dir = "asv_days90_diet25_scale1"))

plots <- list()
legend <- NULL
for(rtype in names(rugs)) {
  order_obj <- order_rug_row_pedigree(rugs[[rtype]])
  row_order <- order_obj$order

  # Compute column order
  rug <- rugs[[rtype]]$rug
  row_labels <- rugs[[rtype]]$hosts
  cluster_obj <- order_obj$hc
  canonical_col_order <- order(colMeans(rug))
  canonical_row_order <- row_order
  rug <- rug[canonical_row_order,canonical_col_order]
  row_labels <- row_labels[canonical_row_order]

  rug <- cbind(1:nrow(rug), rug)
  colnames(rug) <- c("host", paste0(1:(ncol(rug)-1)))
  rug <- pivot_longer(as.data.frame(rug), !host, names_to = "pair", values_to = "correlation")
  rug$pair <- as.numeric(rug$pair)

  p <- ggplot(rug, aes(x = pair, y = host)) +
    geom_raster(aes(fill = correlation)) +
    scale_fill_gradient2(low = "navy", mid = "white", high = "red", midpoint = 0,
                         guide = guide_colorbar(frame.colour = "black",
                                                ticks.colour = "black"))

  if(is.null(row_labels)) {
    row_labels <- 1:length(rug$host)
  }

  title <- "ASVs"
  xlabel <- "ASV-ASV pairs"
  if(rtype == "a") {
    title <- "Phyla"
    xlabel <- "phylum-phylum pairs"
  } else if(rtype == "b") {
    title <- "Families"
    xlabel <- "family-family pairs"
  }

  row_breaks <- 1:length(row_labels)

  # Subset the row labels for readability
  select_vec <- rep(c(FALSE, TRUE, FALSE, FALSE), 56/4)

  p <- p +
    scale_x_continuous(expand = c(0, 0)) +
    labs(fill = "Correlation",
         x = xlabel,
         y = "",
         title = title) +
    scale_y_continuous(breaks = row_breaks[select_vec],
                       labels = row_labels[select_vec],
                       expand = c(0, 0)) +
    theme(axis.text.y = element_text(size = rel(0.95)),
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 14, face = "plain"),
          plot.title = element_text(size = 20, hjust = 0.5),
          legend.title = element_text(size = 14),
          # legend.position = "bottom",
          plot.margin = margin(t = 20, r = 10, b = 10, l = 0),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank())
  # legend.text = element_text(margin = margin(b = -20)))
  if(rtype == "a") {
    p <- p +
      theme(plot.margin = margin(t = 20, r = 10, b = 10, l = 0))
  } else {
    p <- p +
      theme(plot.margin = margin(t = 20, r = 10, b = 10, l = 10))
  }
  if(rtype == "c") {
    legend <- get_legend(p)
  }
  p <- p +
    theme(legend.position = "none")
  if(rtype == "a") {
    dhc <- as.dendrogram(cluster_obj)
    ddata <- dendro_data(dhc, type = "rectangle")
    pd <- ggplot(segment(ddata)) +
      geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
      coord_flip() +
      scale_y_reverse() +
      theme_nothing() +
      theme(plot.margin = margin(t = 0, r = -10, b = 0, l = 10))
    plots[[1]] <- pd
  }
  plots[[length(plots)+1]] <- p
}

# Pad the dendrogram
p1 <- plot_grid(nullGrob(), plots[[1]], nullGrob(), ncol = 1, rel_heights = c(0.08, 1, 0.05))
p2 <- plot_grid(p1, plots[[2]],
                ncol = 2,
                rel_widths = c(0.1, 0.9))
p3 <- plot_grid(p2, plots[[3]], plots[[4]],
                ncol = 3,
                scale = 1,
                rel_widths = c(1, 1.2, 1.2),
                labels = c("a", "b", "c"),
                label_size = 24,
                label_x = 0.02,
                label_y = 0.99)
p4 <- plot_grid(p3, legend, ncol = 2, rel_widths = c(1, 0.1))
ggsave("output/figures/F3.svg",
       p4,
       dpi = 100,
       units = "in",
       height = 6,
       width = 18)

# ------------------------------------------------------------------------------
#   Breakout plots of positively and negative correlated taxa
# ------------------------------------------------------------------------------

pred_obj <- get_predictions_host_list(host_list = c("VEX"),
                                      output_dir = "asv_days90_diet25_scale1",
                                      metadata = data$metadata)

l_idx <- which(rugs$a$hosts == "VEX")
eta_df <- gather_array(pred_obj$predictions$VEX$Eta, val, coord, sample, iteration)

pos_pair <- eta_df %>%
  filter(coord %in% c(2,3)) %>%
  group_by(coord, sample) %>%
  summarize(p2.5 = quantile(val, probs = c(0.025)),
            p25 = quantile(val, probs = c(0.25)),
            mean = mean(val),
            p75 = quantile(val, probs = c(0.75)),
            p97.5 = quantile(val, probs = c(0.975)))

neg_pair <- eta_df %>%
  filter(coord %in% c(31,114)) %>%
  group_by(coord, sample) %>%
  summarize(p2.5 = quantile(val, probs = c(0.025)),
            p25 = quantile(val, probs = c(0.25)),
            mean = mean(val),
            p75 = quantile(val, probs = c(0.75)),
            p97.5 = quantile(val, probs = c(0.975)))

alpha <- 0.6

p <- ggplot() +
  geom_ribbon(data = pos_pair %>% filter(coord == 2),
              mapping = aes(x = sample, ymin = p25, ymax = p75),
              fill = "#fdbf6f",
              alpha = alpha) +
  geom_line(data = pos_pair %>% filter(coord == 2),
            mapping = aes(x = sample, y = mean),
            color = "#ff7f00",
            size = 1.5,
            alpha = alpha) +
  geom_ribbon(data = pos_pair %>% filter(coord == 3),
              mapping = aes(x = sample, ymin = p25, ymax = p75),
              fill = "#a6cee3",
              alpha = alpha) +
  geom_line(data = pos_pair %>% filter(coord == 3),
            mapping = aes(x = sample, y = mean),
            color = "#1f78b4",
            size = 1.5,
            alpha = alpha) +
  theme_bw() +
  labs(x = "days from first sample",
       y = "CLR abundance")

ggsave("output/figures/F3_positive_pair.svg",
       p,
       dpi = 100,
       units = "in",
       height = 2,
       width = 6)

cat(paste0("Median correlation across hosts for this pair: ",
           round(median(rugs$c$rug[,which(rugs$c$tax_idx1 == 2 & rugs$c$tax_idx2 == 3)]), 3)))

p <- ggplot() +
  geom_ribbon(data = neg_pair %>% filter(coord == 31),
              mapping = aes(x = sample, ymin = p25, ymax = p75),
              fill = "#fdbf6f",
              alpha = alpha) +
  geom_line(data = neg_pair %>% filter(coord == 31),
            mapping = aes(x = sample, y = mean),
            color = "#ff7f00",
            size = 1.5,
            alpha = alpha) +
  geom_ribbon(data = neg_pair %>% filter(coord == 114),
              mapping = aes(x = sample, ymin = p25, ymax = p75),
              fill = "#a6cee3",
              alpha = alpha) +
  geom_line(data = neg_pair %>% filter(coord == 114),
            mapping = aes(x = sample, y = mean),
            color = "#1f78b4",
            size = 1.5,
            alpha = alpha) +
  theme_bw() +
  labs(x = "days from first sample",
       y = "CLR abundance")
ggsave("output/figures/F3_negative_pair.svg",
       p,
       dpi = 100,
       units = "in",
       height = 2,
       width = 6)

cat(paste0("Median correlation across hosts for this pair: ",
           round(median(rugs$c$rug[,which(rugs$c$tax_idx1 == 31 & rugs$c$tax_idx2 == 114)]), 3)))

