source("path_fix.R")

library(driver)
library(tidyverse)
library(rulesoflife)
library(RColorBrewer)
library(cowplot)

# ------------------------------------------------------------------------------
#
#   Figure 2 - three "rug" plots at different taxonomic levels plus breakout
#              plots of negatively and positively correlated highly "universal"
#              taxa
#
# ------------------------------------------------------------------------------

rugs <- list(ASV = summarize_Sigmas(output_dir = "asv_days90_diet25_scale1"),
             family = summarize_Sigmas(output_dir = "fam_days90_diet25_scale1"),
             phylum = summarize_Sigmas(output_dir = "phy_days90_diet25_scale1"))

plots <- list()
legend <- NULL
for(rtype in names(rugs)) {
  # Order rows by similarity of baseline composition
  data <- load_data(tax_level = rtype)
  md <- data$metadata
  H <- length(unique(md$sname))
  D <- nrow(data$counts)
  baselines <- matrix(NA, H, D)
  for(i in 1:length(hosts)) {
    # Average of ALR samples for this host is their baseline
    host_samples <- clr_array(data$counts[,which(md$sname == hosts[i])] + 0.5,
                              parts = 1)
    baselines[i,] <- apply(host_samples, 1, mean)
  }

  baseline_distances <- dist(baselines)
  row_order <- hclust(baseline_distances)$order

  # Compute column order
  rug <- rugs[[rtype]]$rug
  canonical_col_order <- order(colMeans(rug))
  canonical_row_order <- row_order
  rug <- rug[canonical_row_order,canonical_col_order]

  rug <- cbind(1:nrow(rug), rug)
  colnames(rug) <- c("host", paste0(1:(ncol(rug)-1)))
  rug <- pivot_longer(as.data.frame(rug), !host, names_to = "pair", values_to = "correlation")
  rug$pair <- as.numeric(rug$pair)

  p <- ggplot(rug, aes(x = pair, y = host)) +
    geom_raster(aes(fill = correlation)) +
    scale_fill_gradient2(low = "navy", mid = "white", high = "red", midpoint = 0,
                         guide = guide_colorbar(frame.colour = "black",
                                                ticks.colour = "black"))

  title <- "ASVs"
  xlabel <- "ASV pairs"
  ylabel <- "hosts"
  if(rtype == "phylum") {
    title <- "Phyla"
    xlabel <- "phylum pairs"
  } else if(rtype == "family") {
    title <- "Families/orders/classes"
    xlabel <- "family/order/class pairs"
  }

  p <- p +
    scale_x_continuous(expand = c(0, 0)) +
    labs(fill = "Correlation",
         x = xlabel,
         y = ylabel,
         title = title) +
    scale_y_continuous(expand = c(0, 0)) +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title = element_text(size = 10, face = "plain"),
          plot.title = element_text(size = 12, hjust = 0.5),
          legend.title = element_text(size = 10),
          # legend.position = "bottom",
          plot.margin = margin(t = 20, r = 10, b = 10, l = 0))
  # legend.text = element_text(margin = margin(b = -20)))
  p <- p +
    theme(plot.margin = margin(t = 20, r = 10, b = 10, l = 10))
  if(rtype == "ASV") {
    legend <- get_legend(p)
  }
  p <- p +
    theme(legend.position = "none")
  plots[[length(plots)+1]] <- p
}

# Plot rugs
p1 <- plot_grid(plots[[1]], NULL, plots[[2]], NULL, plots[[3]],
                ncol = 5,
                rel_widths = c(1.5, 0.1, 1.35, 0.1, 1.25),
                labels = c("A", "", "B", "", "C"),
                label_size = 18,
                label_x = -0.02)

# Append legend to the row
p2 <- plot_grid(p1, legend, ncol = 2, rel_widths = c(1, 0.13))

# ------------------------------------------------------------------------------
#   Supplemental Figure 3
#
#   Note: I think this code has been supplanted by what's in `figure_S3.R`.
#   Need to double-check that and delete this.
# ------------------------------------------------------------------------------

# D_combos <- NULL
# i <- 1 # permuted sample to use
# pdir <- paste0("asv_days90_diet25_scale1_scramble-sample-", sprintf("%02d", i))
# fits <- list.files(file.path("output", "model_fits", pdir, "MAP"), full.names = TRUE)
# for(j in 1:length(fits)) {
#   Sigma <- cov2cor(to_clr(readRDS(fits[j]))$Sigma[,,1])
#   if(is.null(D_combos)) {
#     D <- nrow(Sigma)
#     D_combos <- (D^2 - D)/2
#     rug <- matrix(NA, length(fits), D_combos)
#   }
#   rug[j,] <- Sigma[upper.tri(Sigma)]
# }
#
# rug <- cbind(1:nrow(rug), rug)
# colnames(rug) <- c("host", paste0(1:(ncol(rug)-1)))
# rug <- pivot_longer(as.data.frame(rug), !host, names_to = "pair", values_to = "correlation")
# rug$pair <- as.numeric(rug$pair)
#
# p <- ggplot(rug, aes(x = pair, y = host)) +
#   geom_raster(aes(fill = correlation)) +
#   scale_fill_gradient2(low = "navy", mid = "white", high = "red", midpoint = 0,
#                        guide = guide_colorbar(frame.colour = "black",
#                                               ticks.colour = "black")) +
#   scale_x_continuous(expand = c(0, 0)) +
#   labs(fill = "Correlation",
#        x = "ASV-ASV pairs",
#        y = "hosts") +
#   scale_y_continuous(expand = c(0, 0)) +
#   theme(axis.text.y = element_blank(),
#         axis.text = element_text(size = 10),
#         axis.title = element_text(size = 14, face = "plain"),
#         plot.title = element_text(size = 20, hjust = 0.5),
#         legend.title = element_text(size = 14),
#         # legend.position = "bottom",
#         plot.margin = margin(t = 10, r = 10, b = 10, l = 10),
#         axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         axis.ticks.y = element_blank())
#
# ggsave(file.path("output", "figures", "permuted_rug.png"),
#        p,
#        dpi = 100,
#        units = "in",
#        height = 6,
#        width = 8)

# ------------------------------------------------------------------------------
#   Breakout plots of positively and negative correlated taxa
# ------------------------------------------------------------------------------

pred_obj <- get_predictions_host_list(host_list = c("VEX"),
                                      output_dir = "asv_days90_diet25_scale1",
                                      metadata = data$metadata)

l_idx <- which(rugs$phylum$hosts == "VEX")
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

# Plot positive pair
p3 <- ggplot() +
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

cat(paste0("Median correlation across hosts for this pair: ",
           round(median(rugs$ASV$rug[,which(rugs$ASV$tax_idx1 == 31 & rugs$ASV$tax_idx2 == 114)]), 3)))

# Plot negative pair
p4 <- ggplot() +
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

cat(paste0("Median correlation across hosts for this pair: ",
           round(median(rugs$ASV$rug[,which(rugs$ASV$tax_idx1 == 2 & rugs$ASV$tax_idx2 == 3)]), 3)))

# Render breakout plots
p5 <- plot_grid(p3, p4,
                ncol = 2,
                rel_widths = c(1, 1),
                scale = 0.85,
                labels = c("D", "E"),
                label_size = 18,
                label_x = -0.02)

# Combine top and bottom rows (rugs and breakout plots)
p6 <- plot_grid(p2, NULL, p5,
                ncol = 1,
                rel_heights = c(1, 0.05, 0.5))

ggsave("output/figures/F2.png",
       p6,
       dpi = 100,
       units = "in",
       height = 6,
       width = 12)
