source("path_fix.R")

library(driver)
library(fido)
library(tidyverse)
library(rulesoflife)
library(RColorBrewer)
library(cowplot)

rugs <- list(ASV = summarize_Sigmas(output_dir = "asv_days90_diet25_scale1"),
             family = summarize_Sigmas(output_dir = "fam_days90_diet25_scale1"),
             phylum = summarize_Sigmas(output_dir = "phy_days90_diet25_scale1"))

# ------------------------------------------------------------------------------
#   Enrichment statistics for sign
# ------------------------------------------------------------------------------

# Phylum
filtered_pairs_phy <- filter_joint_zeros(load_data(tax_level = "phylum")$counts, threshold_and = 0.05, threshold_or = 0.5)$threshold
signif(median(c(rugs$phylum$rug[,filtered_pairs_phy])), 3)
phylum_sign <- sign(c(rugs$phylum$rug[,filtered_pairs_phy]))
binom.test(table(phylum_sign)[2], length(phylum_sign), 0.5)
round(table(phylum_sign)[1]/length(phylum_sign), 2)

# Family/order/class
filtered_pairs_fam <- filter_joint_zeros(load_data(tax_level = "family")$counts, threshold_and = 0.05, threshold_or = 0.5)$threshold
signif(median(c(rugs$family$rug[,filtered_pairs_fam])), 3)
family_sign <- sign(c(rugs$family$rug[,filtered_pairs_fam]))
binom.test(table(family_sign)[2], length(family_sign), 0.5)
round(table(family_sign)[1]/length(family_sign), 2)

# ASV
filtered_pairs_asv <- filter_joint_zeros(load_data(tax_level = "ASV")$counts, threshold_and = 0.05, threshold_or = 0.5)$threshold
signif(median(c(rugs$ASV$rug[,filtered_pairs_asv])), 3)
ASV_sign <- sign(c(rugs$ASV$rug[,filtered_pairs_asv]))
binom.test(table(ASV_sign)[2], length(ASV_sign), 0.5)
round(table(ASV_sign)[1]/length(ASV_sign), 2)

filtered_rugs <- list(ASV = rugs[["ASV"]]$rug[,filtered_pairs_asv],
                      family = rugs[["family"]]$rug[,filtered_pairs_fam],
                      phylum = rugs[["phylum"]]$rug[,filtered_pairs_phy])

asv_column_order <- NULL
plots <- list()
legend <- NULL
for(rtype in names(rugs)) {
  # Order rows by similarity of baseline composition
  data <- load_data(tax_level = rtype)
  md <- data$metadata
  hosts <- unique(md$sname)
  H <- length(hosts)
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
  # rug <- rugs[[rtype]]$rug
  rug <- filtered_rugs[[rtype]]
  canonical_col_order <- order(colMeans(rug))
  canonical_row_order <- row_order
  rug <- rug[canonical_row_order,canonical_col_order]

  if(rtype == "ASV") {
    asv_column_order <- canonical_col_order
  }

  rug <- cbind(1:nrow(rug), rug)
  colnames(rug) <- c("host", paste0(1:(ncol(rug)-1)))
  rug <- pivot_longer(as.data.frame(rug), !host, names_to = "pair", values_to = "correlation")
  rug$pair <- as.numeric(rug$pair)

  p <- ggplot(rug, aes(x = pair, y = host)) +
    geom_raster(aes(fill = correlation)) +
    scale_fill_gradientn(limits = c(-1,1), colors = c("navy", "white", "red"),
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
         y = ylabel) +
    # title = title) +
    scale_y_continuous(expand = c(0, 0)) +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title = element_text(size = 11, face = "plain"),
          plot.title = element_text(size = 12, hjust = 0.5),
          legend.title = element_text(size = 10, vjust = 0.85),
          legend.position = "bottom",
          plot.margin = margin(t = 20, r = 10, b = 10, l = 0))
  p <- p +
    theme(plot.margin = margin(t = 20, r = 10, b = 10, l = 10))

  if(rtype == "ASV") {
    legend <- get_legend(p)
  }
  p <- p +
    theme(legend.position = "none")
  plots[[length(plots)+1]] <- p
}

# Plot rugs; omit ASV-level
p1 <- plot_grid(plots[[2]], NULL, plots[[3]],
                ncol = 3,
                rel_widths = c(1.35, 0.1, 1.25),
                labels = c("A", "", "B"),
                label_size = 18,
                label_x = -0.02)

# Append legend to the row
p2 <- plot_grid(p1, legend, ncol = 1, rel_heights = c(1, 0.13))

ggsave(file.path("output", "figures", "Figure_2_Supplement_1.png"),
       p2,
       dpi = 200,
       units = "in",
       height = 5,
       width = 10,
       bg = "white")
