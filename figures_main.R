source("path_fix.R")

library(driver)
library(tidyverse)
library(rulesoflife)
library(RColorBrewer)
library(cowplot)
library(ggdendro)
library(grid)
library(ggridges)

plot_dir <- file.path(check_dir(c("output", "figures")))

data <- load_data(tax_level = "family")

# ------------------------------------------------------------------------------
#
#   Figure 1
#
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
#   F1b - overview of sampling scheme
# ------------------------------------------------------------------------------

hosts_dates <- data$metadata %>%
  select(sname, collection_date)
hosts_dates$sname <- factor(hosts_dates$sname,
                            levels = sort(unique(hosts_dates$sname), decreasing = TRUE))
baseline_time <- min(hosts_dates$collection_date)
hosts_dates$time <- as.numeric(sapply(hosts_dates$collection_date, function(x) difftime(x, baseline_time, units = "days")))

xticks <- seq(from = 0, to = 5000, length = 20)
xlabs <- character(length(xticks))
for(i in 1:length(xticks)) {
  xlabs[i] <- as.character(as.Date(baseline_time) + xticks[i])
}

p <- ggplot(hosts_dates, aes(x = time, y = sname)) +
  geom_point() +
  theme_bw() +
  labs(x = "days from first sample",
       y = "host short name") +
  scale_x_continuous(breaks = xticks, labels = xlabs) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
p
ggsave(file.path(plot_dir, "F1b.svg"),
       p,
       units = "in",
       dpi = 100,
       height = 7,
       width = 8)

# ------------------------------------------------------------------------------
#   F1d - summarized time courses for 5 reference hosts
# ------------------------------------------------------------------------------

relative_abundances <- data$counts
for(j in 1:ncol(relative_abundances)) {
  relative_abundances[,j] <- relative_abundances[,j]/sum(relative_abundances[,j])
}

# Collapse rare stuff for visualization
retain_taxa <- which(apply(relative_abundances, 1, mean) >= 0.01)
collapse_taxa <- setdiff(1:nrow(relative_abundances), retain_taxa)
trimmed_relab <- relative_abundances[retain_taxa,]
trimmed_relab <- rbind(trimmed_relab,
                       colSums(relative_abundances[collapse_taxa,]))

getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
palette <- sample(getPalette(nrow(trimmed_relab)))

# Get all hosts
# ref_hosts <- unique(data$metadata$sname)

# Get the top 20 best-sampled hosts
ref_hosts <- sort(data$metadata %>%
  group_by(sname) %>%
  tally() %>%
  arrange(desc(n)) %>%
  slice(1:20) %>%
  pull(sname))

# Get the best represented in the primary social groups
# ref_hosts <- c("DUI", "DUX", "LIW", "PEB", "VET")

plots <- list()

for(host in ref_hosts) {
  host_relab <- trimmed_relab[,data$metadata$sname == host]

  # Downsample
  ds_idx <- round(seq(1, ncol(host_relab), length.out = min(ncol(host_relab), 20)))
  ds_idx[1] <- 1
  ds_idx[length(ds_idx)] <- ncol(host_relab)
  host_ds <- host_relab[,ds_idx]

  plot_df <- cbind(1:nrow(host_ds), as.data.frame(host_ds))
  colnames(plot_df) <- c("taxon", 1:(ncol(plot_df)-1))
  plot_df <- pivot_longer(plot_df, !taxon, names_to = "sample", values_to = "relative_abundance")
  plot_df$taxon <- factor(plot_df$taxon)
  plot_df$sample <- as.numeric(plot_df$sample)

  p <- ggplot(plot_df, aes(x = sample, y = relative_abundance, fill = taxon)) +
    geom_area() +
    # geom_area(linetype = 1, size = 0.3, color = "black") +
    scale_fill_manual(values = palette) +
    theme_nothing() +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme(plot.margin = margin(t = 10, r = , b = 0, l = 1))
  plots[[length(plots) + 1]] <- p
}

p <- plot_grid(plotlist = plots,
               ncol = length(plots),
               # labels = letters[1:length(plots)],
               labels = ref_hosts,
               scale = 0.95,
               label_x = -0.2,
               label_y = 1,
               label_size = 8)
p
ggsave(file.path(plot_dir, "F1d.svg"),
       p,
       units = "in",
       dpi = 100,
       height = 2,
       width = 10)

# ------------------------------------------------------------------------------
#
#   Figure 2
#
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
#   F2a - three rugs
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
    scale_fill_gradient2(low = "navy", mid = "white", high = "red",
                         midpoint = 0)

  if(is.null(row_labels)) {
    row_labels <- 1:length(rug$host)
  }
  p <- p +
    scale_x_continuous(expand = c(0, 0)) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank()) +
    labs(fill = "Correlation",
         x = "taxon pairs",
         y = "")
  p <- p +
    scale_y_continuous(breaks = 1:length(row_labels),
                       labels = row_labels,
                       expand = c(0, 0)) +
    theme(axis.text.y = element_text(size = rel(0.9)))
  p <- p +
    theme(axis.text = element_text(size = 10),
          axis.title = element_text(size = 14, face = "plain"),
          legend.title = element_text(size = 14),
          plot.margin = margin(t = 20, r = 10, b = 10, l = 0))
  if(rtype == "a") {
    p <- p +
      theme(plot.margin = margin(t = 20, r = 10, b = 10, l = 0))
  } else {
    p <- p +
      theme(plot.margin = margin(t = 20, r = 10, b = 10, l = 10))
  }
  if(is.null(legend)) {
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
p1 <- plot_grid(nullGrob(), plots[[1]], nullGrob(), ncol = 1, rel_heights = c(0.0, 1, 0.02))
p2 <- plot_grid(p1, plots[[2]],
                ncol = 2,
                rel_widths = c(0.1, 0.9))
p3 <- plot_grid(p2, plots[[3]], plots[[4]],
                ncol = 3,
                scale = 0.92,
                rel_widths = c(1, 1.2, 1.2),
                labels = c("a", "b", "c"),
                label_size = 38,
                label_x = 0,
                label_y = 1.01)
ggsave("output/figures/F1a.svg",
       p3,
       dpi = 100,
       units = "in",
       height = 8,
       width = 24)

# ------------------------------------------------------------------------------
#   F2b - breakout plots of positively and negative correlated taxa
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
ggsave("output/figures/F1b.svg",
       p,
       dpi = 100,
       units = "in",
       height = 2,
       width = 6)

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
ggsave("output/figures/F1c.svg",
       p,
       dpi = 100,
       units = "in",
       height = 2,
       width = 6)

# ------------------------------------------------------------------------------
#
#   Figure 3
#
# ------------------------------------------------------------------------------

# These are hard-coded. See `figures_supplemental.R` for the calculation of
# "important" universality scores from permuted data.
thresholds <- data.frame(type = factor(c("Phylum", "Family", "ASV")),
                         x0 = c(0.162, 0.140, 0.136))

# ------------------------------------------------------------------------------
#   F3a - "hockeystick" plot
# ------------------------------------------------------------------------------

asv_scores_pieces <- apply(rugs$c$rug, 2, function(x) calc_universality_score(x, return_pieces = TRUE))
score_df <- data.frame(x = asv_scores_pieces[1,], y = asv_scores_pieces[2,])

# Original version, with a background indicating strength of potential scores
if(FALSE) {
  theoretical_scores <- NULL
  for(x in seq(from = 0.5, to = 1, length.out = 30)) {
    for(y in seq(from = 0, to = 1, length.out = 30)) {
      theoretical_scores <- rbind(theoretical_scores,
                                  data.frame(x = x, y = y, z = x*y))
    }
  }

  p1 <- ggplot() +
    geom_tile(data = theoretical_scores,
              mapping = aes(x = x, y = y, fill = z)) +
    scale_fill_distiller(palette = "Blues",
                         trans = "reverse",
                         breaks = seq(from = 0.1, to = 0.9, length.out = 5),
                         guide = guide_colorbar(frame.colour = "black")) +
    geom_point(data = score_df,
               mapping = aes(x = x, y = y),
               size = 2,
               shape = 21,
               stroke = 0.5,
               color = "#222222",
               fill = "#aaaaaa") +
    theme_bw() +
    theme(legend.title = element_text(margin = margin(b = 5))) +
    labs(x = "proportion shared sign",
         y = "median association strength",
         fill = "Universality score") +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0))
}

# Additional labelings
# 1) Consensus sign
score_df$sign <- apply(rugs$c$rug, 2, calc_consensus_sign)
score_df$sign <- factor(score_df$sign, levels = c(1, -1))
levels(score_df$sign) <- c("positive", "negative")
# 2) Percent significant observations for this taxon pair
score_df$signif <- apply(rugs$c$rug, 2, function(x) {
  sum(x > thresholds %>% filter(type == "ASV") %>% pull(x0))/length(x)
})

p1 <- ggplot(score_df %>% filter(sign != 0)) +
  geom_point(mapping = aes(x = x, y = y, fill = signif, color = sign),
             size = 2,
             shape = 21,
             stroke = 1) +
  scale_fill_distiller(palette = "PuRd", direction = 1) +
  scale_color_manual(values = c("#000000", "#888888")) +
  theme_bw() +
  ylim(c(0,1)) +
  theme(legend.title = element_text(margin = margin(b = 5))) +
  labs(x = "proportion shared sign",
       y = "median association strength",
       fill = "Proportion\nsignificant\nobservations",
       color = "Consensus\ncorrelation sign")

# ------------------------------------------------------------------------------
#   F3b - ggridges plot of universality scores by taxonomic level
# ------------------------------------------------------------------------------

plot_df <- data.frame(score = apply(rugs$a$rug, 2, calc_universality_score),
                      type = "Phylum")
plot_df <- rbind(plot_df,
                 data.frame(score = apply(rugs$b$rug, 2, calc_universality_score),
                            type = "Family"))
plot_df <- rbind(plot_df,
                 data.frame(score = apply(rugs$c$rug, 2, calc_universality_score),
                            type = "ASV"))

p2 <- ggplot(plot_df, aes(x = score, y = type, fill = type)) +
  geom_density_ridges(stat = "binline", bins = 20, scale = 0.9, draw_baseline = TRUE) +
  geom_segment(data = thresholds, aes(x = x0, xend = x0, y = as.numeric(type),
                                      yend = as.numeric(type) + .9),
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
  xlim(c(-0.05, 0.65)) +
  scale_fill_brewer(palette = "Blues")

# Ultimately want to combine these, a laoo
p <- plot_grid(p1, p2, ncol = 2,
               rel_widths = c(1.75, 1),
               labels = c("a", "b"),
               label_size = 18,
               # label_x = -0.02,
               scale = 0.95)

ggsave("output/figures/F3.svg",
       p,
       dpi = 100,
       units = "in",
       height = 4,
       width = 9)

