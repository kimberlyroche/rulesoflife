source("path_fix.R")

library(driver)
library(tidyverse)
library(rulesoflife)
library(RColorBrewer)
library(cowplot)
library(ggdendro)
library(grid)
library(ggridges)
library(ggraph)
library(igraph)
library(doParallel)
library(foreach)
library(fido)

registerDoParallel(6)

plot_dir <- file.path(check_dir(c("output", "figures")))

# ------------------------------------------------------------------------------
#
#   Figure 1
#
# ------------------------------------------------------------------------------

data <- load_data(tax_level = "family")

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
retain_taxa <- which(apply(relative_abundances[1:(nrow(relative_abundances)-1),], 1, mean) >= 0.025)
collapse_taxa <- setdiff(1:nrow(relative_abundances), retain_taxa)
trimmed_relab <- relative_abundances[retain_taxa,]
trimmed_relab <- rbind(trimmed_relab,
                       colSums(relative_abundances[collapse_taxa,]))

labels <- character(length(retain_taxa))
for(i in 1:length(retain_taxa)) {
  taxon_idx <- retain_taxa[i]
  level <- max(which(!is.na(data$taxonomy[taxon_idx,])))
  labels[i] <- paste0(colnames(data$taxonomy)[level],
                      " ",
                      data$taxonomy[taxon_idx,level])
}
labels <- c(labels, "other")

palette_fn <- file.path("output", "timecourse_palette.rds")
# palette_fn <- file.path("output", "family_palette.rds")
if(file.exists(palette_fn)) {
  palette <- readRDS(palette_fn)
} else {
  # Random palette
  getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
  palette <- sample(getPalette(nrow(trimmed_relab)))
  saveRDS(palette, palette_fn)
}

# Get all hosts
ref_hosts <- unique(data$metadata$sname)

# Get the top 20 best-sampled hosts
# ref_hosts <- sort(data$metadata %>%
#   group_by(sname) %>%
#   tally() %>%
#   arrange(desc(n)) %>%
#   slice(1:20) %>%
#   pull(sname))

# Get the best represented in the primary social groups
# ref_hosts <- c("DUI", "DUX", "LIW", "PEB", "VET")

plots <- list()
legend <- NULL
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

  plot_df <- plot_df %>%
    left_join(data.frame(taxon = levels(plot_df$taxon), name = labels), by = "taxon")
  plot_df$taxon <- plot_df$name

  p <- ggplot(plot_df, aes(x = sample, y = relative_abundance, fill = taxon)) +
    geom_area() +
    scale_fill_manual(values = palette) +
    theme(legend.position = "bottom")
  if(is.null(legend)) {
    legend <- get_legend(p)
  }
  p <- p +
    # geom_area(linetype = 1, size = 0.3, color = "black") +
    theme_nothing() +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme(plot.margin = margin(t = 10, r = , b = 0, l = 1))
  plots[[length(plots) + 1]] <- p
}

p <- plot_grid(plotlist = plots,
               # ncol = length(plots),
               nrow = 4,
               # labels = letters[1:length(plots)],
               labels = ref_hosts,
               scale = 0.95,
               label_x = -0.1,
               label_y = 1,
               label_size = 8)
pl <- plot_grid(p, legend, ncol = 1, rel_heights = c(1, 0.15))
ggsave(file.path(plot_dir, "F1d.svg"),
       pl,
       units = "in",
       dpi = 100,
       # height = 2,
       height = 6,
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
ggsave("output/figures/F1a.svg",
       p4,
       dpi = 100,
       units = "in",
       height = 6,
       width = 18)

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
thresholds_scores <- data.frame(type = factor(c("Phylum", "Family", "ASV")),
                                x0 = c(0.162, 0.140, 0.136))

# "Important"/significant correlations.
thresholds <- data.frame(type = factor(c("Phylum", "Family", "ASV")),
                         lower = c(-0.303, -0.256, -0.263),
                         upper = c(0.149, 0.207, 0.254))

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
  sum(x < thresholds %>% filter(type == "ASV") %>% pull(lower) | x > thresholds %>% filter(type == "ASV") %>% pull(upper))/length(x)
})

p1 <- ggplot(score_df %>% filter(sign != 0)) +
  geom_point(mapping = aes(x = x, y = y, fill = signif, color = sign),
             size = 2,
             shape = 21,
             stroke = 1) +
  scale_fill_distiller(palette = "PuRd", direction = 1,
                       guide = guide_colorbar(frame.colour = "black",
                                              ticks.colour = "black")) +
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
  geom_segment(data = thresholds_scores, aes(x = x0, xend = x0, y = as.numeric(type),
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
               rel_widths = c(1.5, 1),
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

# ------------------------------------------------------------------------------
#
#   Figure 4 (in progress)
#
#   Note: I think the labeling of families is wrong here. Many are being pulled
#         as "unknown" at the family level, which isn't right.
#         Also it would be useful to have a single overall family-level palette!
#
# ------------------------------------------------------------------------------

data <- load_data(tax_level = "ASV")

# Note: I'm repeatedly calculating these. It would be much better to do this once
# globally.
scores <- apply(rugs$c$rug, 2, calc_universality_score)
consensus_signs <- apply(rugs$c$rug, 2, calc_consensus_sign)

# palette_fn <- file.path("output", "family_palette.rds")
# if(file.exists(palette_fn)) {
#   family_palette <- readRDS(palette_fn)
# } else {
  unique_families <- unique(data$taxonomy[,6])
  family_palette <- generate_highcontrast_palette(length(unique_families))
  names(family_palette) <- unique_families
  names(family_palette)[is.na(names(family_palette))] <- "Unknown"
  family_palette[names(family_palette) == "Unknown"] <- "#999999"
  # saveRDS(family_palette, file = file.path("output", "family_palette.rds"))
# }

percent <- 2.5

k <- round(length(scores)*(percent/100))
top_pairs <- order(scores, decreasing = TRUE)[1:k]

pair_idx1 <- rugs$c$tax_idx1[top_pairs]
pair_idx2 <- rugs$c$tax_idx2[top_pairs]

map_df <- data.frame(old_idx = unique(c(pair_idx1, pair_idx2)))
map_df$name <- 1:nrow(map_df)

# ------------------------------------------------------------------------------
#   Build node and edge data.frames
# ------------------------------------------------------------------------------

node_df <- map_df
node_df$Family <- sapply(node_df$old_idx, function(x) {
  data$taxonomy[x,6]
})
node_df <- node_df[,2:3]

edge1 <- data.frame(old_idx = pair_idx1)
edge1 <- left_join(edge1, map_df, by = "old_idx")$name
edge2 <- data.frame(old_idx = pair_idx2)
edge2 <- left_join(edge2, map_df, by = "old_idx")$name

edge_df <- data.frame(from = edge1,
                      to = edge2,
                      Sign = factor(consensus_signs[top_pairs]),
                      score = scores[top_pairs])
levels(edge_df$Sign) <- c("negative", "positive")

# ------------------------------------------------------------------------------
#   Build and plot graph object
# ------------------------------------------------------------------------------

na_idx <- which(is.na(node_df$Family))
node_df$Family[na_idx] <- "Unknown"

graph <- graph_from_data_frame(edge_df, node_df, directed = FALSE)

# Not specifying the layout - defaults to "auto"
# fr and kk layouts are ok here
p <- ggraph(graph, layout = "fr") +
  geom_edge_link(aes(color = Sign), width = 2, alpha = 1) +
  geom_node_point(aes(color = Family), size = 5) +
  geom_node_label(aes(label = name), size = 3, repel = TRUE) +
  # scale_colour_manual(values = family_palette) +
  scale_edge_colour_manual(values = c(negative = "gray", positive = "black")) +
  labs(x = "dimension 1",
       y = "dimension 2") +
  theme_bw()
ggsave(file.path(plot_dir, "F4.svg"),
       p,
       units = "in",
       dpi = 100,
       height = 6,
       width = 10)

# ------------------------------------------------------------------------------
#
#   Figure 5
#
# ------------------------------------------------------------------------------

# Get sampling dates for all hosts
md <- data$metadata
hosts <- sort(unique(md$sname))
host_dates <- list()
for(host in hosts) {
  host_dates[[host]] <- md %>%
    filter(sname == host) %>%
    select(collection_date, sample_id)
}

# Get all pairs of hosts
pairs <- combn(1:length(hosts), 2)

overlap <- "daily"

# ------------------------------------------------------------------------------
#   Calculate (or load) all overlaps of a given frequency
# ------------------------------------------------------------------------------

save_fn <- file.path("output", paste0(overlap, "_host-host_overlap.rds"))
if(!file.exists(save_fn)) {
  # This takes almost 90 min. to run from scratch
  overlap_obj <- data.frame(host1 = c(), host2 = c(),
                            overlap_date1 = c(), overlap_date2 = c(),
                            sample_id1 = c(), sample_id2 = c())
  start <- Sys.time()
  for(i in 1:ncol(pairs)) {
    h1 <- hosts[pairs[1,i]]
    h2 <- hosts[pairs[2,i]]
    cat(paste0("Host pair: ", h1, " x ", h2, " (", i, " / ", ncol(pairs), ")\n"))

    d1 <- host_dates[[h1]]$collection_date
    d2 <- host_dates[[h2]]$collection_date
    start <- min(c(d1, d2))
    if(d1[1] == start) {
      # Start with d1
      ref_h <- h1
      alt_h <- h2
    } else {
      # Start with d2
      ref_h <- h2
      alt_h <- h1
    }

    ref_d <- host_dates[[ref_h]]$collection_date
    alt_d <- host_dates[[alt_h]]$collection_date
    ref_s <- host_dates[[ref_h]]$sample_id
    alt_s <- host_dates[[alt_h]]$sample_id

    for(j in 1:length(ref_d)) {
      delta <- abs(sapply(alt_d, function(x) difftime(x, ref_d[j], units = "days")))
      if(overlap == "daily") {
        hits <- unname(which(delta <= 1))
      }
      if(overlap == "weekly") {
        hits <- unname(which(delta <= 7))
      }
      if(overlap == "monthly") {
        hits <- unname(which(delta <= 30))
      }
      for(hit in hits) {
        overlap_obj <- rbind(overlap_obj,
                             data.frame(host1 = ref_h, host2 = alt_h,
                                        overlap_date1 = ref_d[j], overlap_date2 = alt_d[hit],
                                        sample_id1 = ref_s[j], sample_id2 = alt_s[hit]))
      }
    }
  }
  rt <- Sys.time() - start
  cat("Full run time:", rt, attr(rt, "units"), "\n")

  saveRDS(overlap_obj, paste0(overlap, "_host-host_overlap.rds"))
} else {
  overlap_obj <- readRDS(save_fn)
}

# ------------------------------------------------------------------------------
#   For each host-pair (1540), pull an aligned date. For a given taxon, sample
#   a 1540 length vector of these the estimated Etas pairs at this date for
#   these hosts. Calculate the correlation of these to get an idea of
#   "synchrony."
# ------------------------------------------------------------------------------

# Pull host Etas
# This takes < 30 sec.
Etas <- list()
for(h in 1:length(hosts)) {
  host <- hosts[h]
  fit <- readRDS(file.path("output",
                           "model_fits",
                           "asv_days90_diet25_scale1",
                           "MAP",
                           paste0(hosts[h], ".rds")))
  fit.clr <- to_clr(fit)
  Etas[[host]] <- fit.clr$Eta[,,1]
}

n_tax <- nrow(Etas[[1]])

# Build data.frame of overlaps
# This takes < 30 sec.
sampled_overlap <- NULL
for(i in 1:ncol(pairs)) {
  h1 <- hosts[pairs[1,i]]
  h2 <- hosts[pairs[2,i]]

  temp <- overlap_obj %>%
    filter((host1 == h1 & host2 == h2) | (host1 == h2 & host2 == h1))
  if(nrow(temp) > 0) {
    temp <- temp %>%
      arrange(sample(1:nrow(temp))) %>%
      slice(1)
    if(is.null(sampled_overlap)) {
      sampled_overlap <- temp
    } else {
      sampled_overlap <- rbind(sampled_overlap, temp)
    }
  }
}

if(FALSE) {
  # Scrambled/permuted version
  # Haven't run this in a long time; need to check
  for(i in 1:nrow(sampled_overlap)) {
    h1 <- sampled_overlap[i,]$host1
    h2 <- sampled_overlap[i,]$host2
    s1 <- host_dates[[h1]]
    s1 <- s1[sample(1:nrow(s1), size = 1),]
    s2 <- host_dates[[h2]]
    s2 <- s2[sample(1:nrow(s2), size = 1),]
    sampled_overlap[i,]$overlap_date1 <- s1$collection_date
    sampled_overlap[i,]$sample_id1 <- s1$sample_id
    sampled_overlap[i,]$overlap_date2 <- s2$collection_date
    sampled_overlap[i,]$sample_id2 <- s2$sample_id
  }
}

# Parallelized over 6 cores this takes < 1 min.
starts <- seq(from = 1, to = nrow(sampled_overlap), by = 100)
sampled_list <- foreach(k = 1:length(starts), .combine = rbind) %dopar% {
  start <- starts[k]
  end <- min(c(nrow(sampled_overlap), start + 99))
  subset_data <- sampled_overlap[start:end,]
  row_combos <- nrow(subset_data)*(n_tax-1)*2
  sampled_Etas <- data.frame(Eta = numeric(row_combos),
                             tax_idx = numeric(row_combos),
                             partner = numeric(row_combos))
  row_counter <- 1
  for(i in 1:nrow(subset_data)) {
    # cat(paste0("Processing row ", i, " / ", nrow(subset_data), "\n"))
    s_row <- subset_data[i,]
    h1 <- s_row$host1
    h2 <- s_row$host2
    Eta1 <- Etas[[h1]]
    Eta2 <- Etas[[h2]]

    s1_idx <- which(host_dates[[h1]]$collection_date %in% c(s_row$overlap_date1,
                                                            s_row$overlap_date2))
    s2_idx <- which(host_dates[[h2]]$collection_date %in% c(s_row$overlap_date1,
                                                            s_row$overlap_date2))
    s1_idx <- s1_idx[1]
    s2_idx <- s2_idx[1]

    for(j in 1:(n_tax-1)) {
      sampled_Etas[row_counter,] <- data.frame(Eta = Eta1[j,s1_idx],
                                               tax_idx = j,
                                               partner = 1)
      row_counter <- row_counter + 1
      sampled_Etas[row_counter,] <- data.frame(Eta = Eta2[j,s2_idx],
                                               tax_idx = j,
                                               partner = 2)
      row_counter <- row_counter + 1
    }
  }
  sampled_Etas
}

# ------------------------------------------------------------------------------
#   Plot: Mean correlations for each taxon pair x universality scores
# ------------------------------------------------------------------------------

correlations <- c()
for(i in 1:(n_tax-1)) {
  x <- sampled_list %>%
    filter(partner == 1 & tax_idx == i) %>%
    pull(Eta)
  y <- sampled_list %>%
    filter(partner == 2 & tax_idx == i) %>%
    pull(Eta)
  correlations <- c(correlations, cor(x,y))
}

plot_df <- data.frame(synchrony = c(),
                      universality = c())
for(i in 1:length(rugs$c$tax_idx1)) {
  t1 <- rugs$c$tax_idx1[i]
  t2 <- rugs$c$tax_idx2[i]
  plot_df <- rbind(plot_df,
                   data.frame(synchrony = mean(c(correlations[t1], correlations[t2])),
                              universality = scores[i],
                              sign = consensus_signs[i]))
}
plot_df$sign <- factor(plot_df$sign, levels = c(1, -1))
levels(plot_df$sign) <- c("positive", "negative")

p <- ggplot(plot_df %>% filter(!is.na(sign)), aes(x = synchrony, y = universality, fill = sign)) + #, color = sign)) +
  geom_point(size = 2, shape = 21) +
  scale_fill_manual(values = c("#F25250", "#34CCDE")) +
  theme_bw() +
  xlim(c(min(correlations), 0.46)) +
  # xlim(c(min(correlations), 0.1)) + # for scrambled version
  theme_bw() +
  labs(fill = "Consensus\ncorrelation sign")
show(p)
ggsave(file.path(plot_dir, "F5_permuted.svg"),
       p,
       units = "in",
       dpi = 100,
       height = 5,
       width = 8)

cat(paste0("R^2: ", round(cor(plot_df$synchrony, plot_df$universality)^2, 3), "\n"))





