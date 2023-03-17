source("path_fix.R")

library(driver)
library(fido)
library(tidyverse)
library(rulesoflife)
library(RColorBrewer)
library(cowplot)

source("thresholds.R")

data <- load_data(tax_level = "ASV")
md <- data$metadata
hosts <- unique(md$sname)
host_shortlist <- c("DUI", "DUX", "LIW", "PEB", "VET")

host_labels <- read.delim(file.path("output", "host_labels.tsv"),
                          header = TRUE)
host_order <- c(2,1,4,5,3)

# These will be in the order DUX, DUI, PEB, VET, LIW
use_labels <- c()
for(this_host in host_shortlist[host_order]) {
  use_labels <- c(use_labels,
                  host_labels %>% filter(sname == this_host) %>% pull(host_label))
}

# ------------------------------------------------------------------------------
#   Observed ASV-level "rug"
# ------------------------------------------------------------------------------

rug_obj <- summarize_Sigmas(output_dir = "asv_days90_diet25_scale1")

filtered_pairs <- filter_joint_zeros(data$counts, threshold_and = 0.05, threshold_or = 0.5)
filtered_obj <- rug_obj
filtered_obj$tax_idx1 <- filtered_obj$tax_idx1[filtered_pairs$threshold]
filtered_obj$tax_idx2 <- filtered_obj$tax_idx2[filtered_pairs$threshold]
filtered_obj$rug <- filtered_obj$rug[,filtered_pairs$threshold]

rug_obj <- filtered_obj

# Most negatively correlated pair
# idx <- which.min(apply(rug_obj$rug[rug_obj$hosts %in% host_shortlist,], 2, median))
idx <- which.min(apply(rug_obj$rug, 2, median))
rug_obj$tax_idx1[idx]
rug_obj$tax_idx2[idx]
median(rug_obj$rug[,idx])
sd(rug_obj$rug[,idx])

# Most positively correlated pair
# idx <- which.max(apply(rug_obj$rug[rug_obj$hosts %in% host_shortlist,], 2, median))
idx <- which.max(apply(rug_obj$rug, 2, median))
rug_obj$tax_idx1[idx]
rug_obj$tax_idx2[idx]
median(rug_obj$rug[,idx])
sd(rug_obj$rug[,idx])

# ------------------------------------------------------------------------------
#   Enrichment statistics for sign
# ------------------------------------------------------------------------------

asv_column_order <- NULL
legend <- NULL

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
rug <- rug_obj$rug
canonical_col_order <- order(colMeans(rug))
canonical_row_order <- row_order
rug <- rug[canonical_row_order,canonical_col_order]
asv_column_order <- canonical_col_order

# cmeans <- colMeans(rug)
# temp <- data.frame(index = 1:ncol(rug),
#                    cmean = colMeans(rug))
# temp$above_upper <- temp$cmean > thresholds$upper[thresholds$type == "ASV"]
# temp$below_lower <- temp$cmean < thresholds$lower[thresholds$type == "ASV"]
# x_lower <- temp %>%
#   filter(below_lower) %>%
#   arrange(desc(cmean)) %>%
#   slice_head() %>%
#   pull(index)
# x_upper <- temp %>%
#   filter(above_upper) %>%
#   arrange(cmean) %>%
#   slice_head() %>%
#   pull(index)

rug <- cbind(1:nrow(rug), rug)
colnames(rug) <- c("host", paste0(1:(ncol(rug)-1)))
rug <- cbind(host_name = host_labels$host_label, rug)
rug <- pivot_longer(as.data.frame(rug), !c(host_name, host), names_to = "pair", values_to = "correlation")
rug$pair <- as.numeric(rug$pair)
rug$correlation <- as.numeric(rug$correlation)

temp <- rug %>%
  dplyr::select(host, host_name) %>%
  distinct() %>%
  arrange(host_name)
name_map <- unlist(temp$host_name)
names(name_map) <- temp$host

# Mute host names not in the shortlist
name_map <- name_map[name_map %in% use_labels]
plot_breaks <- names(name_map)

p2 <- ggplot(rug, aes(x = pair, y = host)) +
  geom_raster(aes(fill = correlation)) +
  # geom_segment(data = data.frame(x_lower = x_lower),
  #              mapping = aes(x = x_lower, xend = x_lower, y = 1, yend = 56), size = 0.5, linetype = "dashed") +
  # geom_segment(data = data.frame(x_upper = x_upper),
  #              mapping = aes(x = x_upper, xend = x_upper, y = 1, yend = 56), size = 0.5, linetype = "dashed") +
  scale_fill_gradientn(limits = c(-1,1), colors = c("navy", "white", "red"),
                       guide = guide_colorbar(frame.colour = "black",
                                              ticks.colour = "black")) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(labels = name_map, breaks = plot_breaks) +
  labs(fill = "Correlation",
       x = "ASV pairs",
       y = "hosts",
       title = "ASVs") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 12, face = "plain"),
        plot.title = element_blank(),
        plot.margin = margin(t = 20, r = 10, b = 10, l = 10))
legend <- get_legend(p2)
p2 <- p2 +
  theme(legend.position = "none")

rug_na <- rug %>%
  mutate(correlation = case_when(
    correlation >= thresholds$upper[thresholds$type == "ASV"] ~ correlation,
    correlation <= thresholds$lower[thresholds$type == "ASV"] ~ correlation,
    T ~ NA_real_
  ))

p2_na <- ggplot(rug_na, aes(x = pair, y = host)) +
  geom_raster(aes(fill = correlation)) +
  scale_fill_gradientn(limits = c(-1,1), colors = c("navy", "white", "red"),
                       na.value = "black",
                       guide = guide_colorbar(frame.colour = "black",
                                              ticks.colour = "black")) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(labels = name_map, breaks = plot_breaks) +
  labs(fill = "Correlation",
       x = "ASV pairs",
       y = "hosts",
       title = "ASVs") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 12, face = "plain"),
        plot.title = element_blank(),
        plot.margin = margin(t = 20, r = 10, b = 10, l = 10))

ggsave(file.path("output", "figures", "Figure_2_Supplement_4.png"),
       p2_na,
       dpi = 200,
       units = "in",
       height = 5,
       width = 8,
       bg = "white")

# ------------------------------------------------------------------------------
#   Scrambled/permuted "rug"
# ------------------------------------------------------------------------------

D_combos <- NULL
i <- 1 # permuted sample to use
pdir <- paste0("asv_days90_diet25_scale1_scramble-sample-", sprintf("%02d", i))
fits <- list.files(file.path("output", "model_fits", pdir, "MAP"), full.names = TRUE)
for(j in 1:length(fits)) {
  Sigma <- cov2cor(to_clr(readRDS(fits[j]))$Sigma[,,1])
  Sigma <- Sigma[1:(nrow(Sigma)-1),1:(ncol(Sigma)-1)]
  if(is.null(D_combos)) {
    D <- nrow(Sigma)
    D_combos <- (D^2 - D)/2
    rug <- matrix(NA, length(fits), D_combos)
  }
  rug[j,] <- Sigma[upper.tri(Sigma)]
}

# Cut this down to low 0-0 frequency pairs
rug <- rug[,filtered_pairs$threshold]

rug <- cbind(1:nrow(rug), rug[,asv_column_order])
colnames(rug) <- c("host", paste0(1:(ncol(rug)-1)))
rug <- cbind(host_name = host_labels$host_label, rug)
rug <- pivot_longer(as.data.frame(rug), !c(host_name, host), names_to = "pair", values_to = "correlation")
rug$pair <- as.numeric(rug$pair)
rug$correlation <- as.numeric(rug$correlation)

p1 <- ggplot(rug, aes(x = pair, y = host)) +
  geom_raster(aes(fill = correlation)) +
  scale_fill_gradientn(limits = c(-1,1), colors = c("navy", "white", "red"),
                       guide = guide_colorbar(frame.colour = "black",
                                              ticks.colour = "black")) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(labels = name_map, breaks = plot_breaks) +
  labs(fill = "Correlation",
       x = "ASV pairs",
       y = "hosts",
       title = "ASVs") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 12, face = "plain"),
        plot.title = element_blank(),
        plot.margin = margin(t = 20, r = 10, b = 10, l = 10),
        legend.position = "none")

# Plot rugs
p <- plot_grid(p2, NULL, p1, ncol = 3,
               rel_widths = c(1, 0.05, 1),
               labels = c("A", "", "B"),
               label_size = 18,
               label_x = -0.015)

# Append legend to the row
prow1 <- plot_grid(p, legend, ncol = 2, rel_widths = c(1, 0.13))

# ------------------------------------------------------------------------------
#   Breakout plots of positively and negative correlated taxa
# ------------------------------------------------------------------------------

alpha <- 0.6

# Find common baseline
pred_objs <- NULL
min_date <- NULL
max_date <- NULL
for(host in host_shortlist[host_order]) {
  pred_obj <- get_predictions_host_list(host_list = host,
                                        output_dir = "asv_days90_diet25_scale1",
                                        metadata = data$metadata)
  if(is.null(min_date)) {
    min_date <- min(pred_obj$dates[[host]])
    max_date <- max(pred_obj$dates[[host]])
  } else {
    min_date <- min(min_date, min(pred_obj$dates[[host]]))
    max_date <- max(max_date, max(pred_obj$dates[[host]]))
  }
  pred_objs[[host]] <- pred_obj
}

render_trajectories <- function(tax_indices, host_shortlist, host_labels, host_y_offset = 5) {
  plot_df <- NULL
  for(tax_idx in tax_indices) {
    y_offset_counter <- 0
    for(h in 1:length(host_shortlist)) {
      host <- host_shortlist[h]
      label <- host_labels[h]
      pred_obj <- pred_objs[[host]]
      pred_df <- suppressMessages(gather_array(pred_obj$predictions[[host]]$Eta,
                                               val,
                                               coord,
                                               sample,
                                               iteration) %>%
                                    filter(coord == tax_idx) %>%
                                    group_by(coord, sample) %>%
                                    summarize(p25 = quantile(val, probs = c(0.25)),
                                              mean = mean(val),
                                              p75 = quantile(val, probs = c(0.75))))
      pred_df$host <- label

      # Calculate day from (shared) baseline
      addend <- as.numeric(difftime(min(pred_obj$dates[[host]]), min_date, units = "day"))

      day_span <- pred_obj$predictions[[host]]$span
      pred_df <- pred_df %>%
        left_join(data.frame(sample = 1:length(day_span), day = day_span), by = "sample")

      pred_df$day <- pred_df$day + addend

      if(y_offset_counter > 0) {
        pred_df$p25 <- pred_df$p25 + host_y_offset*y_offset_counter
        pred_df$mean <- pred_df$mean + host_y_offset*y_offset_counter
        pred_df$p75 <- pred_df$p75 + host_y_offset*y_offset_counter
      }
      y_offset_counter <- y_offset_counter + 1

      plot_df <- rbind(plot_df, pred_df)
    }
  }

  p3 <- ggplot()
  swap_color <- FALSE
  for(tax_idx in tax_indices) {
    for(this_host in host_labels) {
      p3 <- p3 +
        geom_ribbon(data = plot_df %>% filter(host == this_host & coord == tax_idx),
                    mapping = aes(x = day, ymin = p25, ymax = p75),
                    fill = ifelse(swap_color, "#a6cee3", "#fdbf6f"),
                    alpha = alpha) +
        geom_line(data = plot_df %>% filter(host == this_host & coord == tax_idx),
                  mapping = aes(x = day, y = mean),
                  color = ifelse(swap_color, "#1f78b4", "#ff7f00"),
                  size = 1,
                  alpha = alpha)
    }
    swap_color <- TRUE
  }

  breaks <- plot_df %>%
    group_by(host) %>%
    summarize(y_mean = mean(mean))
  name_map <- unlist(breaks$host)
  names(name_map) <- round(unlist(breaks$y_mean))

  p3 <- p3 +
    theme_bw() +
    labs(x = paste0("days from first sample (", min_date, ")"),
         y = "") +
    scale_x_continuous(expand = c(0.01, 0.01)) +
    scale_y_continuous(breaks = unlist(breaks$y_mean), labels = name_map) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())

  p3
}

# Pick a strongly correlated pair
column_median <- apply(rug_obj$rug, 2, median)
strongest_correlators <- order(column_median, decreasing = T)
positive1 <- rug_obj$tax_idx1[strongest_correlators[10]]
positive2 <- rug_obj$tax_idx2[strongest_correlators[10]]
p3 <- render_trajectories(c(positive1, positive2), host_shortlist[c(2,3,4,5,1)], use_labels[c(1,5,3,4,2)], host_y_offset = 12)

# Pick a strongly negatively correlated pair
strongest_inverse <- order(column_median)
negative1 <- rug_obj$tax_idx1[strongest_inverse[10]]
negative2 <- rug_obj$tax_idx2[strongest_inverse[10]]
p4 <- render_trajectories(c(20, 25), host_shortlist[c(2,3,4,5,1)], use_labels[c(1,5,3,4,2)], host_y_offset = 12)

representation <- represented_taxa(filtered_pairs)
cat(paste0("Strongly positively correlated pair are ASV",
           renumber_taxon(representation, positive1),
           " and ASV",
           renumber_taxon(representation, positive2),
           " (median r = ", round(column_median[strongest_correlators[10]], 3), ")\n"))
cat(paste0("Strongly negatively correlated pair are ASV",
           renumber_taxon(representation, negative1),
           " and ASV",
           renumber_taxon(representation, negative2),
           " (median r = ", round(column_median[strongest_inverse[10]], 3), ")\n"))

prow2 <- plot_grid(p4, NULL, p3, NULL, ncol = 4,
                   labels = c("C", "", "D", ""),
                   rel_widths = c(1, 0.05, 1, 0.27),
                   label_size = 18,
                   label_x = -0.015,
                   scale = 0.98)

p_out <- plot_grid(prow1, prow2, ncol = 1,
                   rel_heights = c(1, 0.75))

ggsave(file.path("output", "figures", "Figure_2.png"),
       p_out,
       # dpi = 50,
       dpi = 200,
       units = "in",
       height = 9,
       width = 12,
       bg = "white")
