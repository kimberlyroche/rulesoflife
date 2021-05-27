source("path_fix.R")

library(rulesoflife)
library(tidyverse)
library(driver)

source("ggplot_fix.R")

plot_dir <- check_dir(c("output", "figures"))
output_dir <- "asv_days90_diet25_scale1"

data <- load_data(tax_level = "ASV")
metadata <- data$metadata

host_filename <- file.path("output", "overlapped_hosts.rds")
if(!file.exists(host_filename)) {
  stop(paste0("File not found: ", host_filename, "\n"))
}
selected_hosts <- readRDS(host_filename)

pred_filename <- file.path("output", "host_mean_predictions_adjusted.rds")
if(!file.exists(pred_filename)) {
  stop(paste0("File not found: ", pred_filename, "\n"))
}
predictions <- readRDS(pred_filename)

rug_filename <- file.path("output", "rug_asv.rds")
if(file.exists(rug_filename)) {
  rug_obj <- readRDS(rug_filename)
} else {
  rug_obj <- summarize_Sigmas(output_dir)
  saveRDS(rug_obj, file = rug_filename)
}

# Subset to selected hosts
subset_idx <- rug_obj$hosts %in% selected_hosts
rug_subset <- rug_obj$rug[subset_idx,]
rug_host_subset <- rug_obj$hosts[subset_idx]

universalities <- apply(rug_subset, 2, calc_universality_score)
consensus_sign <- apply(rug_subset, 2, calc_consensus_sign)

# ------------------------------------------------------------------------------
#   Pull interesting pairs (strong and weak universality)
# ------------------------------------------------------------------------------

positive_idx <- which(consensus_sign > 0)
negative_idx <- which(consensus_sign < 0)
positive_ranks <- order(universalities[positive_idx], decreasing = TRUE)
negative_ranks <- order(universalities[negative_idx], decreasing = TRUE)

# k is 1 here
top_k_positive <- positive_ranks[1]
# Indexed as: universalities[positive_idx[positive_ranks[1:10]]]
top_k_negative <- negative_ranks[1]
# Indexed as: universalities[negative_idx[negative_ranks[1:10]]]
bottom_k <- order(universalities)[1]
# Indexed as: universalities[bottom_k]

# ------------------------------------------------------------------------------
#   PLot distributions for selected pairs
# ------------------------------------------------------------------------------

named_pairs <- list(positive1 = positive_idx[positive_ranks[1]],
                    negative1 = negative_idx[negative_ranks[1]],
                    neutral1 = bottom_k[1])

# Each iteration of this loop takes about 40 sec. if we look at all pairs of
# hosts. If we subset to 100 hosts, it's about 9 sec.
subset_hosts <- TRUE
for(pair_name in names(named_pairs)) {
  cat("Evaluating pair", pair_name, "...\n")
  pair <- named_pairs[[pair_name]]

  coord1 = rug_obj$tax_idx1[pair]
  coord2 = rug_obj$tax_idx2[pair]

  plot_df <- build_between_within_distributions(pair,
                                                coord1 = coord1,
                                                coord2 = coord2,
                                                predictions = predictions %>%
                                                  filter(coord %in% c(coord1, coord2)),
                                                selected_hosts = selected_hosts)

  label1 <- get_tax_label(data$taxonomy, coord1, "clr")
  label2 <- get_tax_label(data$taxonomy, coord2, "clr")
  alpha <- 0.5
  p <- ggplot() +
    geom_histogram(data = plot_df[plot_df$type == "between hosts (1)",],
                   mapping = aes(x = correlation, fill = type),
                   color = "white",
                   alpha = alpha) +
    geom_histogram(data = plot_df[plot_df$type == "between hosts (2)",],
                   mapping = aes(x = correlation, fill = type),
                   color = "white",
                   alpha = alpha) +
    geom_histogram(data = plot_df[plot_df$type == "within hosts",],
                   mapping = aes(x = correlation, fill = type),
                   color = "white",
                   alpha = alpha) +
    scale_fill_manual(values = c("#a84a32", "#a84a32", "#326da8")) +
    labs(fill = "Correlation type",
         x = "CLR correlation",
         title = paste0("Between- and within-host correlation of series\n",
                        "(1) ", label1,
                        "\n",
                        "(2) ", label2,
                        "\n"))
  ggsave(file.path(plot_dir, paste0("within-between_", pair_name, ".png")),
         plot = p,
         units = "in",
         dpi = 100,
         height = 4.5,
         width = 6)
  show(p)
}
