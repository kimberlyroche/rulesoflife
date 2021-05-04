library(rulesoflife)
library(tidyverse)

plot_dir <- check_dir(c("output", "figures"))
output_dir <- "asv_days90_diet25_scale1"

# ------------------------------------------------------------------------------
#   Visualize overlap in hosts; select a subset with good overlap in time
# ------------------------------------------------------------------------------

data <- load_data(tax_level = "ASV")
metadata <- data$metadata

# Get max and min day per host, relative to a common baseline
common_baseline <- min(metadata$collection_date)
all_hosts <- unique(metadata$sname)
host_days <- data.frame(min_day = c(),
                        max_day = c(),
                        median_day = c(),
                        host = c())
for(host in all_hosts) {
  days <- sapply(metadata[metadata$sname == host,]$collection_date, function(x) {
    round(difftime(x, common_baseline, units = "days"))
  })
  host_days <- rbind(host_days,
                     data.frame(min_day = min(days),
                                max_day = max(days),
                                median_day = min(days) + (max(days) - min(days))/2,
                                host = host))
}

host_days <- host_days %>%
  arrange(median_day)
host_days$index <- 1:nrow(host_days)

host_days$type <- factor(sapply(1:nrow(host_days), function(i) {
  host_days$min_day[i] <= 750 && host_days$max_day[i] >= 4000
}))
levels(host_days$type) <- c("excluded", "included")

ggplot(host_days, aes(x = min_day,
                      xend = max_day,
                      y = index,
                      yend = index,
                      color = type)) +
  geom_segment(size = 1) +
  scale_color_manual(values = c("#aaaaaa", "#000000")) +
  scale_y_continuous(breaks = host_days$index,
                   labels = host_days$host) +
  theme(axis.text.y = element_text(size = 6)) +
  labs(x = "sampled days",
       y = "host",
       color = "Host status")
ggsave(file.path(plot_dir, "host_day_overlap.png"),
       units = "in",
       dpi = 100,
       height = 5,
       width = 8)

selected_hosts <- host_days[host_days$type == "included",]$host
cat("No. selected hosts:", length(selected_hosts), "\n")

# ------------------------------------------------------------------------------
#   Render any predictions that haven't already been run for selected hosts
# ------------------------------------------------------------------------------

pred_dir <- check_dir(c("output", "model_fits", output_dir, "predictions"))
for(host in selected_hosts) {
  if(!file.exists(file.path(pred_dir, paste0(host, ".rds")))) {
    cat("Predicting on", host, "...\n")
    temp <- predict_trajectory(host, output_dir, 5)
  }
}

# Runtime ~30 sec.
pred_obj <- get_predictions_host_list(selected_hosts, output_dir, metadata)

# ------------------------------------------------------------------------------
#   Functions
# ------------------------------------------------------------------------------

pull_series <- function(filtered_obj, sname, idx) {
  filtered_obj %>%
    filter(host == sname) %>%
    arrange(day) %>%
    filter(coord == idx) %>%
    pull(centered_series)
}

# ------------------------------------------------------------------------------
#   Pull interesting pairs (strong and weak universality)
# ------------------------------------------------------------------------------

# At the ASV level this takes 2-3 min.
rug_obj <- summarize_Sigmas(output_dir)

# Subset to selected hosts
subset_idx <- rug_obj$hosts %in% selected_hosts
rug_subset <- rug_obj$rug[subset_idx,]
rug_host_subset <- rug_obj$hosts[subset_idx]

universalities <- apply(rug_subset, 2, calc_universality_score)
consensus_sign <- apply(rug_subset, 2, function(x) {
  sign(sum(sign(x)))
})

positive_idx <- which(consensus_sign > 0)
negative_idx <- which(consensus_sign < 0)
positive_ranks <- order(universalities[positive_idx], decreasing = TRUE)
negative_ranks <- order(universalities[negative_idx], decreasing = TRUE)

k <- 3

top_k_positive <- positive_ranks[1:k]
# universalities[positive_idx[positive_ranks[1:10]]]
top_k_negative <- negative_ranks[1:k]
# universalities[negative_idx[negative_ranks[1:10]]]
bottom_k <- order(universalities)[1:k]
# universalities[bottom_k]

# ------------------------------------------------------------------------------
#   Align and fully interpolate host series
# ------------------------------------------------------------------------------

named_pairs <- list(positive1 = positive_idx[positive_ranks[1]],
                    negative1 = negative_idx[negative_ranks[1]],
                    neutral1 = bottom_k[1],
                    positive2 = positive_idx[positive_ranks[2]],
                    negative2 = negative_idx[negative_ranks[2]],
                    neutral2 = bottom_k[2],
                    positive3 = positive_idx[positive_ranks[3]],
                    negative3 = negative_idx[negative_ranks[3]],
                    neutral3 = bottom_k[3])

for(pair_name in names(named_pairs)) {
  cat("Evaluating pair", pair_name, "...\n")
  pair <- named_pairs[[pair_name]]
  coord1 <- rug_obj$tax_idx1[pair]
  coord2 <- rug_obj$tax_idx2[pair]

  cat("Aligning and interpolating host series ...\n")
  # This takes 1-2 min.
  aligned_obj <- NULL
  for(host in selected_hosts) {
    aligned_obj_host <- get_paired_trajectories(pred_obj$dates[[host]],
                                            coord1,
                                            coord2,
                                            host,
                                            pred_obj$common_baseline,
                                            pred_obj$predictions[[host]],
                                            center = NULL,
                                            mean_only = TRUE)
    # Linearly interpolate any missing days
    coord1_days <- aligned_obj_host %>%
      filter(coord == coord1) %>%
      arrange(day) %>%
      pull(day)
    coord1_mean <- aligned_obj_host %>%
      filter(coord == coord1) %>%
      arrange(day) %>%
      pull(mean)
    coord1_days_interp <- min(coord1_days):max(coord1_days)
    coord1_mean_interp <- approx(coord1_days, coord1_mean, coord1_days_interp)$y
    coord2_days <- aligned_obj_host %>%
      filter(coord == coord2) %>%
      arrange(day) %>%
      pull(day)
    coord2_mean <- aligned_obj_host %>%
      filter(coord == coord2) %>%
      arrange(day) %>%
      pull(mean)
    coord2_days_interp <- min(coord2_days):max(coord2_days)
    coord2_mean_interp <- approx(coord2_days, coord2_mean, coord2_days_interp)$y

    # Center the series before adding them to the full data.frame
    aligned_obj_host_coord1 <- data.frame(host = host,
                                          coord = coord1,
                                          day = coord1_days_interp,
                                          mean = coord1_mean_interp)
    aligned_obj_host_coord2 <- data.frame(host = host,
                                          coord = coord2,
                                          day = coord2_days_interp,
                                          mean = coord2_mean_interp)

    if(is.null(aligned_obj)) {
      aligned_obj <- rbind(aligned_obj_host_coord1, aligned_obj_host_coord2)
    } else {
      aligned_obj <- rbind(aligned_obj,
                           rbind(aligned_obj_host_coord1, aligned_obj_host_coord2))
    }
  }

  # Filtering to days in common across retained hosts
  shared_days <- aligned_obj[aligned_obj$host == selected_hosts[1],]$day
  for(host in selected_hosts[2:length(selected_hosts)]) {
    shared_days <- intersect(shared_days,
                             aligned_obj[aligned_obj$host == host,]$day)
  }

  filtered_obj <- aligned_obj %>%
    filter(day %in% shared_days) %>%
    arrange(day)

  # Center the series after subsetting for days-in-common
  filtered_obj <- aligned_obj %>%
    filter(day %in% shared_days) %>%
    arrange(day) %>%
    group_by(host, coord) %>%
    mutate(offset = mean(mean)) %>%
    mutate(centered_series = mean - offset) %>%
    select(host, coord, day, centered_series)

  # ------------------------------------------------------------------------------
  #   Estimate within- and between-host correlation distributions
  # ------------------------------------------------------------------------------

  cat("Building within- and between-host distributions ...\n")

  # Within hosts
  within_distro <- c()
  for(host in selected_hosts) {
    series1 <- pull_series(filtered_obj, host, coord1)
    series2 <- pull_series(filtered_obj, host, coord2)
    within_distro <- c(within_distro, cor(series1, series2))
  }

  # Across hosts
  host_combos <- combn(selected_hosts, m = 2)
  between_distro_1 <- sapply(1:ncol(host_combos), function(i) {
    cor(pull_series(filtered_obj, host_combos[1,i], coord1),
        pull_series(filtered_obj, host_combos[2,i], coord1))
  })
  between_distro_2 <- sapply(1:ncol(host_combos), function(i) {
    cor(pull_series(filtered_obj, host_combos[1,i], coord2),
        pull_series(filtered_obj, host_combos[2,i], coord2))
  })

  plot_df <- data.frame(correlation = within_distro,
                        type = "within hosts")
  plot_df <- rbind(plot_df,
                   data.frame(correlation = between_distro_1,
                              type = "between hosts (1)"))
  plot_df <- rbind(plot_df,
                   data.frame(correlation = between_distro_2,
                              type = "between hosts (2)"))

  label1 <- get_tax_label(data$taxonomy, coord1, "clr")
  label2 <- get_tax_label(data$taxonomy, coord2, "clr")
  alpha <- 0.5
  ggplot() +
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
                        label1,
                        " x ",
                        label2,
                        "\n"))
  ggsave(file.path(plot_dir, paste0("within-between_", pair_name, ".png")),
         units = "in",
         dpi = 100,
         height = 4,
         width = 6)
}
