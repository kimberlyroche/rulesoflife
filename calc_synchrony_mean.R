library(rulesoflife)
library(tidyverse)
library(driver)

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
#   Pull interesting pairs (strong and weak universality)
# ------------------------------------------------------------------------------

# At the ASV level this takes 2-3 min.
rug_obj <- summarize_Sigmas(output_dir)

# Subset to selected hosts
subset_idx <- rug_obj$hosts %in% selected_hosts
rug_subset <- rug_obj$rug[subset_idx,]
rug_host_subset <- rug_obj$hosts[subset_idx]

universalities <- apply(rug_subset, 2, calc_universality_score)
consensus_sign <- apply(rug_subset, 2, calc_consensus_sign)

positive_idx <- which(consensus_sign > 0)
negative_idx <- which(consensus_sign < 0)
positive_ranks <- order(universalities[positive_idx], decreasing = TRUE)
negative_ranks <- order(universalities[negative_idx], decreasing = TRUE)

k <- 3

top_k_positive <- positive_ranks[1:k]
# Indexed as: universalities[positive_idx[positive_ranks[1:10]]]
top_k_negative <- negative_ranks[1:k]
# Indexed as: universalities[negative_idx[negative_ranks[1:10]]]
bottom_k <- order(universalities)[1:k]
# Indexed as: universalities[bottom_k]

# ------------------------------------------------------------------------------
#   Render all host predictions once
#
#   Code pulled from basset prediction: predict mean of Lambda
# ------------------------------------------------------------------------------

# Runtime ~10 sec. per host
predict_GP_mean <- function(output_dir, host) {
  cat("Predicting mean Lambda for host", host, "...\n")
  fit <- readRDS(file.path("output", "model_fits", output_dir, "full_posterior", paste0(host, ".rds")))

  # Observed data
  X_o <- fit$X

  # Build unobserved data set
  first_day <- min(X_o[1,])
  last_day <- max(X_o[1,])
  n_days <- last_day
  span <- seq(from = first_day, to = last_day)
  X_u <- matrix(NA, nrow(X_o), length(span))
  X_u[1,] <- span
  # Linearly interpolate the diet PCs for lack of a better option
  x <- fit$X[1,]
  for(cov_idx in 2:4) {
    y <- fit$X[cov_idx,]
    X_u[cov_idx,] <- approx(x = x, y = y, xout = span, ties = "ordered")$y
  }

  Gamma <- fit$Gamma(cbind(X_o, X_u))
  obs <- c(rep(TRUE, ncol(fit$X)), rep(FALSE, ncol(X_u)))

  # Predict Lambda
  Gamma_oo <- Gamma[obs, obs, drop=F]
  Gamma_ou <- Gamma[obs, !obs, drop=F]

  Theta_o <- fit$Theta(X_o)
  Theta_u <- fit$Theta(X_u)
  Lambda_o <- apply(fit$Lambda, c(1,2), mean)
  mean_Lambda_pred <- Theta_u + (Lambda_o-Theta_o)%*%solve(Gamma_oo, Gamma_ou)

  # convert to CLR
  mean_Lambda.clr <- clr_array(alrInv_array(mean_Lambda_pred, fit$alr_base, 1), 1)

  return(list(predictions = mean_Lambda.clr, span = span))
}

pred_filename <- file.path("output", "host_mean_predictions.rds")
if(file.exists(pred_filename)) {
  predictions <- readRDS(pred_filename)
} else {
  # Runtime: ~7 min.
  predictions <- data.frame(host = c(),
                            coord = c(),
                            day = c(),
                            clr_abundance = c())
  for(host in selected_hosts) {
    host_pred_obj <- predict_GP_mean(output_dir, host)
    pred_obj <- host_pred_obj$predictions
    host_dates <- data$metadata[data$metadata$sname == host,]$collection_date
    host_days <- round(unname(sapply(host_dates, function(x) {
      difftime(x, common_baseline, units = "days")
    }))) + 1
    # Recalculate the "days" index for this host given the new baseline
    offset <- host_days[1] - 1
    days_obj <- host_pred_obj$span + offset
    # Roll this into a data.frame

    pred_obj <- cbind(1:nrow(pred_obj), pred_obj)
    colnames(pred_obj) <- c("coord", days_obj)
    pred_long <- pivot_longer(as.data.frame(pred_obj),
                              !coord,
                              names_to = "day",
                              values_to = "clr_abundance")
    pred_long <- cbind(host = host, pred_long)
    predictions <- rbind(predictions, pred_long)
  }
  saveRDS(predictions, file = pred_filename)
}

cat("Filtering to days in common across retained hosts...\n")
shared_days <- predictions[predictions$host == selected_hosts[1],]$day
for(host in selected_hosts[2:length(selected_hosts)]) {
  shared_days <- intersect(shared_days,
                           predictions[predictions$host == host,]$day)
}

filtered_predictions <- predictions %>%
  filter(day %in% shared_days) %>%
  arrange(day)

cat("Centering the series...\n")
centered_predictions <- filtered_predictions %>%
  group_by(host, coord) %>%
  mutate(offset = mean(clr_abundance)) %>%
  # mutate(centered_clr = clr_abundance - offset) %>%
  mutate(centered_clr = scale(clr_abundance)) %>%
  select(host, coord, day, centered_clr)

# ------------------------------------------------------------------------------
#   Align and fully interpolate host series
# ------------------------------------------------------------------------------

pull_series <- function(df, sname, idx) {
  df %>%
    filter(host == sname) %>%
    arrange(day) %>%
    filter(coord == idx) %>%
    pull(centered_clr)
}

named_pairs <- list(positive1 = positive_idx[positive_ranks[1]],
                    negative1 = negative_idx[negative_ranks[1]],
                    neutral1 = bottom_k[1],
                    positive2 = positive_idx[positive_ranks[2]],
                    negative2 = negative_idx[negative_ranks[2]],
                    neutral2 = bottom_k[2],
                    positive3 = positive_idx[positive_ranks[3]],
                    negative3 = negative_idx[negative_ranks[3]],
                    neutral3 = bottom_k[3])

# Each iteration of this loop takes about 40 sec. if we look at all pairs of
# hosts. If we subset to 100 hosts, it's about 9 sec.
subset_hosts <- FALSE
for(pair_name in names(named_pairs)) {
  cat("Evaluating pair", pair_name, "...\n")
  pair <- named_pairs[[pair_name]]
  coord1 <- rug_obj$tax_idx1[pair]
  coord2 <- rug_obj$tax_idx2[pair]

  cat("Subsetting predictions to coordinates of interest...\n")
  subset_predictions <- centered_predictions %>%
    filter(coord %in% c(coord1, coord2))

  # ------------------------------------------------------------------------------
  #   Estimate within- and between-host correlation distributions
  # ------------------------------------------------------------------------------

  cat("Building within-host distribution...\n")
  within_distro <- c()
  for(host in selected_hosts) {
    series1 <- pull_series(subset_predictions, host, coord1)
    series2 <- pull_series(subset_predictions, host, coord2)
    within_distro <- c(within_distro, cor(series1, series2))
  }

  cat("Building between-host distribution...\n")
  host_combos <- combn(selected_hosts, m = 2)
  if(subset_hosts) {
    host_combos <- host_combos[,sample(1:ncol(host_combos), size = 100)]
  }
  between_distros <- sapply(1:ncol(host_combos), function(i) {
    c(cor(pull_series(subset_predictions, host_combos[1,i], coord1),
          pull_series(subset_predictions, host_combos[2,i], coord1)),
      cor(pull_series(subset_predictions, host_combos[1,i], coord2),
          pull_series(subset_predictions, host_combos[2,i], coord2)))
  })

  plot_df <- data.frame(correlation = within_distro,
                        type = "within hosts")
  plot_df <- rbind(plot_df,
                   data.frame(correlation = between_distros[1,],
                              type = "between hosts (1)"))
  plot_df <- rbind(plot_df,
                   data.frame(correlation = between_distros[2,],
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
                        "(1) ", label1,
                        "\n",
                        "(2) ", label2,
                        "\n")) #+
    # theme(plot.title = element_text(size = 12))
  ggsave(file.path(plot_dir, paste0("within-between_", pair_name, ".png")),
         units = "in",
         dpi = 100,
         height = 4.5,
         width = 6)
}
