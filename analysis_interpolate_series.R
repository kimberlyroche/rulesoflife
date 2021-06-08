source("path_fix.R")

library(rulesoflife)
library(tidyverse)
library(driver)

source("ggplot_fix.R")

plot_dir <- check_dir(c("output", "figures"))
output_dir <- "asv_days90_diet25_scale1"

data <- load_data(tax_level = "ASV")
metadata <- data$metadata

# Get max and min day per host, relative to a common baseline
common_baseline <- min(metadata$collection_date)

host_filename <- file.path("output", "overlapped_hosts.rds")
if(!file.exists(host_filename)) {
  stop(paste0("File not found: ", host_filename, "\n"))
}
selected_hosts <- readRDS(host_filename)

pred_filename <- file.path("output", "host_mean_predictions.rds")
pred_adj_filename <- file.path("output", "host_mean_predictions_adjusted.rds")
if(!file.exists(pred_adj_filename)) {
  if(!file.exists(pred_filename)) {
    # Runtime: ~7 min.
    predictions <- data.frame(host = c(),
                              coord = c(),
                              day = c(),
                              clr_abundance = c())
    for(host in selected_hosts) {
      cat(paste0("Processing host ", host, "...\n"))
      host_pred_obj <- predict_GP_mean(output_dir = output_dir, host = host,
                                       interpolation = "none")
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
  } else {
    predictions <- readRDS(pred_filename)
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
    mutate(centered_clr = clr_abundance - offset) %>%
    # mutate(centered_clr = scale(clr_abundance)) %>%
    select(host, coord, day, centered_clr)

  saveRDS(centered_predictions, file = pred_adj_filename)
}
