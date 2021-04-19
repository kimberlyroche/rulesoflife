rm(list = ls())

library(tidyverse)
library(rulesoflife)
library(fido)

# Load
data <- load_data(tax_level = "family")

# Fit sample model
# fit <- fit_GP(sname = "VAP", counts = data$counts, metadata = data$metadata, point_est = FALSE)
# fit <- fit_DLM(sname = "VAP", counts = data$counts, metadata = data$metadata, point_est = FALSE, smoothed = TRUE)

# plot_trajectory(fit,
#                 coord = coord,
#                 coord_label = get_tax_label(data$taxonomy, 1, "CLR"),
#                 show_observations = TRUE,
#                 save_file = FALSE)

plot_aligned_trajectories(host_list = c("ACA", "ALE", "CAI", "COB", "COO", "DAG"),
                          tax_idx1 = 1,
                          tax_idx2 = 2,
                          metadata = data$metadata,
                          save_file = TRUE)

fit <- fit_GP(sname = "VAP",
              counts = data$counts,
              metadata = data$metadata,
              point_est = FALSE,
              output_dir = "MAP_days30",
              diet_weight = 0,
              days_to_min_autocorrelation = 30)

sensitivity_sweep(output_dir_list = c("MAP_diet0", "MAP_diet25", "MAP_diet50"))
