rm(list = ls())

library(rulesoflife)

# Load
data <- load_data(tax_level = "family")

# Fit sample model
# fit <- fit_GP(sname = "VAP", counts = data$counts, metadata = data$metadata, point_est = FALSE)

fit <- readRDS("output/GP_fits/VAP.rds")
fit$sname

predictions <- predict_trajectory(fit, resolution = 10)

taxon <- 1
taxon_label <- paste0("family ", data$taxonomy[taxon,]$family)
plot_trajectory(fit, predictions, taxon, taxon_label = taxon_label,
                            show_observations = TRUE, save_file = TRUE)
