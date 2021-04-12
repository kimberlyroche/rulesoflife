rm(list = ls())

library(rulesoflife)
library(fido)

# Load
data <- load_data(tax_level = "family")

# Chop down VAP's series for DLM test
VAP_idx <- which(data$metadata$sname == "VAP")
dates <- data$metadata$collection_date[VAP_idx]
# unname(sapply(dates, function(x) difftime(x, min(dates), units = "days")))
VAP_idx <- VAP_idx[4:15]

data$counts <- data$counts[,VAP_idx]
data$metadata <- data$metadata[VAP_idx,]

counts <- data$counts
metadata <- data$metadata
sname <- "VAP"
point_est <- FALSE
smoothed <- FALSE

# Fit sample model
# fit <- fit_GP(sname = "VAP", counts = data$counts, metadata = data$metadata, point_est = FALSE)
fit <- fit_DLM(sname = "VAP", counts = data$counts, metadata = data$metadata, point_est = FALSE, smoothed = TRUE)

str(fit)

coord <- 1

fit <- readRDS("output/GP_fits/VAP.rds")
plot_trajectory(fit,
                coord = coord,
                coord_label = get_tax_label(data$taxonomy, 1, "CLR"),
                show_observations = TRUE,
                save_file = FALSE)

fit <- readRDS("output/DLM_fits/VAP.rds")
plot_trajectory(fit,
                coord = 1,
                coord_label = get_tax_label(data$taxonomy, 1, "ALR"),
                show_observations = TRUE,
                save_file = FALSE)
