source("path_fix.R")

library(rulesoflife)
library(tidyverse)
library(driver)
library(optparse)

source("ggplot_fix.R")

option_list = list(
  make_option(c("--start"), type = "numeric", default = NULL,
              help = "taxon pair index at which to begin", metavar = "numeric"),
  make_option(c("--end"), type = "numeric", default = NULL,
              help = "taxon pair index at which to end", metavar = "numeric")
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

if(is.null(opt$start) | is.null(opt$end)) {
  stop("Start and end indices must be provided!")
}

plot_dir <- check_dir(c("output", "figures"))
output_dir <- "asv_days90_diet25_scale1"
rug_filename <- file.path("output", "rug_asv.rds")
if(file.exists(rug_filename)) {
  rug_obj <- readRDS(rug_filename)
} else {
  rug_obj <- summarize_Sigmas(output_dir)
  saveRDS(rug_obj, file = rug_filename)
}
data <- load_data(tax_level = "ASV")

if(opt$start < 1 | opt$start > opt$end) {
  stop("Invalid start and end indices!")
}
start <- opt$start
end <- opt$end

# ------------------------------------------------------------------------------
#   This script assumes the `host_mean_predictions_adjusted.rds` file already
#   exists. This file has aligned, centered series for all taxa in all (37)
#   selected hosts.
# ------------------------------------------------------------------------------

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

# Subset to selected hosts
subset_idx <- rug_obj$hosts %in% selected_hosts
rug_subset <- rug_obj$rug[subset_idx,]
rug_host_subset <- rug_obj$hosts[subset_idx]

universalities <- apply(rug_subset, 2, calc_universality_score)
# consensus_sign <- apply(rug_subset, 2, calc_consensus_sign)

subset_hosts <- TRUE
distro_df <- NULL
for(pair in start:end) {
  temp_df <- build_between_within_distributions(pair,
                                                coord1 = rug_obj$tax_idx1[pair],
                                                coord2 = rug_obj$tax_idx2[pair],
                                                predictions = predictions %>%
                                                  filter(coord %in% c(coord1, coord2)),
                                                selected_hosts = selected_hosts)

  if(is.null(distro_df)) {
    distro_df <- temp_df
  } else {
    distro_df <- bind_rows(distro_df, temp_df)
  }
}

saveRDS(distro_df,
        file = file.path("output",
                         paste0("within_between_distros_",
                                start,
                                "-",
                                end,
                                ".rds")))
