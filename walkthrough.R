source("path_fix.R")

library(tidyverse)
library(rulesoflife)
library(fido)

# Load
data <- load_data(tax_level = "family")

# Fit sample model
fit <- fit_GP(sname = "ALE",
              counts = data$counts,
              metadata = data$metadata,
              output_dir = "no_diet",
              MAP = TRUE,
              diet_weight = 0,
              days_to_min_autocorrelation = 30,
              var_scale_taxa = 2,
              var_scale_samples = 2)

sname <- "ALE"
output_dir <- "no_diet"

fit <- readRDS(paste0("output/model_fits/", output_dir, "/MAP/", sname, ".rds"))
fit <- readRDS(paste0("output/model_fits/", output_dir, "/full_posterior/", sname, ".rds"))
fit <- to_clr(fit)
Sigma <- fit$Sigma
for(i in 1:fit$iter) {
  Sigma[,,i] <- cov2cor(Sigma[,,i])
}
Sigma <- apply(Sigma, c(1,2), mean)
plot_kernel_or_cov_matrix(Sigma, save_name = NULL)

fit_pred_obj <- predict_trajectory(sname, output_dir)
plot_trajectory(sname, output_dir, coord = 1, coord_label = NULL,
                            show_observations = TRUE, save_file = FALSE)

rug_obj <- summarize_Sigmas(output_dir = "diet_50", use_proportionality = FALSE)
ordering <- plot_rug(rug_obj$rug, canonical_col_order = NULL,
                     canonical_row_order = NULL, save_name = NULL)

plot_aligned_trajectories(host_list = c("ALE", "VAP"),
                          output_dir = "no_diet",
                          tax_idx1 = 1,
                          tax_idx2 = 2,
                          metadata = data$metadata,
                          save_file = FALSE)

plot_clr_vs_diet(sname, tax_idx = 1, data$counts, data$metadata)

sensitivity_sweep(output_dir_list = c("no_diet", "diet_50"))

# Downsampling plot
data <- readRDS("output/simdata_results.rds")
View(data)
ggplot(data, aes(x = factor(sample_no), y = matched_CLR_sign)) +
  geom_boxplot() +
  xlab("sample number") +
  ylab("proportion estimates matched for CLR sign")
ggsave()













