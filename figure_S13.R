source("path_fix.R")

library(rulesoflife)
library(fido)
library(tidyverse)
library(dplyr)

# ------------------------------------------------------------------------------
#
#   Supplemental Figure S13 - reproducibility of estimates using two different
#                             ALR references
#
# ------------------------------------------------------------------------------

save_fn <- file.path("output", "alr_sanity_check.rds")
if(!file.exists(save_fn)) {
  data <- readRDS(file.path("input", "processed_ASV_1_20.rds"))
  data_alr <- readRDS(file.path("input", "processed_ASV_1_20_ALRmedian.rds"))
  map <- left_join(data.frame(tax = 1:nrow(data$taxonomy),
                              OTU = data$taxonomy$OTU),
                   data.frame(tax_alr = 1:nrow(data_alr$taxonomy),
                              OTU = data_alr$taxonomy$OTU),
                   by = "OTU")
  map <- map %>%
    select(tax, tax_alr)

  output_dir <- "asv_days90_diet25_scale1"
  output_dir_alr <- "asv_days90_diet25_scale1_ALRmedian"

  # Pull the averaged estimated correlation across CLR taxa for a random host and
  # plot estimates made using ALR reference #1 ("other" category) to estimates of
  # the same pair made using ALR reference #2 (median CoV taxon).
  metadata <- load_data()$metadata
  hosts <- unique(metadata$sname)
  x = c()
  y = c()
  for(host in hosts) {
    cat(paste0("Parsing data for host ", host, "\n"))
    fit <- readRDS(file.path("output",
                             "model_fits",
                             output_dir,
                             "full_posterior",
                             paste0(host, ".rds")))
    fit <- to_clr(fit)
    fit_alr <- readRDS(file.path("output",
                                 "model_fits",
                                 output_dir_alr,
                                 "full_posterior",
                                 paste0(host, ".rds")))
    fit_alr <- to_clr(fit_alr)

    # Convert estimated covariance to correlation and average over posterior
    # samples.
    Sigma <- fit$Sigma
    Sigma_alr <- fit_alr$Sigma
    for(i in 1:fit$iter) {
      Sigma[,,i] <- cov2cor(Sigma[,,i])
      Sigma_alr[,,i] <- cov2cor(Sigma_alr[,,i])
    }
    Sigma <- apply(Sigma, c(1,2), mean)
    Sigma_alr <- apply(Sigma_alr, c(1,2), mean)

    # For this host, grab 100 randomly sampled pairs of CLR taxa and pull
    # estimates from both model paramterizations.
    for(i in 1:100) {
      idx <- sample(1:nrow(map), size = 2)
      map1 <- map[idx[1],]
      map2 <- map[idx[2],]
      x <- c(x, Sigma[map1$tax,map2$tax])
      y <- c(y, Sigma_alr[map1$tax_alr,map2$tax_alr])
    }
  }
  saveRDS(list(alr_other = x, alr_median = y), save_fn)
} else {
  save_obj <- readRDS(save_fn)
  alr_other <- save_obj$alr_other
  alr_median <- save_obj$alr_median
}

plot_df <- data.frame(x = alr_other, y = alr_median)
p <- ggplot(plot_df, aes(x = x, y = y)) +
  geom_point(size = 2, shape = 21, fill = "#888888") +
  xlab("CLR correlation (ALR ref. 'other')") +
  ylab("CLR correlation (ALR ref. median)") +
  theme_bw()

ggsave(file.path("output", "figures", "S13.png"),
       p,
       units = "in",
       height = 4,
       width = 4)

cat("Correlation of estimates:", round(cor(x, y), 3), "\n")
cat("R^2:", round(summary(lm(y ~ x, data = plot_df))$r.squared, 3), "\n")
