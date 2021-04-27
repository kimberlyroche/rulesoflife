# This script calculates the mean and median absolute CLR correlation across
# taxa in the observed and permuted data at each taxonomic level.
#
# It also calculates an approximate 95% interval associated with spurious
# correlations.

source("path_fix.R")

library(rulesoflife)

output_dirs <- c("phy_days90_diet25_scale1",
                 "fam_days90_diet25_scale1",
                 "asv_days90_diet25_scale1")

for(output_dir in output_dirs) {
  cat("Evaluating", output_dir, "\n")
  rug_obj <- summarize_Sigmas(output_dir = output_dir)
  rug <- rug_obj$rug

  cat("Mean absolute CLR correlation:", round(mean(abs(rug)), 3), "\n")
  cat("Median absolute CLR correlation:", round(median(abs(rug)), 3), "\n")

  # Permuted data
  rug_obj_p <- summarize_Sigmas(output_dir = paste0(output_dir, "_scrambled"))
  rug_p <- rug_obj_p$rug

  lower <- mean(rug_p) - 2*sd(rug_p)
  upper <- mean(rug_p) + 2*sd(rug_p)

  cat("Lower, upper spurious correlation bounds:",
      round(lower, 3),
      ",",
      round(upper, 3),
      "\n")

  cat("Proportion of observed correlations exceeding spurious threshold:",
      round(sum(rug < lower | rug > upper) / (nrow(rug)*ncol(rug)), 3)*100,
      "\n")
}
