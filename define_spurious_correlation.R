# This script calculates the mean and median absolute CLR correlation across
# taxa in the observed and permuted data at each taxonomic level.
#
# It also calculates an approximate 95% interval associated with spurious
# correlations.

library(rulesoflife)

rug_obj <- summarize_Sigmas(output_dir = "fam_days90_diet25_scale1")
rug <- rug_obj$rug

cat("Mean absolute CLR correlation:", round(mean(abs(rug)), 3), "\n")
cat("Median absolute CLR correlation:", round(median(abs(rug)), 3), "\n")

# Permuted data
rug_obj <- summarize_Sigmas(output_dir = "fam_days90_diet25_scale1")
rug <- rug_obj$rug

cat("Lower, upper spurious correlation bounds:",
    round(mean(rug) - 2*sd(rug), 3),
    ",",
    round(mean(rug) + 2*sd(rug), 3),
    "\n")
