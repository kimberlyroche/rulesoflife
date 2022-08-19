library(rulesoflife)
library(fido)

# ------------------------------------------------------------------------------
#   Distribution of days between samples
# ------------------------------------------------------------------------------

data <- load_data(level = "ASV")
md <- data$metadata

# Get within-host daily differences
md %<>%
  group_by(sname) %>%
  mutate(days_indexed = as.numeric(difftime(collection_date, min(collection_date), units = "day"))) %>%
  arrange(days_indexed) %>%
  mutate(days_diff = c(NA, diff(days_indexed))) %>%
  ungroup()

x <- md$days_diff
x <- x[complete.cases(x)]

# IQR of daily differences
quantile(x, probs = c(0.25, 0.5, 0.75))

# ------------------------------------------------------------------------------
#   Credible intervals for key taxon-taxon correlations
# ------------------------------------------------------------------------------

model_dir <- "asv_days90_diet25_scale1"
rug_obj <- summarize_Sigmas(output_dir = model_dir)

# Strongly positive pair
asv1 <- 2
asv2 <- 3

# Strongly negative pair
# asv1 <- 25
# asv2 <- 107

combo_idx <- which(sapply(1:ncol(rug_obj$rug), function(i) rug_obj$tax_idx1[i] == asv1 & rug_obj$tax_idx2[i] == asv2))

sd(rug_obj$rug[,combo_idx])




