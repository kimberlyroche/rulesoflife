# This script compares estimates of CLR correlations made using two different
# taxa as the ALR reference in the basset model: the "other" taxon and the taxon
# with the median coefficient of variation in terms of relative abundances (this
# appears to be all taxa collapsed to the class Mollicutes).

library(rulesoflife)
library(fido)
library(tidyverse)

# rug_obj <- summarize_Sigmas(output_dir = "fam_days90_diet25_scale1")
# rug_obj_alr <- summarize_Sigmas(output_dir = "fam_days90_diet25_scale1_ALRmedian/")
#
# ordering <- plot_rug(rug = rug_obj$rug)

# Build a map of taxon indices from one model (other-as-ref) to the other
# (median-as-ref).
data <- readRDS(file.path("input", "processed_family_1_20.rds"))
data_alr <- readRDS(file.path("input", "processed_family_1_20_ALRmedian.rds"))
map <- left_join(data.frame(tax = 1:nrow(data$taxonomy),
                            OTU = data$taxonomy$OTU),
                 data.frame(tax_alr = 1:nrow(data_alr$taxonomy),
                            OTU = data_alr$taxonomy$OTU),
                 by = "OTU")
map <- map %>%
  select(tax, tax_alr)

output_dir <- "fam_days90_diet25_scale1"
output_dir_alr <- "fam_days90_diet25_scale1_ALRmedian"


# Canonical ordering
# order_df <- data.frame(col_in_rug = 1:ncol(rug_obj$rug),
#                        tax1 = rug_obj$tax_idx1[ordering$col_order],
#                        tax2 = rug_obj$tax_idx2[ordering$col_order])
# order_df$OTU1 <- data$taxonomy$OTU[rug_obj$tax_idx1[ordering$col_order]]
# order_df$OTU2 <- data$taxonomy$OTU[rug_obj$tax_idx2[ordering$col_order]]
# head(order_df)
#
# alt_df <- data.frame(col_alr = 1:ncol(rug_obj_alr$rug),
#                      tax1_alr = rug_obj_alr$tax_idx1,
#                      tax2_alr = rug_obj_alr$tax_idx2)
# alt_df$OTU1 <- data_alr$taxonomy$OTU[order_df$tax1]
# alt_df$OTU2 <- data_alr$taxonomy$OTU[order_df$tax2]
#
# head(alt_df)
#
# temp <- left_join(order_df, alt_df, by = c("OTU1", "OTU2"))
# View(temp[1:10,])
#
# temp[1,]
# plot(data$counts[1,], data_alr$counts[1,])
#
#
# order_df$tax1_alr <- sapply(order_df$tax1, function(x) {
#   map %>% filter(tax == x) %>% pull(tax_alr)
# })
# order_df$tax2_alr <- sapply(order_df$tax2, function(x) {
#   map %>% filter(tax == x) %>% pull(tax_alr)
# })
#
# head(order_df)
#
# # Test agreement
# par(mfrow = c(1,2))
# idx <- sample(1:nrow(order_df), size = 1)
# plot(data$counts[order_df[idx,]$tax1,], data_alr$counts[order_df[idx,]$tax1_alr,])
# plot(data$counts[order_df[idx,]$tax2,], data_alr$counts[order_df[idx,]$tax2_alr,])
#
# # Translate this into a new column ordering based on the identities of the pairs
# alr_pairs_to_idx <- data.frame(pair = 1:length(rug_obj_alr$tax_idx1),
#                                tax1_ref = rug_obj_alr$tax_idx1,
#                                tax2_ref = rug_obj_alr$tax_idx2)
# head(alr_pairs_to_idx)
#
# map2 <- left_join(order_df, alr_pairs_to_idx, by = c("tax1_alr" = "tax1_ref",
#                                                      "tax2_alr" = "tax2_ref"))
# map2 %>%
#   filter(!is.na(pair)) %>%
#   pull(pair)

# Pull the averaged estimated correlation across CLR taxa for a random host and
# plot estimates made using ALR reference #1 ("other" category) to estimates of
# the same pair made using ALR reference #2 (median CoV taxon).
metadata <- load_data()$metadata
hosts <- unique(metadata$sname)
x = c()
y = c()
for(host in hosts) {
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

plot_df <- data.frame(x = x, y = y)
p <- ggplot(plot_df, aes(x = x, y = y)) +
  geom_point() +
  xlab("CLR correlation (ALR ref. 'other')") +
  ylab("CLR correlation (ALR ref. median)")
show(p)
plot_dir <- check_dir(c("output", "figures"))
ggsave(file.path(plot_dir, "ALR_ref_comparison.png"),
       p,
       units = "in",
       height = 4,
       width = 4)

cat("Correlation of estimates:", round(cor(x, y), 3), "\n")
cat("R^2:", round(summary(lm(y ~ x, data = plot_df))$r.squared, 3), "\n")



