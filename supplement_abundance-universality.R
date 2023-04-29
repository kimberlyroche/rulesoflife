source("path_fix.R")

library(tidyverse)
library(rulesoflife)
library(driver)
library(RColorBrewer)
library(ggridges)
library(cowplot)

# Pull top 2.5% most universal pairs
data <- load_data(tax_level = "ASV")
rug_asv <- summarize_Sigmas(output_dir = "asv_days90_diet25_scale1")
filtered_pairs <- filter_joint_zeros(data$counts, threshold_and = 0.05, threshold_or = 0.5)
scores <- apply(rug_asv$rug[,filtered_pairs$threshold], 2, calc_universality_score)
scores_piecewise <- apply(rug_asv$rug[,filtered_pairs$threshold],
                          2,
                          function(x) calc_universality_score(x = x, return_pieces = TRUE))

# Get average CLR abundance for all taxa
clr.counts <- clr_array(data$counts + 0.5, parts = 1)
clr.means <- rowMeans(clr.counts)

plot_df <- data.frame(pair_idx = 1:sum(filtered_pairs$threshold),
                      tax_idx1 = rug_asv$tax_idx1[filtered_pairs$threshold],
                      tax_idx2 = rug_asv$tax_idx2[filtered_pairs$threshold],
                      median_r = apply(rug_asv$rug[,filtered_pairs$threshold], 2, median),
                      score = scores,
                      prop_agree = scores_piecewise[1,],
                      mcs = scores_piecewise[2,])
plot_df$mean1 <- clr.means[plot_df$tax_idx1]
plot_df$mean2 <- clr.means[plot_df$tax_idx2]

threshold <- quantile(scores, probs = c(0.975))
plot_df$top <- sapply(scores, function(x) x > threshold)

# Rearrange taxa 1/2 to be lower-higher
for(i in 1:nrow(plot_df)) {
  if(plot_df$mean2[i] < plot_df$mean1[i]) {
    t1 <- plot_df$tax_idx1[i]
    plot_df$tax_idx1[i] <- plot_df$tax_idx2[i]
    plot_df$tax_idx2[i] <- t1
    m1 <- plot_df$mean1[i]
    plot_df$mean1[i] <- plot_df$mean2[i]
    plot_df$mean2[i] <- m1
  }
}

# ------------------------------------------------------------------------------
#
#   Does average abundance of either partner microbe predict the pair's
#   correlation?
#
# ------------------------------------------------------------------------------

# Scale everybody so we can report r
plot_df_scaled <- plot_df %<>%
  mutate(median_r = scale(median_r),
         mean1 = scale(mean1),
         mean2 = scale(mean2))
cor(plot_df_scaled$mean1, plot_df_scaled$median_r)
cor(plot_df_scaled$mean2, plot_df_scaled$median_r)

fit <- lm(median_r ~ mean1 + mean2, data = plot_df_scaled)
cat(paste0("Mean CLR abundance of partner 1 is positively associated with median correlation strength: ",
           "beta = ", round(coef(summary(fit))[2,1], 3), ", p-value = ", round(coef(summary(fit))[2,4], 3), "\n"))
total_var <- var(scale(plot_df$median_r))
var_expl <- 1 - var(fit$residuals)/total_var
cat(paste0("\tMean 1, r = ", round(coef(summary(fit))[2,1], 3), "\n"))
cat(paste0("\tMean 2, r = ", round(coef(summary(fit))[3,1], 3), "\n"))
cat(paste0("\tPercent variance explained: ", round(var_expl*100,3), "%\n"))

# ------------------------------------------------------------------------------
#
#   Does average abundance of either partner microbe predict the pair's
#   correlation strength?
#
# ------------------------------------------------------------------------------

# fit <- lm(mcs ~ mean1*mean2, data = plot_df)
# cat(paste0("Mean CLR abundance of partner 1 is positively associated with median correlation strength: ",
#            "beta = ", round(coef(summary(fit))[2,1], 3), ", p-value = ", round(coef(summary(fit))[2,4], 3), "\n"))
# var_expl <- (var(plot_df$mcs) - var(fit$residuals))/var(plot_df$mcs)
# cat(paste0("\tPercent variance explained: ", round(var_expl*100,3), "%\n"))





