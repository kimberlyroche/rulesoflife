source("path_fix.R")

library(rulesoflife)
library(tidyverse)

source("ggplot_fix.R")

data <- load_data(tax_level = "ASV")
output_dir <- "asv_days90_diet25_scale1"
rug_obj <- summarize_Sigmas(output_dir)

scores <- apply(rug_obj$rug, 2, calc_universality_score)
consensus_signs <- apply(rug_obj$rug, 2, calc_consensus_sign)

# Supply a cutoff for "universality" scores.
# These can be derived from `figures_universality_histograms.R`.
#   phylum: ~0.164
#   family: ~0.139
#      ASV: ~0.129
cutoff <- 0.15

# ----------------------------------------------------------------------------
#   ID the weakest passing "significant" pair
# ----------------------------------------------------------------------------

# Utility function
plot_trajectories_from_idx <- function(idx, rug_obj, data, file_tag) {
  tax_idx1 <- rug_obj$tax_idx1[idx]
  tax_idx2 <- rug_obj$tax_idx2[idx]
  tax_label1 <- get_tax_label(data$taxonomy, tax_idx1, "clr")
  tax_label2 <- get_tax_label(data$taxonomy, tax_idx2, "clr")

  # Note: this takes ~10 min. at the ASV level!
  p <- plot_aligned_trajectories(output_dir,
                            tax_idx1 = tax_idx1,
                            tax_idx2 = tax_idx2,
                            tax_label1 = tax_label1,
                            tax_label2 = tax_label2,
                            metadata = data$metadata,
                            return_plot = FALSE,
                            file_tag = file_tag)
}

scores_df <- data.frame(index = 1:length(scores),
                        score = scores,
                        sign = consensus_signs) %>%
  mutate(above_cutoff = ifelse(score > cutoff,
                               TRUE,
                               FALSE)) %>%
  filter(above_cutoff == TRUE)

low_passing <- scores_df %>%
  arrange(score) %>%
  slice(1) %>%
  pull(index)

highest_positive <- scores_df %>%
  filter(sign > 0) %>%
  arrange(desc(score)) %>%
  slice(1) %>%
  pull(index)

highest_negative <- scores_df %>%
  filter(sign < 0) %>%
  arrange(desc(score)) %>%
  slice(1) %>%
  pull(index)

plot_trajectories_from_idx(low_passing, rug_obj, data, "weak")
cat(paste0("Plotting selected pair: ", tax_label1, " x ", tax_label2,
           "; sign: ",
           scores_df$sign[low_passing],
           "\n"))

plot_trajectories_from_idx(highest_positive, rug_obj, data, "positive")
plot_trajectories_from_idx(highest_negative, rug_obj, data, "negative")
