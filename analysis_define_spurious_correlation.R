# NEED TO TEST THIS!

source("path_fix.R")

library(rulesoflife)

source("ggplot_fix.R")

plot_dir <- check_dir(c("output", "figures"))

map <- data.frame(level = c("phylum", "family", "ASV"),
                  short_name = c("phy", "fam", "asv"))

for(level in map$level) {
  output_dir <- paste0(map[map$level == level,]$short_name,
                       "_days90_diet25_scale1")
  cat("Evaluating", level, "\n")
  rug_obj <- summarize_Sigmas(output_dir = output_dir)
  rug <- rug_obj$rug

  # ----------------------------------------------------------------------------
  #   Visualize observed vs. null distributions
  # ----------------------------------------------------------------------------

  # Permuted data
  rug_obj_p <- summarize_Sigmas(output_dir = paste0(output_dir, "_scrambled"))
  rug_p <- rug_obj_p$rug

  vector_observed <- c(rug)
  vector_null <- c(rug_p)
  plot_df <- data.frame(correlation = c(vector_observed, vector_null),
                        type = c(rep("observed", length(vector_observed)),
                                 rep("null", length(vector_null))))
  p <- ggplot() +
    geom_histogram(data = plot_df[plot_df$type != "null",],
                   aes(x = correlation, fill = type),
                   alpha = 0.5,
                   color = "white") +
    geom_histogram(data = plot_df[plot_df$type == "null",],
                   aes(x = correlation, fill = type),
                   alpha = 0.5,
                   color = "white") +
    scale_fill_manual(values = c(null = "red", observed = "blue")) +
    labs(fill = "Distribution",
         x = "universality score")
  ggsave(file.path(plot_dir, paste0("correlation_distributions_",
                                    level,
                                    ".png")),
         plot = p,
         units = "in",
         dpi = 100,
         height = 4,
         width = 5)
  show(p)

  # ----------------------------------------------------------------------------
  #   Define a 95% cutoff
  # ----------------------------------------------------------------------------

  cutoff <- quantile(vector_null, probs = c(0.95))
  cat(paste0("95% cutoff at: ", round(cutoff, 3), "\n"))

  # How many hits?
  cat(paste0("Values above 95% cutoff: ",
             round(sum(vector_observed > cutoff) / length(vector_observed), 3)*100,
             "%\n"))

  # ----------------------------------------------------------------------------
  #   Define a cutoff in terms of standard deviations from the mean
  # ----------------------------------------------------------------------------

  lower <- mean(vector_null) - 2*sd(vector_null)
  upper <- mean(vector_null) + 2*sd(vector_null)

  cat("Lower, upper spurious correlation bounds:",
      round(lower, 3),
      ",",
      round(upper, 3),
      "\n")

  cat("Proportion of observed correlations exceeding spurious threshold:",
      round(sum(vector_observed < lower | vector_observed > upper) / length(vector_observed), 3)*100,
      "\n")
}
