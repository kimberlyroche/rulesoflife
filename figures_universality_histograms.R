source("path_fix.R")

library(rulesoflife)
library(tidyverse)

source("ggplot_fix.R")

plot_dir <- check_dir(c("output", "figures"))

map <- data.frame(level = c("phylum", "family", "ASV"),
                  short_name = c("phy", "fam", "asv"))

for(level in map$level) {
  cat("Rendering BASELINE histogram for", level, "\n")
  output_dir <- paste0(map[map$level == level,]$short_name,
                       "_days90_diet25_scale1")
  rug_obj <- summarize_Sigmas(output_dir)

  # ----------------------------------------------------------------------------
  #   Histogram of CLR correlation
  # ----------------------------------------------------------------------------

  p <- ggplot() +
    geom_histogram(data = data.frame(correlation = c(rug_obj$rug)),
                   aes(x = correlation),
                   alpha = 1,
                   color = "white")+
    labs(x = "CLR correlation")
  ggsave(file.path(plot_dir, paste0("CLR_correlation_distribution_",
                                    level,
                                    ".png")),
         plot = p,
         units = "in",
         dpi = 100,
         height = 4,
         width = 5)
  show(p)

  # ----------------------------------------------------------------------------
  #   Histogram of observed/null universality scores
  # ----------------------------------------------------------------------------

  scores <- apply(rug_obj$rug, 2, calc_universality_score)

  scrambled_rug <- t(apply(rug_obj$rug, 1, sample))
  scores_scrambled <- apply(scrambled_rug, 2, calc_universality_score)

  plot_df <- data.frame(score = c(scores, scores_scrambled),
                        type = c(rep("observed", length(scores)),
                                 rep("null", length(scores_scrambled))))
  p <- ggplot() +
    geom_histogram(data = plot_df[plot_df$type != "null",],
                   aes(x = score, fill = type),
                   alpha = 0.5,
                   color = "white") +
    geom_histogram(data = plot_df[plot_df$type == "null",],
                   aes(x = score, fill = type),
                   alpha = 0.5,
                   color = "white") +
    scale_fill_manual(values = c(null = "red", observed = "blue")) +
    labs(fill = "Distribution",
         x = "universality score")
  ggsave(file.path(plot_dir, paste0("universality_distributions_",
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

  cutoff <- quantile(scores_scrambled, probs = c(0.95))
  cat(paste0("95% cutoff at: ", round(cutoff, 3), "\n"))

  # How many hits?
  cat(paste0("Scores above 95% cutoff: ",
             round(sum(scores > cutoff) / length(scores), 3)*100,
             "%\n"))

}
