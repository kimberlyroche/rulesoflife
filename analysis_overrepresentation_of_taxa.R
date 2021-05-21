source("path_fix.R")

library(rulesoflife)
library(tidyverse)
library(ape)

source("ggplot_fix.R")

# ------------------------------------------------------------------------------
#   Pull universality scores and get top k taxa
# ------------------------------------------------------------------------------

data <- load_data(tax_level = "ASV")
plot_dir <- check_dir(c("output", "figures"))
rug_obj <- summarize_Sigmas(output_dir = "asv_days90_diet25_scale1")
scores_df <- data.frame(score = apply(rug_obj$rug, 2, calc_universality_score),
                        sign = apply(rug_obj$rug, 2, calc_consensus_sign),
                        tax_idx1 = rug_obj$tax_idx1,
                        tax_idx2 = rug_obj$tax_idx2)
clr.means <- get_mean_clr_abundance(tax_level = "ASV")

percents <- c(5)
for(percent in percents) {
  k <- percent_to_k(percent, nrow(scores_df))
  top_pairs <- scores_df %>%
    arrange(desc(score)) %>%
    slice(1:k)

  # Get unique taxa in these interactions
  unique_tax_idx <- c(top_pairs$tax_idx1, top_pairs$tax_idx2)
  plot_df <- data.frame(taxon = unique_tax_idx) %>%
    group_by(taxon) %>%
    count() %>%
    arrange(desc(n)) %>%
    mutate(family = sapply(taxon,
                           function(x) get_tax_label(data$taxonomy, x, "clr"))) %>%
    mutate(label = paste0("Taxon ", taxon, " - ", family))
  plot_df <- plot_df %>%
    filter(n > median(plot_df$n))
  plot_df$clr_mean <- clr.means[plot_df$taxon]
  p <- ggplot(plot_df, aes(x = reorder(label, n), y = n)) +
    geom_bar(aes(fill = clr_mean), stat = "identity") +
    scale_fill_gradient2(low = "#50A3A4", mid = "#FCAF38", high = "#F95335",
                         midpoint = (max(plot_df$clr_mean) - min(plot_df$clr_mean))/2) +
    labs(x = "",
         y = paste0("counts in top ", percent, "% pairs"),
         fill = "Mean CLR\nabundance") +
    # theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    coord_flip()
  ggsave(file.path(plot_dir, paste0("representation_top_", percent, ".png")),
         p,
         units = "in",
         dpi = 100,
         height = 8,
         width = 8)
  show(p)
}

# ------------------------------------------------------------------------------
#   Moran's I - does phylogenetic distance predict coocurrence in universal
#     pairs?
# ------------------------------------------------------------------------------

# Statistic resource:
# https://rfunctions.blogspot.com/2014/02/measuring-phylogenetic-signal-in-r.html#:~:text=Moran's%20I%20is%20known%20as,instead%20of%20'spatial%20proximity'.&text=When%20I%20equals%200%2C%20species,under%20a%20brownian%20motion%20model.

## weights w[i,j] = 1/d[i,j]:

w <- 1/cophenetic(tr)
## set the diagonal w[i,i] = 0 (instead of Inf...):
diag(w) <- 0
Moran.I(x, w)
Moran.I(x, w, alt = "l")
Moran.I(x, w, alt = "g")
Moran.I(x, w, scaled = TRUE) # usualy the same
