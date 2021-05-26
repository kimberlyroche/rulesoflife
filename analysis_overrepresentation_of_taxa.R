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

library(ape)

# Statistic resource:
# https://rfunctions.blogspot.com/2014/02/measuring-phylogenetic-signal-in-r.html#:~:text=Moran's%20I%20is%20known%20as,instead%20of%20'spatial%20proximity'.&text=When%20I%20equals%200%2C%20species,under%20a%20brownian%20motion%20model.

d_mat <- rug_phylogenetic_distances(rug_obj, data$taxonomy, as_matrix = TRUE)
plot_kernel_or_cov_matrix(d_mat) +
  labs(x = "ASV 1",
       y = "ASV 2",
       fill = "phylogenetic\ndistance")

d_vec <- c(d_mat)
weight_vec <- sapply(d_vec, function(d) 1/d)
weight_mat <- matrix(weight_vec, nrow(d_mat), ncol(d_mat), byrow = FALSE)
diag(weight_mat) <- 0

# Scale the weight matrix
weight_mat <- weight_mat / max(weight_mat)

plot_kernel_or_cov_matrix(weight_mat) +
  labs(x = "ASV 1",
       y = "ASV 2",
       fill = "1 / phylogenetic\ndistance")

# From the 'ape' package example
## weights w[i,j] = 1/d[i,j]:
# w <- 1/cophenetic(tr)
# diag(w) <- 0

# We need to consolidate associations within an ASV, since x must have the same
# dimesion as the weight matrix. This is why Johannes' analysis didn't yield
# much: universality averaged over lots of pairs (for a given ASV) is pretty
# weak usually.

# ------------------------------------------------------------------------------
#   Unsigned Moran's I
# ------------------------------------------------------------------------------

scores_tax <- calc_universality_score_taxon(rug_obj, collapse = "unsigned")
combined_scores <- as.data.frame(scores_tax %>% arrange(tax1))

# Not all taxa are present after filtering for sign. Chop down the weight matrix
# accordingly.
include_tax <- intersect(1:nrow(weight_mat), combined_scores$tax1)

# Exclude the "other" category
Moran.I(combined_scores$mean_score,
        weight_mat[include_tax, include_tax],
        scaled = TRUE)

# ------------------------------------------------------------------------------
#   Signed Moran's I
# ------------------------------------------------------------------------------

scores_tax <- calc_universality_score_taxon(rug_obj, collapse = "signed")
pos_scores <- as.data.frame(scores_tax %>% filter(sign > 0) %>% arrange(tax1))
neg_scores <- as.data.frame(scores_tax %>% filter(sign < 0) %>% arrange(tax1))

# Not all taxa are present after filtering for sign. Chop down the weight matrix
# accordingly.
include_tax <- intersect(1:nrow(weight_mat), pos_scores$tax1)

# Exclude the "other" category
Moran.I(pos_scores$mean_score,
        weight_mat[include_tax, include_tax],
        scaled = TRUE)

include_tax <- intersect(1:nrow(weight_mat), neg_scores$tax1)

# Exclude the "other" category
Moran.I(neg_scores$mean_score,
        weight_mat[include_tax, include_tax],
        scaled = TRUE)

# ------------------------------------------------------------------------------
#   Omitting extremely closely related species
# ------------------------------------------------------------------------------

closely_related <- which(d_mat < quantile(d_mat, probs = c(0.008)),
                         arr.ind = TRUE)
closely_related <- closely_related[apply(closely_related,
                                         1,
                                         function(x) x[1] != x[2]),]

omit_taxa <- unname(apply(closely_related, 1, min))

pos_scores <- as.data.frame(scores_tax %>%
                              filter(sign > 0) %>%
                              arrange(tax1) %>%
                              filter(!(tax1 %in% omit_taxa)))
# Not all taxa are present after filtering for sign. Chop down the weight matrix
# accordingly.
include_tax <- intersect(1:nrow(weight_mat), pos_scores$tax1)

# Exclude the "other" category
Moran.I(pos_scores$mean_score,
        weight_mat[include_tax, include_tax],
        scaled = TRUE)

# ------------------------------------------------------------------------------
#   Box plots
# ------------------------------------------------------------------------------

pos_scores_tax <- pos_scores
pos_scores_tax$family <- data$taxonomy[pos_scores_tax$tax1,6]
pos_scores_tax <- pos_scores_tax %>%
  group_by(family) %>%
  mutate(n = n()) %>%
  filter(n >= 5) %>%
  filter(!is.na(family))

ggplot(pos_scores_tax, aes(x = family, y = mean_score)) +
  geom_boxplot() +
  geom_jitter(width = 0.1) +
  labs(y = "average ASV universality score")
