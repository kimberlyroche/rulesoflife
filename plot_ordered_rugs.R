# This script plots the "rug" heatmap over hosts and pairwise correlation over
# taxa using various row and column orderings.

source("path_fix.R")

library(rulesoflife)
library(tidyverse)

# Pull reference data and generate rug object.
data <- load_data(tax_level = "ASV")

output_dir = "asv_days90_diet25_scale1"
plot_dir <- file.path("output", "figures")

rug_obj <- summarize_Sigmas(output_dir = output_dir)

# ------------------------------------------------------------------------------
#   Column (pair) reorderings
# ------------------------------------------------------------------------------

# Column order: distance between taxon pair mean CLR abundances

# order_obj <- order_rug_cols_mean_abundance(rug_obj = rug_obj, counts = data$counts)
# temp <- plot_rug(rug_obj$rug,
#          canonical_col_order = order_obj$order,
#          canonical_row_order = NULL,
#          save_name = paste0(output_dir, "_colorder-abundance"))
#
# p <- ggplot(data.frame(x = 1:length(order_obj$difference),
#                   y = order_obj$difference[order_obj$order]),
#        aes(x = x, y = y)) +
#   geom_point(size = 0.5) +
#   xlab("sample index") +
#   ylab("mean absolute CLR difference") +
#   theme_minimal()
# save_name <- paste0(output_dir, "_colorder-abundance_differences")
# ggsave(file.path(plot_dir, paste0(save_name, ".png")),
#        p,
#        units = "in",
#        height = 2,
#        width = 8)

# ------------------------------------------------------------------------------
#   Row (host) reorderings
# ------------------------------------------------------------------------------

# Distance between baseline compositions

row_order <- order_rug_row_baseline(rug_obj = rug_obj,
                                    counts = data$counts,
                                    metadata = data$metadata)
temp <- plot_rug(rug_obj$rug,
                 canonical_col_order = NULL,
                 canonical_row_order = row_order,
                 row_labels = rug_obj$hosts[row_order],
                 save_name = paste0(output_dir, "_roworder-baseline"))

# Distance between mean diversity (slow to run)

row_order <- order_rug_row_diversity(rug_obj = rug_obj)
temp <- plot_rug(rug_obj$rug,
                 canonical_col_order = NULL,
                 canonical_row_order = row_order,
                 row_labels = rug_obj$hosts[row_order],
                 save_name = paste0(output_dir, "_roworder-diversity"))

# By social groups

row_order <- order_rug_row_group(rug_obj = rug_obj)
temp <- plot_rug(rug_obj$rug,
                 canonical_col_order = NULL,
                 canonical_row_order = row_order,
                 row_labels = rug_obj$hosts[row_order],
                 save_name = paste0(output_dir, "_roworder-group"))

# By pedigree
row_order <- order_rug_row_pedigree(rug_obj)
temp <- plot_rug(rug_obj$rug,
                 canonical_col_order = NULL,
                 canonical_row_order = row_order,
                 row_labels = rug_obj$hosts[row_order],
                 save_name = paste0(output_dir, "_roworder-pedigree"))

# ------------------------------------------------------------------------------
#   Column (pair) reorderings
# ------------------------------------------------------------------------------

# By phylogenetic distance between taxon pair

canonical_col_order <- order_rug_row_pedigree(rug_obj, data$taxonomy)

ordering <- plot_rug(rug_obj$rug,
                     canonical_col_order = canonical_col_order,
                     canonical_row_order = NULL,
                     row_labels = rug_obj$hosts,
                     save_name = paste0(output_dir, "_colorder-taxonomic"))

# Render some extra plots

plot_df <- data.frame(x = tax_distances)
ggplot(plot_df, aes(x = x)) +
  geom_histogram(color = "white") +
  labs(x = "phylogenetic distance")
ggsave(file.path(plot_dir, "distro_phylogenetic_distances.png"),
       units = "in",
       dpi = 100,
       height = 3.5,
       width = 4)

universalities <- apply(rug_obj$rug, 2, calc_universality_score)
consensus_sign <- apply(rug_obj$rug, 2, function(x) {
  sign(sum(sign(x)))
})

plot_df <- data.frame(universality = universalities[canonical_col_order],
                      universality_sign = consensus_sign[canonical_col_order],
                      tax_distance = tax_distances[canonical_col_order])
plot_df$universality_sign <- factor(plot_df$universality_sign)

ggplot(plot_df, aes(x = tax_distance, y = universality, color = universality_sign)) +
  geom_point(alpha = 0.66) +
  scale_color_manual(values = c("#5680e3", "#bbbbbb", "#e35656")) +
  labs(x = "phylogenetic distance",
       color = "Consensus\nsign")
ggsave(file.path(plot_dir, "phylogenetic_distance_vs_universality.png"),
       units = "in",
       dpi = 100,
       height = 5,
       width = 7)
