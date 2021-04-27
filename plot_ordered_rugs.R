# This script plots the "rug" heatmap over hosts and pairwise correlation over
# taxa using various row and column orderings

source("path_fix.R")

library(rulesoflife)
library(ggplot2)

# Pull reference data and generate rug object.
data <- load_data(tax_level = "family")
output_dir <- "fam_days90_diet25_scale1"
rug_obj <- summarize_Sigmas(output_dir = output_dir)

# Column order: distance between taxon pair mean CLR abundances
order_obj <- order_rug_cols_mean_abundance(rug_obj = rug_obj, counts = data$counts)
temp <- plot_rug(rug_obj$rug,
         canonical_col_order = order_obj$order,
         canonical_row_order = NULL,
         save_name = paste0(output_dir, "_colorder-abundance"))

p <- ggplot(data.frame(x = 1:length(order_obj$difference),
                  y = order_obj$difference[order_obj$order]),
       aes(x = x, y = y)) +
  geom_point(size = 0.5) +
  xlab("sample index") +
  ylab("mean absolute CLR difference") +
  theme_minimal()
save_name <- paste0(output_dir, "_colorder-abundance_differences")
plot_dir <- check_dir(c("output", "figures"))
ggsave(file.path(plot_dir, paste0(save_name, ".png")),
       p,
       units = "in",
       height = 2,
       width = 8)

# Row order: distance between baselines

row_order <- order_rug_row_baseline(rug_obj = rug_obj,
                                    counts = data$counts,
                                    metadata = data$metadata)
temp <- plot_rug(rug_obj$rug,
                 canonical_col_order = NULL,
                 canonical_row_order = row_order,
                 row_labels = rug_obj$hosts[row_order],
                 save_name = paste0(output_dir, "_roworder-baseline"))

# Row order: distance between mean diversity (slow to run)

row_order <- order_rug_row_diversity(rug_obj = rug_obj)
temp <- plot_rug(rug_obj$rug,
                 canonical_col_order = NULL,
                 canonical_row_order = row_order,
                 row_labels = rug_obj$hosts[row_order],
                 save_name = paste0(output_dir, "_roworder-diversity"))

# Row order: in order of social groups

row_order <- order_rug_row_group(rug_obj = rug_obj)
temp <- plot_rug(rug_obj$rug,
                 canonical_col_order = NULL,
                 canonical_row_order = row_order,
                 row_labels = rug_obj$hosts[row_order],
                 save_name = paste0(output_dir, "_roworder-group"))







