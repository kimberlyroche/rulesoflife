source("path_fix.R")

library(rulesoflife)
library(tidyverse)

source("ggplot_fix.R")

# ------------------------------------------------------------------------------
#   Baseline and scrambled rugs
# ------------------------------------------------------------------------------

map <- data.frame(level = c("phylum", "family", "ASV"),
                  short_name = c("phy", "fam", "asv"))

for(level in map$level) {
  # Baseline parameterization
  cat("Rendering BASELINE parameterization for", level, "\n")
  output_dir <- paste0(map[map$level == level,]$short_name,
                       "_days90_diet25_scale1")
  rug_obj <- summarize_Sigmas(output_dir = output_dir)
  ordering <- plot_rug(rug = rug_obj$rug,
                       row_labels = rug_obj$hosts,
                       save_name = output_dir)

  cat("Rendering on COMPOSITION for", level, "\n")
  data <- load_data(tax_level = level)
  order_obj <- order_rug_row_baseline(rug_obj = rug_obj,
                                      counts = data$counts,
                                      metadata = data$metadata)
  row_order <- order_obj$order
  discard <- plot_rug(rug = rug_obj$rug,
                      canonical_col_order = ordering$col_order,
                      canonical_row_order = row_order,
                      row_labels = rug_obj$hosts,
                      cluster_obj = order_obj$hc,
                      save_name = paste0(output_dir, "_roworder-baseline"))

  # Scrambled rug
  cat("Rendering SCRAMBLED parameterization for", level, "\n")
  output_dir <- paste0(map[map$level == level,]$short_name,
                       "_days90_diet25_scale1_scrambled")
  rug_obj <- summarize_Sigmas(output_dir = output_dir)
  discard <- plot_rug(rug = rug_obj$rug,
                      canonical_col_order = ordering$col_order,
                      canonical_row_order = ordering$row_order,
                      row_labels = rug_obj$hosts,
                      save_name = output_dir)
}

# ------------------------------------------------------------------------------
#   Alternative model parameterizations (ASV)
# ------------------------------------------------------------------------------

# We can use the `ordering` variable, as it should still have the canonical
#   ASV-level ordering

alt_parameterizations <- c("asv_days30_diet0_scale1",
                           "asv_days90_diet0_scale1",
                           "asv_days90_diet50_scale1",
                           "asv_days90_diet0_scale2")
for(param in alt_parameterizations) {
  cat("Rendering ALTERNATIVE parameterization", param, "\n")
  rug_obj <- summarize_Sigmas(output_dir = param)
  ordering <- plot_rug(rug_obj$rug,
                       canonical_col_order = ordering$col_order,
                       canonical_row_order = ordering$row_order,
                       row_labels = rug_obj$hosts,
                       save_name = output_dir)
}

# ------------------------------------------------------------------------------
#   Host/row reorderings
# ------------------------------------------------------------------------------

cat("Rendering on PEDIGREE for", level, "\n")
order_obj <- order_rug_row_pedigree(rug_obj)
row_order <- order_obj$order
discard <- plot_rug(rug = rug_obj$rug,
                    canonical_col_order = ordering$col_order,
                    canonical_row_order = row_order,
                    row_labels = rug_obj$hosts,
                    cluster_obj = order_obj$hc,
                    save_name = paste0(output_dir, "_roworder-pedigree"))

cat("Rendering on PRIMARY SOCIAL GROUP for", level, "\n")
order_obj <- order_rug_row_group(rug_obj = rug_obj)
discard <- plot_rug(rug = rug_obj$rug,
                    canonical_col_order = ordering$col_order,
                    canonical_row_order = order_obj$order,
                    row_labels = order_obj$label,
                    save_name = paste0(output_dir, "_roworder-group"))

cat("Rendering on ALPHA-DIVERSITY for", level, "\n")
order_obj <- order_rug_row_diversity(rug_obj = rug_obj)
row_order <- order_obj$order
discard <- plot_rug(rug = rug_obj$rug,
                    canonical_col_order = ordering$col_order,
                    canonical_row_order = row_order,
                    row_labels = rug_obj$hosts[row_order],
                    cluster_obj = order_obj$hc,
                    save_name = paste0(output_dir, "_roworder-diversity"))

# ------------------------------------------------------------------------------
#   Pair/column reorderings
# ------------------------------------------------------------------------------

cat("Rendering on TAXON PAIR PHYLOGENETIC DISTANCE for", level, "\n")
order_obj <- order_rug_col_phylogenetic(rug_obj, data$taxonomy)
ordering <- plot_rug(rug = rug_obj$rug,
                     canonical_col_order = order_obj$order,
                     canonical_row_order = ordering$row_order,
                     row_labels = rug_obj$hosts,
                     save_name = paste0(output_dir, "_colorder-taxonomic"))
