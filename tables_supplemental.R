source("path_fix.R")

library(rulesoflife)
library(tidyverse)
library(ape)

# source("ggplot_fix.R")

# ------------------------------------------------------------------------------
#   Print table of all phyla / families / ASVs
# ------------------------------------------------------------------------------

tax_phy <- load_data(tax_level = "phylum")$taxonomy
tax_phy <- cbind(id = 1:nrow(tax_phy), tax_phy[,c(2:3,1)])
write.table(tax_phy,
            file.path("output", "phy_table.tsv"),
            sep = "\t",
            quote = F,
            row.names = F)

tax_fam <- load_data(tax_level = "family")$taxonomy
tax_fam <- cbind(id = 1:nrow(tax_fam), tax_fam[,c(2:6,1)])
write.table(tax_fam,
            file.path("output", "fam_table.tsv"),
            sep = "\t",
            quote = F,
            row.names = F)

tax_asv <- load_data(tax_level = "ASV")$taxonomy
tax_asv <- cbind(id = 1:nrow(tax_asv), tax_asv[,c(2:7,1)])
write.table(tax_asv,
            file.path("output", "asv_table.tsv"),
            sep = "\t",
            quote = F,
            row.names = F)

# ------------------------------------------------------------------------------
#   Pre-calculate distances across sequences
# ------------------------------------------------------------------------------

sequence_identity <- sequence_distance(distance_type = "raw")
mismatches <- sequence_distance(distance_type = "N")
K80 <- sequence_distance(distance_type = "K80")

# ------------------------------------------------------------------------------
#   Print table of top universal pairs
# ------------------------------------------------------------------------------

data <- load_data(tax_level = "ASV")
output_dir <- "asv_days90_diet25_scale1"
rug_obj <- summarize_Sigmas(output_dir)

scores <- apply(rug_obj$rug, 2, calc_universality_score)
piecewise_scores <- apply(rug_obj$rug, 2, function(x) {
  calc_universality_score(x, return_pieces = TRUE)
})
consensus_signs <- apply(rug_obj$rug, 2, calc_consensus_sign)

scores_df <- data.frame(index = 1:length(scores),
                        score = scores,
                        prop_agree = piecewise_scores[1,],
                        med_assoc_strength = piecewise_scores[2,],
                        consensus_sign = consensus_signs)

percents <- c(100, 5, 2.5, 1)
for(percent in percents) {
  k <- percent_to_k(percent, nrow(scores_df))
  if(percent == 100) {
    k <- nrow(scores_df)
  }
  table_df <- scores_df %>%
    arrange(desc(score)) %>%
    slice(1:k) %>%
    arrange(desc(score))
  table_df$rank <- 1:k
  table_df$tax_index1 <- sapply(table_df$index, function(x) {
    rug_obj$tax_idx1[x]
  })
  table_df$tax_index2 <- sapply(table_df$index, function(x) {
    rug_obj$tax_idx2[x]
  })
  table_df$taxonomy_1 <- sapply(table_df$tax_index1, function(x) {
    tax_levels <- colnames(data$taxonomy)[2:ncol(data$taxonomy)]
    paste(paste(tax_levels,
                data$taxonomy[x,2:ncol(data$taxonomy)]), collapse = " / ")
  })
  table_df$sequence_1 <- sapply(table_df$tax_index1, function(x) {
    data$taxonomy[x,1]
  })
  table_df$taxonomy_2 <- sapply(table_df$tax_index2, function(x) {
    tax_levels <- colnames(data$taxonomy)[2:ncol(data$taxonomy)]
    paste(paste(tax_levels,
                data$taxonomy[x,2:ncol(data$taxonomy)]), collapse = " / ")
  })
  table_df$sequence_2 <- sapply(table_df$tax_index2, function(x) {
    data$taxonomy[x,1]
  })

  # Add phylogenetic distances from ape::dist.dna()
  table_df$percent_identity <- NA
  table_df$mismatches <- NA
  table_df$K80_dist <- NA
  for(i in 1:nrow(table_df)) {
    table_df$percent_identity[i] <- sequence_identity[table_df$tax_index1[i],
                                                      table_df$tax_index2[i]]
    table_df$mismatches[i] <- mismatches[table_df$tax_index1[i],
                                         table_df$tax_index2[i]]
    table_df$K80_dist[i] <- K80[table_df$tax_index1[i],
                                table_df$tax_index2[i]]
  }

  write.table(table_df %>% dplyr::select(rank, tax_index1, tax_index2,
                                         consensus_sign, prop_agree, med_assoc_strength, score,
                                         percent_identity, mismatches, K80_dist,
                                         taxonomy_1, taxonomy_2, sequence_1, sequence_2),
              file = file.path("output", paste0("top_", percent, "_universal.tsv")),
              sep = "\t",
              quote = FALSE,
              row.names = FALSE)

  # print(k)
  # print(table(table_df$consensus_sign) / k)
}
