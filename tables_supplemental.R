source("path_fix.R")

library(rulesoflife)
library(tidyverse)

source("ggplot_fix.R")

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
#   Print table of top 250 universal pairs
# ------------------------------------------------------------------------------

output_dir <- "asv_days90_diet25_scale1"
rug_obj <- summarize_Sigmas(output_dir)

scores <- apply(rug_obj$rug, 2, calc_universality_score)
consensus_signs <- apply(rug_obj$rug, 2, calc_consensus_sign)

scores_df <- data.frame(index = 1:length(scores),
                        score = scores,
                        sign = consensus_signs)

percents <- c(5, 2.5, 1)
for(percent in percents) {
  k <- percent_to_k(percent, nrow(scores_df))
  table_df <- scores_df %>%
    arrange(desc(score)) %>%
    slice(1:k) %>%
    select(index, score) %>%
    arrange(desc(score))
  table_df$rank <- 1:k
  table_df$consensus_sign <- sapply(table_df$index, function(x) {
    consensus_signs[x]
  })
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

  write.table(table_df %>% select(rank, score, tax_index1, tax_index2,
                                  consensus_sign, taxonomy_1, taxonomy_2,
                                  sequence_1, sequence_2),
              file = file.path("output", paste0("top_", percent, "_universal.tsv")),
              sep = "\t",
              quote = FALSE,
              row.names = FALSE)

  print(k)
  print(table(table_df$consensus_sign) / k)
}
