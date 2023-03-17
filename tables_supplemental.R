source("path_fix.R")

library(rulesoflife)
library(tidyverse)
library(ape)

# source("ggplot_fix.R")

# ------------------------------------------------------------------------------
#   Print table of all phyla / families / ASVs
# ------------------------------------------------------------------------------

# Renumber (TBD)
# representation <- represented_taxa(filtered_pairs)
# scores_df$tax_index1_renumbered <- sapply(scores_df$tax_index1, function(x) renumber_taxon(representation, x))
# scores_df$tax_index2_renumbered <- sapply(scores_df$tax_index2, function(x) renumber_taxon(representation, x))

data_phy <- load_data(tax_level = "phylum")
filtered_pairs <- filter_joint_zeros(data_phy$counts, threshold_and = 0.05, threshold_or = 0.5)
representation <- represented_taxa(filtered_pairs)
tax_phy <- data_phy$taxonomy
tax_phy <- tax_phy[1:(nrow(tax_phy)-1),]
tax_phy <- cbind(id_raw = 1:nrow(tax_phy), tax_phy[,c(2:3,1)])
tax_phy$id <- sapply(tax_phy$id_raw, function(x) renumber_taxon(representation, x))
tax_phy <- tax_phy %>%
  filter(!is.na(id))
write.table(tax_phy %>% select(id, domain, phylum),
            file.path("output", "Table_S2.tsv"),
            sep = "\t",
            quote = F,
            row.names = F)

data_fam <- load_data(tax_level = "family")
filtered_pairs <- filter_joint_zeros(data_fam$counts, threshold_and = 0.05, threshold_or = 0.5)
representation <- represented_taxa(filtered_pairs)
tax_fam <- data_fam$taxonomy
tax_fam <- tax_fam[1:(nrow(tax_fam)-1),]
tax_fam <- cbind(id_raw = 1:nrow(tax_fam), tax_fam[,c(2:6,1)])
tax_fam$id <- sapply(tax_fam$id_raw, function(x) renumber_taxon(representation, x))
tax_fam <- tax_fam %>%
  filter(!is.na(id))
write.table(tax_fam %>% select(id, domain, phylum, class, order, family),
            file.path("output", "Table_S3.tsv"),
            sep = "\t",
            quote = F,
            row.names = F)


tax_asv <- load_data(tax_level = "ASV")$taxonomy
tax_asv <- cbind(id = 1:nrow(tax_asv), tax_asv[,c(2:7,1)])
write.table(tax_asv,
            file.path("output", "Table_S1.tsv"),
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

filtered_pairs <- filter_joint_zeros(data$counts, threshold_and = 0.05, threshold_or = 0.5)

scores <- apply(rug_obj$rug, 2, calc_universality_score)
piecewise_scores <- apply(rug_obj$rug, 2, function(x) {
  calc_universality_score(x, return_pieces = TRUE)
})
consensus_signs <- apply(rug_obj$rug, 2, calc_consensus_sign)

scores_df <- data.frame(index = 1:length(scores),
                        score = scores,
                        threshold = filtered_pairs$threshold,
                        prop_agree = piecewise_scores[1,],
                        med_assoc_strength = piecewise_scores[2,],
                        consensus_sign = consensus_signs)

scores_df %<>%
  filter(threshold) %>%
  arrange(desc(score))

k5 <- percent_to_k(5, nrow(scores_df))
k10 <- percent_to_k(10, nrow(scores_df))
scores_df$percent <- "top 100%"
scores_df$percent[1:k10] <- "top 10%"
scores_df$percent[1:k5] <- "top 5%"

scores_df$rank <- 1:nrow(scores_df)
scores_df$percent_identity <- NA
scores_df$mismatches <- NA
scores_df$K80_dist <- NA

scores_df %<>%
  filter(percent %in% c("top 5%", "top 10%"))

for(i in 1:nrow(scores_df)) {
  scores_df$tax_index1 <- sapply(scores_df$index, function(x) {
    rug_obj$tax_idx1[x]
  })
  scores_df$tax_index2 <- sapply(scores_df$index, function(x) {
    rug_obj$tax_idx2[x]
  })
  scores_df$taxonomy_1 <- sapply(scores_df$tax_index1, function(x) {
    tax_levels <- colnames(data$taxonomy)[2:ncol(data$taxonomy)]
    paste(paste(tax_levels,
                data$taxonomy[x,2:ncol(data$taxonomy)]), collapse = " / ")
  })
  scores_df$sequence_1 <- sapply(scores_df$tax_index1, function(x) {
    data$taxonomy[x,1]
  })
  scores_df$taxonomy_2 <- sapply(scores_df$tax_index2, function(x) {
    tax_levels <- colnames(data$taxonomy)[2:ncol(data$taxonomy)]
    paste(paste(tax_levels,
                data$taxonomy[x,2:ncol(data$taxonomy)]), collapse = " / ")
  })
  scores_df$sequence_2 <- sapply(scores_df$tax_index2, function(x) {
    data$taxonomy[x,1]
  })

  # Add phylogenetic distances from ape::dist.dna()
  scores_df$percent_identity[i] <- sequence_identity[scores_df$tax_index1[i],
                                                     scores_df$tax_index2[i]]
  scores_df$mismatches[i] <- mismatches[scores_df$tax_index1[i],
                                        scores_df$tax_index2[i]]
  scores_df$K80_dist[i] <- K80[scores_df$tax_index1[i],
                               scores_df$tax_index2[i]]
}

# Renumber
representation <- represented_taxa(filtered_pairs)
scores_df$tax_index1_renumbered <- sapply(scores_df$tax_index1, function(x) renumber_taxon(representation, x))
scores_df$tax_index2_renumbered <- sapply(scores_df$tax_index2, function(x) renumber_taxon(representation, x))

write.table(scores_df %>%
              dplyr::select(rank, percent, tax_index1_renumbered, tax_index2_renumbered,
                            consensus_sign, prop_agree, med_assoc_strength,
                            percent_identity, mismatches, K80_dist,
                            taxonomy_1, taxonomy_2, sequence_1, sequence_2),
            file = file.path("output", "Table_S4.tsv"),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

# Is there an enrichment of sign in the top 5% of pairs by universality score?
binom.test(table(scores_df %>% filter(percent == "top 5%") %>% pull(consensus_sign)), k5, 0.5)
