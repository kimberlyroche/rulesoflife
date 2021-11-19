source("path_fix.R")

library(tidyverse)
library(rulesoflife)
library(driver)
library(cowplot)
library(fido)

# ------------------------------------------------------------------------------
#
#   Supplemental Figure S6
#
# ------------------------------------------------------------------------------

data <- load_data(tax_level = "ASV")
rug_asv <- summarize_Sigmas(output_dir = "asv_days90_diet25_scale1")

phy_dist <- rug_phylogenetic_distances(rug_asv, data$taxonomy, as_matrix = FALSE)
scores <- apply(rug_asv$rug, 2, calc_universality_score)
signs <- apply(rug_asv$rug, 2, calc_consensus_sign)
median_assoc_strength <- apply(rug_asv$rug, 2, function(x) {
  median(abs(x))
})

plot_df <- data.frame(d = phy_dist,
                      score = scores,
                      sign = signs,
                      mas = median_assoc_strength)
plot_df$sign <- factor(plot_df$sign, levels = c(1, 0, -1))
levels(plot_df$sign) <- c("positive", NA, "negative")

phylo_neg_mean <- mean(phy_dist[signs < 0])
phylo_pos_mean <- mean(phy_dist[signs > 0])

p <- ggplot(plot_df %>% filter(!is.na(sign)), aes(x = d, y = score, fill = factor(sign))) +
  geom_segment(aes(x = phylo_neg_mean, xend = phylo_neg_mean, y = 0.26, yend = 0.8),
               color = "#34CCDE",
               size = 2,
               linetype = "dashed") +
  geom_segment(aes(x = phylo_pos_mean, xend = phylo_pos_mean, y = 0.26, yend = 0.8),
               color = "#F25250",
               size = 2,
               linetype = "dashed") +
  geom_point(size = 2, shape = 21) +
  scale_fill_manual(values = c("#F25250", "#34CCDE")) +
  theme_bw() +
  labs(x = "phylogenetic distance",
       y = "universality score",
       fill = "Consensus\ncorrelation sign")

ggsave("output/figures/SF6a.svg",
       p,
       dpi = 100,
       units = "in",
       height = 5,
       width = 6.5)

p <- ggplot(plot_df %>% filter(!is.na(sign)), aes(x = d, y = mas, fill = factor(sign))) +
  geom_segment(aes(x = phylo_neg_mean, xend = phylo_neg_mean, y = 0.29, yend = 0.8),
               color = "#34CCDE",
               size = 2,
               linetype = "dashed") +
  geom_segment(aes(x = phylo_pos_mean, xend = phylo_pos_mean, y = 0.29, yend = 0.8),
               color = "#F25250",
               size = 2,
               linetype = "dashed") +
  geom_point(size = 2, shape = 21) +
  scale_fill_manual(values = c("#F25250", "#34CCDE")) +
  theme_bw() +
  labs(x = "phylogenetic distance",
       y = "median association strength",
       fill = "Consensus\ncorrelation sign")

ggsave("output/figures/SF6b.svg",
       p,
       dpi = 100,
       units = "in",
       height = 5,
       width = 6.5)

# ------------------------------------------------------------------------------
#   Enrichment of top left-hand pairs
# ------------------------------------------------------------------------------

topleft_pairs <- which(phy_dist < 0.2 & median_assoc_strength > 0.5)

all_pairs <- data.frame(idx1 = rug_asv$tax_idx1,
                        idx2 = rug_asv$tax_idx2,
                        topleft = FALSE)
all_pairs$tax1 <- sapply(1:nrow(all_pairs), function(x) {
  paste0(data$taxonomy[all_pairs$idx1[x],6], collapse = "/")
})
all_pairs$tax2 <- sapply(1:nrow(all_pairs), function(x) {
  paste0(data$taxonomy[all_pairs$idx2[x],6], collapse = "/")
})

for(i in 1:length(topleft_pairs)) {
  all_pairs$topleft[all_pairs$idx1 == rug_asv$tax_idx1[topleft_pairs[i]] &
                      all_pairs$idx2 == rug_asv$tax_idx2[topleft_pairs[i]]] <- TRUE
}

# 6105 of 8911 family-family pairs don't involve "NA" families
all_pairs_noNA <- all_pairs %>%
  filter(tax1 != "NA" & tax2 != "NA")

all_pairs_noNA <- all_pairs_noNA %>%
  mutate(taxpair = paste0(tax1, " - ", tax2))

# Calculate a baseline frequency for all family-family pairs
frequencies <- table(all_pairs_noNA$taxpair)

# Define a huge color palette over all observed family-family pairs
palette_fn <- "output/family-family_palette.rds"
if(file.exists(palette_fn)) {
  fam_pair_palette <- readRDS(palette_fn)
} else {
  fam_pair_palette <- generate_highcontrast_palette(length(frequencies))
  names(fam_pair_palette) <- names(frequencies)
  saveRDS(fam_pair_palette, palette_fn)
}

# ------------------------------------------------------------------------------
#   Visualize differences in relative abundance of family pairs for top center
#   clump of pairs
# ------------------------------------------------------------------------------

# Stacked bar: plot observed relative representation of family-family pairs
frequencies_subset <- table(all_pairs_noNA$taxpair[all_pairs_noNA$topleft == TRUE])

plot_enrichment(frequencies_subset,
                frequencies,
                plot_height = 6,
                plot_width = 8,
                legend_topmargin = 60,
                legend_leftmargin = -0.5,
                save_name = "SF6c.svg")
