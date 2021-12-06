source("path_fix.R")

library(tidyverse)
library(rulesoflife)
library(RColorBrewer)
library(cowplot)
library(grid)
library(ggridges)
library(ggraph)
library(igraph)

# ------------------------------------------------------------------------------
#
#   Figure 4 - "hockey stick" plot of universality scores, histograms of
#              universality scores, and barplot enrichment visualizations of
#              family-family pairs in the top 2.5% of universal pairs
#
# ------------------------------------------------------------------------------

# These are hard-coded. See `figures_supplemental.R` for the calculation of
# "important" universality scores from permuted data.
thresholds_scores <- data.frame(type = factor(c("Phylum", "Family", "ASV")),
                                x0 = c(0.162, 0.140, 0.136))

# "Important"/significant correlations.
thresholds <- data.frame(type = factor(c("Phylum", "Family", "ASV")),
                         lower = c(-0.303, -0.256, -0.263),
                         upper = c(0.149, 0.207, 0.254))

rug_phy <- summarize_Sigmas(output_dir = "phy_days90_diet25_scale1")
rug_fam <- summarize_Sigmas(output_dir = "fam_days90_diet25_scale1")
rug_asv <- summarize_Sigmas(output_dir = "asv_days90_diet25_scale1")

# ------------------------------------------------------------------------------
#   "Hockey stick" plot
# ------------------------------------------------------------------------------

asv_scores_pieces <- apply(rug_asv$rug, 2, function(x) calc_universality_score(x, return_pieces = TRUE))
score_df <- data.frame(x = asv_scores_pieces[1,], y = asv_scores_pieces[2,])

# Additional labelings
# 1) Consensus sign
score_df$sign <- apply(rug_asv$rug, 2, calc_consensus_sign)
score_df$sign <- factor(score_df$sign, levels = c(1, -1))
levels(score_df$sign) <- c("positive", "negative")
# 2) Percent significant observations for this taxon pair
score_df$signif <- apply(rug_asv$rug, 2, function(x) {
  sum(x < thresholds %>% filter(type == "ASV") %>% pull(lower) | x > thresholds %>% filter(type == "ASV") %>% pull(upper))/length(x)
})

p1 <- ggplot(score_df %>% filter(sign != 0)) +
  geom_point(mapping = aes(x = x, y = y, fill = signif, color = sign),
             size = 2,
             shape = 21,
             stroke = 1) +
  scale_fill_distiller(palette = "PuRd", direction = 1,
                       guide = guide_colorbar(frame.colour = "black",
                                              ticks.colour = "black")) +
  scale_color_manual(values = c("#000000", "#888888")) +
  theme_bw() +
  ylim(c(0,1)) +
  theme(legend.title = element_text(margin = margin(b = 5))) +
  labs(x = "proportion shared sign",
       y = "median association strength",
       fill = "Proportion\nsignificant\nobservations",
       color = "Consensus\ncorrelation sign")

# ------------------------------------------------------------------------------
#   ggridges stacked histograms of universality scores by taxonomic level
# ------------------------------------------------------------------------------

plot_df <- data.frame(score = apply(rug_phy$rug, 2, calc_universality_score),
                      type = "Phylum")
plot_df <- rbind(plot_df,
                 data.frame(score = apply(rug_fam$rug, 2, calc_universality_score),
                            type = "Family"))
plot_df <- rbind(plot_df,
                 data.frame(score = apply(rug_asv$rug, 2, calc_universality_score),
                            type = "ASV"))

p2 <- ggplot(plot_df, aes(x = score, y = type, fill = type)) +
  geom_density_ridges(stat = "binline", bins = 20, scale = 0.9, draw_baseline = TRUE) +
  geom_segment(data = thresholds_scores, aes(x = x0, xend = x0, y = as.numeric(type),
                                             yend = as.numeric(type) + .9),
               color = "black",
               size = 1,
               linetype = "dashed") +
  theme_bw() +
  labs(x = "universality score",
       y = "") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme_ridges(grid = FALSE, center_axis_labels = TRUE) +
  theme(legend.position = "none") +
  xlim(c(-0.05, 0.65)) +
  scale_fill_brewer(palette = "Blues")

p <- plot_grid(p1, p2, ncol = 2,
               rel_widths = c(1.5, 1),
               labels = c("a", "b"),
               label_size = 18,
               # label_x = -0.02,
               scale = 0.95)

ggsave("output/figures/F4ab.svg",
       p,
       dpi = 100,
       units = "in",
       height = 4,
       width = 9)

# ------------------------------------------------------------------------------
#   Enrichment bar plots
# ------------------------------------------------------------------------------

data <- load_data(tax_level = "ASV")

# Note: I'm repeatedly calculating these. It would be much better to do this once
# globally.
scores <- apply(rug_asv$rug, 2, calc_universality_score)
consensus_signs <- apply(rug_asv$rug, 2, calc_consensus_sign)

# Pull top k percent most universal associations
percent <- 2.5
k <- round(length(scores)*(percent/100))
top_pairs <- order(scores, decreasing = TRUE)[1:k]

# Get the taxon indices of each partner in these top pairs
pair_idx1 <- rug_asv$tax_idx1[top_pairs]
pair_idx2 <- rug_asv$tax_idx2[top_pairs]

# ------------------------------------------------------------------------------
#   Family version
# ------------------------------------------------------------------------------

families_top <- c(data$taxonomy[pair_idx1,6], data$taxonomy[pair_idx2,6])
families_top <- families_top[!is.na(families_top)]

families_all <- c(data$taxonomy[rug_asv$tax_idx1,6], data$taxonomy[rug_asv$tax_idx2,6])
families_all <- families_all[!is.na(families_all)]

frequencies_subset <- table(families_top)
frequencies <- table(families_all)

# Test for enrichment statistically
signif <- c()
for(fam in names(frequencies_subset)) {
  fam_in_sample <- unname(unlist(frequencies_subset[fam]))
  sample_size <- unname(unlist(sum(frequencies_subset)))
  fam_in_bg <- unname(unlist(frequencies[fam]))
  bg_size <- unname(unlist(sum(frequencies)))
  ctab <- matrix(c(fam_in_sample,
                   sample_size - fam_in_sample,
                   fam_in_bg,
                   bg_size - fam_in_bg),
                 2, 2, byrow = TRUE)
  prob <- fisher.test(ctab, alternative = "greater")$p.value
  if(prob < 0.05) {
    signif <- c(signif, fam)
    cat(paste0("ASV family: ", fam, ", p-value: ", round(prob, 3), "\n"))
  }
}

plot_enrichment(frequencies_subset1 = frequencies_subset,
                frequencies_subset2 = NULL,
                frequencies = frequencies,
                significant_families1 = signif,
                significant_families2 = NULL,
                plot_height = 6,
                plot_width = 5,
                legend_topmargin = 100,
                use_pairs = FALSE,
                rel_widths = c(1, 0.35, 1, 0.3, 2),
                labels = c("overall", "2.5% most universal ASVs"),
                save_name = "F4-enrichment.png")

# ------------------------------------------------------------------------------
#   Family-family version
# ------------------------------------------------------------------------------

fam1 <- data$taxonomy[pair_idx1,6]
fam2 <- data$taxonomy[pair_idx2,6]
retain_idx <- !is.na(fam1) & !is.na(fam2)
fam1 <- fam1[retain_idx]
fam2 <- fam2[retain_idx]
family_pairs_top <- paste0(fam1, " - ", fam2)

fam1 <- data$taxonomy[rug_asv$tax_idx1,6]
fam2 <- data$taxonomy[rug_asv$tax_idx2,6]
retain_idx <- !is.na(fam1) & !is.na(fam2)
fam1 <- fam1[retain_idx]
fam2 <- fam2[retain_idx]

family_pairs_all <- paste0(fam1, " - ", fam2)

frequencies_subset <- table(family_pairs_top)
frequencies <- table(family_pairs_all)

# Test for enrichment statistically
signif <- c()
for(fam in names(frequencies_subset)) {
  fam_in_sample <- unname(unlist(frequencies_subset[fam]))
  sample_size <- unname(unlist(sum(frequencies_subset)))
  fam_in_bg <- unname(unlist(frequencies[fam]))
  bg_size <- unname(unlist(sum(frequencies)))
  ctab <- matrix(c(fam_in_sample,
                   sample_size - fam_in_sample,
                   fam_in_bg,
                   bg_size - fam_in_bg),
                 2, 2, byrow = TRUE)
  prob <- fisher.test(ctab, alternative = "greater")$p.value
  if(prob < 0.05) {
    signif <- c(signif, fam)
    cat(paste0("ASV family: ", fam, ", p-value: ", round(prob, 3), "\n"))
  }
}

plot_enrichment(frequencies_subset1 = frequencies_subset,
                frequencies_subset2 = NULL,
                frequencies = frequencies,
                significant_families1 = signif,
                significant_families2 = NULL,
                plot_height = 6,
                plot_width = 6,
                legend_topmargin = 100,
                use_pairs = TRUE,
                rel_widths = c(1, 0.35, 1, 0.3, 2.75),
                labels = c("overall", "2.5% most universal ASVs"),
                save_name = "F4-enrichment-pair.png")

# ------------------------------------------------------------------------------
#   Network plot
# ------------------------------------------------------------------------------

# Build a mapping of taxon indices to new labels 1..n
map_df <- data.frame(taxon_idx = unique(c(pair_idx1, pair_idx2)))
map_df$name <- 1:nrow(map_df)

# Build a node data.frame with columns name (index 1..n) and family
node_df <- map_df
node_df$Family <- sapply(node_df$taxon_idx, function(x) {
  highest_tax_level <- max(which(!is.na(data$taxonomy[x,])))
  if(highest_tax_level >= 6) {
    data$taxonomy[x,6]
  } else {
    "Unknown"
  }
})
node_df <- node_df[,c(2:3,1)]

# Build an edge data.frame with columns from, to, sign, and score
edge1 <- data.frame(taxon_idx = pair_idx1)
edge1 <- left_join(edge1, map_df, by = "taxon_idx")$name
edge2 <- data.frame(taxon_idx = pair_idx2)
edge2 <- left_join(edge2, map_df, by = "taxon_idx")$name

edge_df <- data.frame(from = edge1,
                      to = edge2,
                      Sign = factor(consensus_signs[top_pairs], levels = c("1", "-1")),
                      score = scores[top_pairs])
levels(edge_df$Sign) <- c("positive", "negative")

fam_df <- edge_df %>%
  left_join(node_df, by = c("from" = "name")) %>%
  select(to, Family1 = Family)
fam_df <- fam_df %>%
  left_join(node_df, by = c("to" = "name")) %>%
  select(Family1, Family2 = Family)

graph <- graph_from_data_frame(edge_df, node_df, directed = FALSE)

# Not specifying the layout - defaults to "auto"
# fr and kk layouts are ok here
p <- ggraph(graph, layout = "fr") +
  geom_edge_link(aes(color = Sign), width = 2, alpha = 1) +
  geom_node_point(aes(color = Family), size = 5) +
  geom_node_label(aes(label = taxon_idx), size = 3, repel = TRUE) +
  scale_colour_manual(values = family_palette[order(names(family_palette))]) +
  scale_edge_colour_manual(values = c(negative = "gray", positive = "black")) +
  labs(x = "dimension 1",
       y = "dimension 2") +
  theme_bw() +
  guides(edge_color = guide_legend(title = "Consensus\ncorrelation sign"))

ggsave(file.path("output", "figures", "F5.svg"),
       p,
       units = "in",
       dpi = 100,
       height = 7,
       width = 9)
