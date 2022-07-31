source("path_fix.R")

library(tidyverse)
library(rulesoflife)
library(RColorBrewer)
library(cowplot)
library(grid)
library(ggridges)
library(ggraph)
library(igraph)
library(scales)

# source("thresholds.R")

rug_obj <- summarize_Sigmas(output_dir = "asv_days90_diet25_scale1")
rug <- rug_obj$rug

# ------------------------------------------------------------------------------
#   "Hockey stick" plot
# ------------------------------------------------------------------------------

scores_pieces <- apply(rug, 2, function(x) calc_universality_score(x, return_pieces = TRUE))
scores <- apply(scores_pieces, 2, function(x) x[1]*x[2])
score_df <- data.frame(x = scores_pieces[1,], y = scores_pieces[2,])

# Check number of pairs with 100% sign agreement
sign_agree <- which(scores_pieces[1,] == 1)
sign_disagree <- which(scores_pieces[1,] == 0.5)
cat(paste0(length(sign_agree), " pairs have 100% sign agreement\n"))
cat(paste0(round(median(abs(rug[,sign_agree])), 3), " median assocation strength (100% agree)\n"))
cat(paste0(round(median(abs(rug[,sign_disagree])), 3), " median assocation strength (<100% agree)\n"))

# Additional labelings
# 1) Consensus sign
score_df$sign <- apply(rug, 2, calc_consensus_sign)
score_df$sign <- factor(score_df$sign, levels = c(1, -1))
levels(score_df$sign) <- c("positive", "negative")
# 2) Percent significant observations for this taxon pair
score_df$signif <- apply(rug, 2, function(x) {
  sum(x < thresholds %>% filter(type == "ASV") %>% pull(lower) | x > thresholds %>% filter(type == "ASV") %>% pull(upper))/length(x)
})

p1a <- ggplot(score_df %>% filter(sign == "negative") %>% mutate(overall = "Consensus negatively associated pairs")) +
  # geom_point(mapping = aes(x = x, y = y, fill = y),
  geom_point(mapping = aes(x = x, y = y),
             size = 2,
             shape = 21,
             fill = muted("navy")) +
  theme_bw() +
  ylim(c(0,1)) +
  facet_wrap(~ overall) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 11),
        axis.title.x = element_text(size = 13),
        axis.text.y = element_text(size = 11),
        axis.title.y = element_text(size = 13),
        strip.text.x = element_text(size = 12)) +
  labs(x = "proportion shared sign",
       y = "median correlation strength",
       color = "Consensus\ncorrelation sign")

p1b <- ggplot(score_df %>% filter(sign == "positive") %>% mutate(overall = "Consensus positively associated pairs")) +
  geom_point(mapping = aes(x = x, y = y),
             size = 2,
             shape = 21,
             fill = "red") +
  theme_bw() +
  ylim(c(0,1)) +
  facet_wrap(~ overall) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 11),
        axis.title.x = element_text(size = 13),
        axis.text.y = element_text(size = 11),
        axis.title.y = element_text(size = 13),
        strip.text.x = element_text(size = 12)) +
  labs(x = "proportion shared sign",
       y = "median correlation strength",
       color = "Consensus\ncorrelation sign")

# ------------------------------------------------------------------------------
#   ggridges stacked histograms of universality scores at the ASV level
# ------------------------------------------------------------------------------

plot_df <- data.frame(score = scores, type = "ASV")

p2 <- ggplot(data.frame(score = apply(rug, 2, calc_universality_score)),
             aes(x = score, y = 1)) +
  geom_density_ridges(stat = "binline", bins = 20, scale = 0.9, draw_baseline = TRUE) +
  theme_bw() +
  labs(x = "universality score",
       y = "") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme_ridges(grid = FALSE, center_axis_labels = TRUE) +
  theme(legend.position = "none") +
  xlim(c(-0.05, 0.65)) +
  scale_fill_brewer(palette = "Blues") +
  theme(axis.text = element_text(size = 9),
        axis.title = element_text(size = 11))

# ------------------------------------------------------------------------------
#   Enrichment bar plots
# ------------------------------------------------------------------------------

data <- load_data(tax_level = "ASV")

consensus_signs <- apply(rug, 2, calc_consensus_sign)

# Pull top k percent most universal associations
percent <- 2.5
k <- round(length(scores)*(percent/100))
top_pairs <- order(scores, decreasing = TRUE)[1:k]

# Get the taxon indices of each partner in these top pairs
pair_idx1 <- rug_obj$tax_idx1[top_pairs]
pair_idx2 <- rug_obj$tax_idx2[top_pairs]

# ------------------------------------------------------------------------------
#   Family version
# ------------------------------------------------------------------------------

families_top <- c(data$taxonomy[pair_idx1,6], data$taxonomy[pair_idx2,6])
families_top <- families_top[!is.na(families_top)]

families_all <- c(data$taxonomy[rug_obj$tax_idx1,6], data$taxonomy[rug_obj$tax_idx2,6])
families_all <- families_all[!is.na(families_all)]

frequencies_subset <- table(families_top)
frequencies <- table(families_all)

enrichment <- NULL

e_obj <- plot_enrichment(frequencies,
                         frequencies_subset,
                         type_label = "family",
                         location_label = "top ASVs",
                         enrichment = enrichment,
                         cap_size = 5,
                         pt_sz = 3,
                         title_text_sz = 14,
                         text_sz = 12,
                         stroke_sz = 1.25)
p4 <- e_obj$p
enrichment <- e_obj$enrichment

# ------------------------------------------------------------------------------
#   Family-family version
# ------------------------------------------------------------------------------

fam1 <- data$taxonomy[pair_idx1,6]
fam2 <- data$taxonomy[pair_idx2,6]
retain_idx <- !is.na(fam1) & !is.na(fam2)
fam1 <- fam1[retain_idx]
fam2 <- fam2[retain_idx]
family_pairs_top <- paste0(fam1, " - ", fam2)

fam1 <- data$taxonomy[rug_obj$tax_idx1,6]
fam2 <- data$taxonomy[rug_obj$tax_idx2,6]
retain_idx <- !is.na(fam1) & !is.na(fam2)
fam1 <- fam1[retain_idx]
fam2 <- fam2[retain_idx]

family_pairs_all <- paste0(fam1, " - ", fam2)

frequencies_subset <- table(family_pairs_top)
frequencies <- table(family_pairs_all)

e_obj <- plot_enrichment(frequencies,
                         frequencies_subset,
                         type_label = "family-pair",
                         location_label = "top ASVs",
                         enrichment = enrichment,
                         cap_size = 5,
                         pt_sz = 3,
                         title_text_sz = 14,
                         text_sz = 12,
                         stroke_sz = 1.25)
p5 <- e_obj$p
enrichment <- e_obj$enrichment

enrichment %<>%
  arrange(type, name) %>%
  dplyr::select(name, type, location, oddsratio, lower95, upper95, pvalue, qvalue) %>%
  rename(`ASV family or pair name` = name,
         `Type` = type,
         `Enrichment evaluated in` = location,
         `Odds ratio` = oddsratio,
         `Lower 95% CI` = lower95,
         `Upper 95% CI` = upper95,
         `P-value (Fisher's exact test)` = pvalue,
         `Adj. p-value (Benjamini-Hochberg)` = qvalue)
write.table(enrichment,
            file = file.path("output", "enrichment_top-pairs.tsv"),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

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
  dplyr::select(to, Family1 = Family)
fam_df <- fam_df %>%
  left_join(node_df, by = c("to" = "name")) %>%
  dplyr::select(Family1, Family2 = Family)

graph <- graph_from_data_frame(edge_df, node_df, directed = FALSE)

family_palette <- readRDS(file.path("output", "family_palette.rds"))

# Not specifying the layout - defaults to "auto"
# fr and kk layouts are ok here
p3 <- ggraph(graph, layout = "fr") +
  geom_edge_link(aes(color = Sign), width = 1.5, alpha = 1) +
  geom_node_point(aes(color = Family), size = 5) +
  geom_node_label(aes(label = taxon_idx), size = 4, label.size = NA, repel = TRUE, label.padding = unit(0.08, "lines")) +
  scale_colour_manual(values = family_palette[order(names(family_palette))]) +
  scale_edge_colour_manual(values = c(negative = "gray", positive = "black")) +
  labs(x = "dimension 1",
       y = "dimension 2") +
  theme_bw() +
  guides(edge_color = "none") +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank())

col1 <- plot_grid(p1b, p1a, p2,
                  ncol = 1,
                  rel_heights = c(1, 1, 0.8),
                  labels = c("A", "B", "C"),
                  label_size = 18,
                  label_y = 1.02,
                  scale = 1)
row1 <- plot_grid(p4, p5,
                  ncol = 2,
                  rel_widths = c(1, 1),
                  labels = c("E", "F"),
                  label_size = 18,
                  label_y = 1.02,
                  label_x = -0.02,
                  scale = 0.95)
p4_padded <- plot_grid(p3,
                       labels = c("D"),
                       label_size = 18,
                       label_y = 1.01,
                       label_x = -0.02,
                       scale = 0.95)
col2 <- plot_grid(p4_padded, row1,
                  ncol = 1,
                  rel_heights = c(1, 0.75))
p <- plot_grid(col1, NULL, col2,
               ncol = 3,
               rel_widths = c(1, 0.05, 2))

ggsave(file.path("output", "figures", "network.svg"),
       p,
       units = "in",
       dpi = 100,
       height = 11,
       width = 12)

