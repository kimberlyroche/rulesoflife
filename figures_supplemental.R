source("path_fix.R")

library(tidyverse)
library(rulesoflife)
library(driver)
library(cowplot)
library(fido)
library(scales)

data <- load_data(tax_level = "ASV")
clr_counts <- clr_array(data$counts + 0.5, parts = 1)

# ------------------------------------------------------------------------------
#
#   Supplemental Figure S2
#
# ------------------------------------------------------------------------------

PCA <- prcomp(t(clr_counts))
PoV <- PCA$sdev^2 / sum(PCA$sdev^2)

rgb_colors <- matrix(c(240, 7, 7,
                       227, 125, 16,
                       68, 194, 10,
                       90, 141, 242,
                       200, 90, 242), 5, 3, byrow = TRUE)/255
hex_colors <- apply(rgb_colors, 1, function(x) {
  rgb(x[1], x[2], x[3])
})

ref_hosts <- c("DUI", "DUX", "LIW", "PEB", "VET")
year_range <- seq(from = 2003, to = 2008, by = 1)

labels_years <- data$metadata %>%
  mutate(year = substr(collection_date, 1, 4)) %>%
  mutate(label = ifelse(sname %in% ref_hosts & year %in% year_range,
                        sname,
                        "other")) %>%
  select(label, year)
labels <- factor(labels_years$label, levels = c(ref_hosts, "other"))
years <- as.numeric(labels_years$year)

plot_df <- data.frame(x = PCA$x[,1],
                      y = PCA$x[,2],
                      label = labels,
                      year = years)
plot_df <- plot_df[plot_df$year %in% year_range,]

p <- ggplot(plot_df) +
  geom_point(data = plot_df[plot_df$label == "other",], mapping = aes(x = x, y = y),
             size = 2,
             shape = 21,
             fill = "#bbbbbb",
             color = "#888888") +
  geom_point(data = plot_df[plot_df$label != "other",], mapping = aes(x = x, y = y, fill = label),
             size = 3,
             shape = 21) +
  scale_fill_manual(values = hex_colors) +
  facet_wrap(. ~ year) +
  labs(title = "Selected host samples (2004-2009)",
       fill = "Host",
       x = paste0("PC1 (", round(PoV[1], 3)*100 ,"% variance expl.)"),
       y = paste0("PC2 (", round(PoV[2], 3)*100 ,"% variance expl.)")) +
  theme_bw()
ggsave("output/figures/SF1.svg",
       plot = p,
       units = "in",
       dpi = 100,
       height = 6,
       width = 10)

# ------------------------------------------------------------------------------
#
#   Supplemental Figure S4
#
# ------------------------------------------------------------------------------

# These are hard-coded. See `figures_supplemental.R` for the calculation of
# "important"/significant correlations from permuted data.
thresholds <- data.frame(type = factor(c("Phylum", "Family", "ASV")),
                         lower = c(-0.303, -0.256, -0.263),
                         upper = c(0.149, 0.207, 0.254))

fam_scores_pieces <- apply(rug_fam$rug, 2, function(x) calc_universality_score(x, return_pieces = TRUE))
fam_score_df <- data.frame(x = fam_scores_pieces[1,], y = fam_scores_pieces[2,])

# Code below copied from `figures_main.R`
fam_score_df$sign <- apply(rug_fam$rug, 2, calc_consensus_sign)
fam_score_df$sign <- factor(fam_score_df$sign, levels = c(1, -1))
levels(fam_score_df$sign) <- c("positive", "negative")
# 2) Percent significant observations for this taxon pair
fam_score_df$signif <- apply(rug_fam$rug, 2, function(x) {
  sum(x < thresholds %>% filter(type == "ASV") %>% pull(lower) | x > thresholds %>% filter(type == "ASV") %>% pull(upper))/length(x)
})

phy_scores_pieces <- apply(rug_phy$rug, 2, function(x) calc_universality_score(x, return_pieces = TRUE))
phy_score_df <- data.frame(x = phy_scores_pieces[1,], y = phy_scores_pieces[2,])

# Code below copied from `figures_main.R`
phy_score_df$sign <- apply(rug_phy$rug, 2, calc_consensus_sign)
phy_score_df$sign <- factor(phy_score_df$sign, levels = c(1, -1))
levels(phy_score_df$sign) <- c("positive", "negative")
# 2) Percent significant observations for this taxon pair
phy_score_df$signif <- apply(rug_phy$rug, 2, function(x) {
  sum(x < thresholds %>% filter(type == "ASV") %>% pull(lower) | x > thresholds %>% filter(type == "ASV") %>% pull(upper))/length(x)
})

score_df <- rbind(cbind(fam_score_df, type = "Family"),
                  cbind(phy_score_df, type = "Phylum"))

p <- ggplot(score_df %>% filter(sign != 0)) +
  geom_point(mapping = aes(x = x, y = y, fill = signif, color = sign),
             size = 2,
             shape = 21,
             stroke = 1) +
  scale_fill_distiller(palette = "PuRd", direction = 1,
                       guide = guide_colorbar(frame.colour = "black",
                                              ticks.colour = "black")) +
  scale_color_manual(values = c("#000000", "#888888")) +
  facet_wrap(. ~ type) +
  theme_bw() +
  ylim(c(0,1)) +
  theme(legend.title = element_text(margin = margin(b = 5))) +
  labs(x = "proportion shared sign",
       y = "median association strength",
       fill = "Proportion\nsignificant\nobservations",
       color = "Consensus\ncorrelation sign")
ggsave("output/figures/SF3.svg",
       p,
       dpi = 100,
       units = "in",
       height = 4,
       width = 8)

# ------------------------------------------------------------------------------
#
#   Supplemental Figure S5
#
# ------------------------------------------------------------------------------

clr_mean_abundance <- rowMeans(clr_counts)

pairs <- read.table(file.path("output", "top_5_universal.tsv"), sep = "\t", header = TRUE)

# Remove sequences themselves
pairs <- pairs %>%
  select(-c(sequence_1, sequence_2))

# Pull taxon IDs and readable names from each pair only
pairs <- rbind(pairs %>% select(a = tax_index1, b = taxonomy_1),
               pairs %>% select(a = tax_index2, b = taxonomy_2))

# Tally the appearances of each taxon
pairs <- pairs %>%
  group_by(a) %>%
  mutate(count = n()) %>%
  distinct() %>%
  arrange(desc(count))

# Generate short names for each taxon
pairs$b <- unname(sapply(pairs$b, function(x) {
  tax_pieces <- str_split(x, " / ")[[1]]
  for(i in length(tax_pieces):1) {
    pieces <- str_split(tax_pieces[i], " ")[[1]]
    if(pieces[2] != "NA") {
      return(paste0("ASV in ", tax_pieces[i]))
    }
  }
  return(NA)
}))

# Attach the mean CLR abundance to the appropriate taxon
pairs <- pairs %>%
  left_join(data.frame(a = 1:length(clr_mean_abundance),
                       clr_mean = clr_mean_abundance),
            by = "a")

# Fix the order in which taxa will plot (by frequency; `x`)
pairs <- pairs %>%
  arrange(count, b)
pairs$x <- 1:nrow(pairs)

# Filter to taxa with at least a few occurrences
pairs <- pairs %>%
  filter(count > 5)

p <- ggplot(pairs, aes(x = factor(x), y = count, fill = clr_mean)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  scale_fill_distiller(palette = "PiYG") +
  scale_x_discrete(breaks = pairs$x,
                   labels = pairs$b) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x = "",
       y = "frequency in top 250 pairs",
       fill = "Mean CLR abundance") +
  coord_flip()
p
ggsave("output/figures/SF3.svg",
       p,
       dpi = 100,
       units = "in",
       height = 8,
       width = 10)

# ------------------------------------------------------------------------------
#
#   Supplemental Figure SX (defunct barplot)
#
#   Note: This requires the list of top pairs to already have been generated
#         and saved to `output`.
#
# ------------------------------------------------------------------------------

# TBD pull in network figures from saved .svg (etc.)

top_pairs <- read.delim(file.path("output", "top_2.5_universal.tsv"),
                        header = TRUE,
                        sep = "\t")

same_fam_distro <- c()
diff_fam_distro <- c()
for(i in unique(top_pairs$tax_index1)) {
  same_fam <- which(data$taxonomy$family == data$taxonomy$family[i])
  diff_fam <- which(data$taxonomy$family != data$taxonomy$family[i])
  same_fam_distro <- c(same_fam_distro,
                       sum(top_pairs$tax_index1 == i & top_pairs$tax_index2 %in% same_fam))
  diff_fam_distro <- c(diff_fam_distro,
                       sum(top_pairs$tax_index1 == i & top_pairs$tax_index2 %in% diff_fam))
}

# Note: These results are just for taxa which have been resolved to at least the
# family level (112 / 223)!

rug_asv <- summarize_Sigmas(output_dir = "asv_days90_diet25_scale1")

fam1 <- sapply(rug_asv$tax_idx1, function(x) {
  data$taxonomy$family[x]
})
fam2 <- sapply(rug_asv$tax_idx2, function(x) {
  data$taxonomy$family[x]
})

x1 <- sum(fam1 == fam2, na.rm = TRUE)
x2 <- sum(fam1 != fam2, na.rm = TRUE)

cat(paste0("Overall proportion of family-family pairs: ",
           round(x1 / (x1 + x2), 3),
           "\n"))
cat(paste0("OBSERVED proportion of family-family pairs: ",
           round(sum(same_fam_distro) / (sum(same_fam_distro) + sum(diff_fam_distro)), 3),
           "\n"))
