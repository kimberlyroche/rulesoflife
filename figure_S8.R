source("path_fix.R")

library(tidyverse)
library(rulesoflife)
library(driver)
library(cowplot)
library(RColorBrewer)
library(fido)

# ------------------------------------------------------------------------------
#
#   Supplemental Figure S9 - scatterplots of phylogenetic distance vs. univer-
#                            sality score or median association strength, plus
#                            barplots of enrichment of closely related and
#                            strongly universal pairs
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
                      mas = median_assoc_strength,
                      tax_idx1 = rug_asv$tax_idx1,
                      tax_idx2 = rug_asv$tax_idx2)
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

overall <- lm(y ~ x, data.frame(x = plot_df$d, y = plot_df$score))
cat(paste0("Overall trend: ", round(summary(overall)$coef[2,1], 3), ", p-value: ", round(summary(overall)$coef[2,4], 3), "\n"))

pos <- plot_df %>%
  filter(sign == "positive")
pos <- lm(y ~ x, data.frame(x = pos$d, y = pos$score))
cat(paste0("Positive trend: ", round(summary(pos)$coef[2,1], 3), ", p-value: ", round(summary(pos)$coef[2,4], 3), "\n"))

neg <- plot_df %>%
  filter(sign == "negative")
neg <- lm(y ~ x, data.frame(x = neg$d, y = neg$score))
cat(paste0("Negative trend: ", round(summary(neg)$coef[2,1], 3), ", p-value: ", round(summary(neg)$coef[2,4], 3), "\n"))

ggsave("output/figures/SF6a.png",
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

ggsave("output/figures/SF6b.png",
       p,
       dpi = 100,
       units = "in",
       height = 5,
       width = 6.5)

# ------------------------------------------------------------------------------
#   Label the phylogenetic distance vs. median association strength plot by the
#   number of sequence mismatches
# ------------------------------------------------------------------------------

if(FALSE) {
  d <- sequence_distance(distance_type = "N")
  plot_df$mismatches <- NA
  for(i in 1:nrow(plot_df)) {
    plot_df$mismatches[i] <- d[plot_df$tax_idx1[i],
                               plot_df$tax_idx2[i]]
  }

  plot_df$mismatch_factor <- as.factor(plot_df$mismatches)
  show_max <- 8
  levels(plot_df$mismatch_factor) <- c(1:show_max,
                                       rep(paste0(show_max+1, "+"),
                                           length(levels(plot_df$mismatch_factor))-show_max))

  p <- ggplot(plot_df,
              aes(x = d, y = mas, fill = factor(mismatch_factor))) +
    geom_point(size = 3, shape = 21) +
    theme_bw() +
    scale_fill_manual(values = c(brewer.pal(show_max, "Spectral"), "#dddddd")) +
    labs(x = "phylogenetic distance",
         y = "median association strength",
         fill = "Sequence mismatches")

  ggsave("output/figures/SF6a-mismatches.png",
         p,
         dpi = 100,
         units = "in",
         height = 5,
         width = 8)

  plot_df$fam1 <- sapply(plot_df$tax_idx1, function(x) data$taxonomy[x,6])
  plot_df$fam2 <- sapply(plot_df$tax_idx2, function(x) data$taxonomy[x,6])

  plot_df %>%
    filter(mismatches < 10) %>%
    arrange(mismatches) %>%
    select(mismatches, fam1, fam2)
}

# ------------------------------------------------------------------------------
#   Enrichment of top left-hand pairs
# ------------------------------------------------------------------------------

# Score enrichment of family pairs or families themselves?
use_pairs <- FALSE

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

if(use_pairs) {
  frequencies <- table(all_pairs_noNA$taxpair)

  frequencies_subset <- table(all_pairs_noNA$taxpair[all_pairs_noNA$topleft == TRUE])

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
                  rel_widths = c(1, 0.35, 1, 0.3, 3),
                  labels = c("overall", "top left"),
                  save_name = "SF8-enrichment-pair.png")

} else {
  # Calculate a baseline frequency for all family-family pairs
  frequencies <- table(c(all_pairs_noNA$tax1, all_pairs_noNA$tax2))

  # Stacked bar: plot observed relative representation of family-family pairs
  all_pairs_noNA_tl <- all_pairs_noNA %>%
    filter(topleft == TRUE)
  frequencies_subset <- table(c(all_pairs_noNA_tl$tax1, all_pairs_noNA_tl$tax2))

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
                  labels = c("overall", "top left"),
                  save_name = "SF8-enrichment.png")

}
