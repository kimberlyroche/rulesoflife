source("path_fix.R")

library(tidyverse)
library(rulesoflife)
library(driver)
library(cowplot)
library(RColorBrewer)
library(fido)
library(ape)
library(scales)

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

p1 <- ggplot(plot_df %>% filter(sign == "positive"), aes(x = d, y = score)) +
  geom_point(size = 2, shape = 21, fill = "#888888") +
  theme_bw() +
  labs(x = "phylogenetic distance",
       y = "universality score",
       fill = "Consensus\ncorrelation sign")

p1 <- ggplot(plot_df %>% filter(sign == "positive"), aes(x = d, y = score, fill = score)) +
  geom_point(size = 2, shape = 21) +
  xlim(c(min(plot_df$d), max(plot_df$d))) +
  ylim(c(min(plot_df$score), max(plot_df$score))) +
  scale_fill_gradient2(low = "white", high = "red") +
  theme_bw() +
  labs(x = "phylogenetic distance",
       y = "universality score",
       fill = "Consensus\ncorrelation sign") +
  theme(legend.position = "none")

p2 <- ggplot(plot_df %>% filter(sign == "negative"), aes(x = d, y = score, fill = score)) +
  geom_point(size = 2, shape = 21) +
  xlim(c(min(plot_df$d), max(plot_df$d))) +
  ylim(c(min(plot_df$score), max(plot_df$score))) +
  scale_fill_gradient2(low = "white", high = muted("navy")) +
  theme_bw() +
  labs(x = "phylogenetic distance",
       y = "universality score",
       fill = "Consensus\ncorrelation sign") +
  theme(legend.position = "none")

legend <- get_legend(ggplot(data.frame(x = 1:2, y = 1:2, sign = factor(c("1", "-1"), levels = c("1", "-1"))),
                            aes(x = x, y = y, fill = sign)) +
                       geom_point(size = 2, shape = 21) +
                       scale_fill_manual(values = c("red", muted("navy")), labels = c("positive", "negative")) +
                       theme_bw() +
                       labs(fill = "Consensus sign"))

# ------------------------------------------------------------------------------
#   Enrichment of closely related, strongly universal pairs
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

# ------------------------------------------------------------------------------
#   Family enrichment
# ------------------------------------------------------------------------------

enrichment <- NULL

# Calculate a baseline frequency for all family-family pairs
frequencies <- table(c(all_pairs_noNA$tax1, all_pairs_noNA$tax2))

# Stacked bar: plot observed relative representation of family-family pairs
all_pairs_noNA_tl <- all_pairs_noNA %>%
  filter(topleft == TRUE)
frequencies_subset <- table(c(all_pairs_noNA_tl$tax1, all_pairs_noNA_tl$tax2))

for(fam in names(frequencies_subset)) {
  fam_in_sample <- unname(unlist(frequencies_subset[fam]))
  sample_size <- unname(unlist(sum(frequencies_subset)))
  fam_in_bg <- unname(unlist(frequencies[fam]))
  bg_size <- unname(unlist(sum(frequencies)))
  # ctab <- matrix(c(fam_in_sample,
  #                  sample_size - fam_in_sample,
  #                  fam_in_bg,
  #                  bg_size - fam_in_bg),
  #                2, 2, byrow = TRUE)
  ctab <- matrix(c(fam_in_sample,
                   fam_in_bg - fam_in_sample,
                   sample_size - fam_in_sample,
                   bg_size - fam_in_bg),
                 2, 2, byrow = TRUE)
  prob <- fisher.test(ctab, alternative = "greater")$p.value
  enrichment <- rbind(enrichment,
                      data.frame(name = fam,
                                 type = "family",
                                 location = "Low phylogenetic distance, high median association strength",
                                 pvalue = prob,
                                 qvalue = NA))
}

# Multiple test correction
sel_idx <- which(enrichment$type == "family" & enrichment$location == "Low phylogenetic distance, high median association strength")
enrichment$qvalue[sel_idx] <- p.adjust(enrichment$pvalue[sel_idx], method = "BH")

signif <- c()
for(i in sel_idx) {
  q <- enrichment$qvalue[i]
  if(q < 0.05) {
    signif <- c(signif, enrichment$name[i])
    cat(paste0("ASV family: ", enrichment$name[i], ", adj. p-value: ", round(q, 3), "\n"))
  }
}

p3 <- plot_enrichment(frequencies_subset1 = frequencies_subset,
                      frequencies_subset2 = NULL,
                      frequencies = frequencies,
                      significant_families1 = signif,
                      significant_families2 = NULL,
                      plot_height = 6,
                      plot_width = 3,
                      legend_topmargin = 100,
                      use_pairs = FALSE,
                      rel_widths = c(1, 0.35, 1, 0.4, 2),
                      labels = c("overall\n", "low distance\nhigh univ."),
                      save_name = NULL,
                      suppress_y = TRUE)

# ------------------------------------------------------------------------------
#   Family-pair enrichment
# ------------------------------------------------------------------------------

frequencies <- table(all_pairs_noNA$taxpair)

frequencies_subset <- table(all_pairs_noNA$taxpair[all_pairs_noNA$topleft == TRUE])

for(fam in names(frequencies_subset)) {
  fam_in_sample <- unname(unlist(frequencies_subset[fam]))
  sample_size <- unname(unlist(sum(frequencies_subset)))
  fam_in_bg <- unname(unlist(frequencies[fam]))
  bg_size <- unname(unlist(sum(frequencies)))
  # ctab <- matrix(c(fam_in_sample,
  #                  sample_size - fam_in_sample,
  #                  fam_in_bg,
  #                  bg_size - fam_in_bg),
  #                2, 2, byrow = TRUE)
  ctab <- matrix(c(fam_in_sample,
                   fam_in_bg - fam_in_sample,
                   sample_size - fam_in_sample,
                   bg_size - fam_in_bg),
                 2, 2, byrow = TRUE)
  prob <- fisher.test(ctab, alternative = "greater")$p.value
  enrichment <- rbind(enrichment,
                      data.frame(name = fam,
                                 type = "family-pair",
                                 location = "Low phylogenetic distance, high median association strength",
                                 pvalue = prob,
                                 qvalue = NA))
}

# Multiple test correction
sel_idx <- which(enrichment$type == "family-pair" & enrichment$location == "Low phylogenetic distance, high median association strength")
enrichment$qvalue[sel_idx] <- p.adjust(enrichment$pvalue[sel_idx], method = "BH")

signif <- c()
for(i in sel_idx) {
  q <- enrichment$qvalue[i]
  if(q < 0.05) {
    signif <- c(signif, enrichment$name[i])
    cat(paste0("ASV family: ", enrichment$name[i], ", adj. p-value: ", round(q, 3), "\n"))
  }
}

p4 <- plot_enrichment(frequencies_subset1 = frequencies_subset,
                      frequencies_subset2 = NULL,
                      frequencies = frequencies,
                      significant_families1 = signif,
                      significant_families2 = NULL,
                      plot_height = 6,
                      plot_width = 6,
                      legend_topmargin = 100,
                      use_pairs = TRUE,
                      rel_widths = c(1, 0.3, 1, 1, 2.8),
                      labels = c("overall\n", "low distance\nhigh univ."),
                      save_name = NULL,
                      suppress_y = TRUE)

# ------------------------------------------------------------------------------
#   Plot all panels
# ------------------------------------------------------------------------------

common_scale <- 0.95

prow1 <- plot_grid(p1,
                   p2,
                   # p1 + ggtitle("Consensus positively associated ASVs"),
                   # p2 + ggtitle("Consensus negatively associated ASVs"),
                   legend,
                   ncol = 3,
               labels = c("A", "B", ""),
               label_size = 18,
               scale = common_scale,
               rel_widths = c(1, 1, 0.3))
prow2 <- plot_grid(NULL, p3, NULL, p4, NULL, ncol = 5,
               labels = c("", "C", "", "D", ""),
               label_size = 18,
               label_x = -0.03,
               label_y = 1.02,
               scale = common_scale,
               rel_widths = c(0.3, 0.9, 0.1, 1.1, 0.3))
p <- plot_grid(prow1, prow2, ncol = 1,
               rel_heights = c(1, 0.9))
ggsave(file.path("output", "figures", "phylogenetic.svg"),
       p,
       dpi = 100,
       units = "in",
       height = 8,
       width = 11)

# ------------------------------------------------------------------------------
#   Write out enrichment results
# ------------------------------------------------------------------------------

enrichment <- enrichment %>%
  arrange(location, type, name)
colnames(enrichment) <- c("ASV family or pair name",
                          "Type",
                          "Enrichment evaluated in",
                          "P-value (Fisher's exact test)",
                          "Adj. p-value (Benjamini-Hochberg)")
write.table(enrichment,
            file = file.path("output", "FigS8_table.tsv"),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

# ------------------------------------------------------------------------------
#   Evaluate trends: phylogenetic distance x universality score
# ------------------------------------------------------------------------------

overall <- lm(y ~ x, data.frame(x = plot_df$d, y = plot_df$score))
cat(paste0("Overall trend: ",
           round(summary(overall)$coef[2,1], 3),
           ", p-value: ",
           round(summary(overall)$coef[2,4], 3),
           "\n"))

pos <- plot_df %>%
  filter(sign == "positive")
pos <- lm(y ~ x, data.frame(x = pos$d, y = pos$score))
cat(paste0("Positive trend: ",
           round(summary(pos)$coef[2,1], 3),
           ", p-value: ",
           round(summary(pos)$coef[2,4], 3),
           "\n"))

# For positive pairs, we would expect a decrease in phylogenetic distance of
# 0.1 to increase universality by

neg <- plot_df %>%
  filter(sign == "negative")
neg <- lm(y ~ x, data.frame(x = neg$d, y = neg$score))
cat(paste0("Negative trend: ",
           round(summary(neg)$coef[2,1], 3),
           ", p-value: ",
           round(summary(neg)$coef[2,4], 3),
           "\n"))

# ------------------------------------------------------------------------------
#   Estimate effect sizes
# ------------------------------------------------------------------------------

test_df <- plot_df %>%
  filter(!is.na(sign))
test_df$sign <- factor(test_df$sign, levels = c("positive", "negative"))
levels(test_df$sign) <- c(0, 1)
lr_fit <- glm(sign ~ d, data = test_df, family = binomial(link = "logit"))
summary(lr_fit)

lo <- unname(coef(lr_fit)[2]) # Log odds associated with phylogenetic distance
exp(lo)
exp(lo)/(1+exp(lo))

cat(paste0("P-value for phylogenetic distance on sign: ",
           round(coef(summary(lr_fit))[2,4], 8), "\n"))
cat(paste0("Fold change in log odds of positive sign with unit increase in phylogenetic distance: ",
           round(lo, 3), "\n"))
cat(paste0("Odds of positive sign with unit increase in phylogenetic distance: ",
           round(exp(lo)/(1+exp(lo)), 3), "\n"))

