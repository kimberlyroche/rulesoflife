library(rulesoflife)
library(tidyverse)

output_dir <- "asv_days90_diet25_scale1"
plot_dir <- file.path("output", "figures")
data <- load_data(tax_level = "ASV")

# ------------------------------------------------------------------------------
#   Get universality distribution from canonical ASV rug
# ------------------------------------------------------------------------------

rug_obj <- summarize_Sigmas(output_dir)
scores <- apply(rug_obj$rug, 2, calc_universality_score)

# ------------------------------------------------------------------------------
#   Get universality distribution from scrambled ASV rug
# ------------------------------------------------------------------------------

scrambled_rug <- t(apply(rug_obj$rug, 1, sample))
scores_scrambled <- apply(scrambled_rug, 2, calc_universality_score)

# ------------------------------------------------------------------------------
#   Plot these distributions together
# ------------------------------------------------------------------------------

plot_df <- data.frame(score = c(scores, scores_scrambled),
                      type = c(rep("observed", length(scores)),
                               rep("null", length(scores_scrambled))))
ggplot() +
  geom_histogram(data = plot_df[plot_df$type != "null",],
                 aes(x = score, fill = type),
                 alpha = 0.5,
                 color = "white") +
  geom_histogram(data = plot_df[plot_df$type == "null",],
                 aes(x = score, fill = type),
                 alpha = 0.5,
                 color = "white") +
  scale_fill_manual(values = c(null = "red", observed = "navy")) +
  labs(fill = "Distribution")
ggsave(file.path(plot_dir, "universality_distributions.png"),
       units = "in",
       dpi = 100,
       height = 4,
       width = 5)

# ------------------------------------------------------------------------------
#   Defining a 95% cutoff
# ------------------------------------------------------------------------------

cutoff <- quantile(scores_scrambled, probs = c(0.95))
cat(paste0("95% cutoff at: ", round(cutoff, 3), "\n"))

# How many hits?
cat(paste0("Scores above 95% cutoff: ",
           round(sum(scores > cutoff) / length(scores), 3)*100,
           "%\n"))

# ------------------------------------------------------------------------------
#   Plotting weakest the weakest "hit"
# ------------------------------------------------------------------------------

scores_df <- data.frame(index = 1:length(scores), score = scores) %>%
  mutate(above_cutoff = ifelse(score > cutoff,
                               TRUE,
                               FALSE)) %>%
  filter(above_cutoff == TRUE)
lowest_passing <- scores_df %>%
  arrange(score) %>%
  slice(1) %>%
  pull(index)

tax_idx1 <- rug_obj$tax_idx1[lowest_passing]
tax_idx2 <- rug_obj$tax_idx2[lowest_passing]

# Note: this takes ~10 min.!
plot_aligned_trajectories(output_dir,
                          tax_idx1 = tax_idx1,
                          tax_idx2 = tax_idx2,
                          tax_label1 = get_tax_label(data$taxonomy, tax_idx1, "clr"),
                          tax_label2 = get_tax_label(data$taxonomy, tax_idx2, "clr"),
                          metadata = data$metadata,
                          save_file = TRUE)

cat("This taxon has net sign:",
    sign(sum(sign(rug_obj$rug[,lowest_passing]))),
    "\n")

# ------------------------------------------------------------------------------
#   Print table of top 250 universal pairs
# ------------------------------------------------------------------------------

table_df <- scores_df %>%
  arrange(desc(score)) %>%
  slice(1:250) %>%
  select(index, score)
table_df$rank <- 1:250
table_df$net_sign <- sapply(table_df$index, function(x) {
  sign(sum(sign(rug_obj$rug[,x])))
})
table_df$idx1 <- sapply(table_df$index, function(x) {
  rug_obj$tax_idx1[x]
})
table_df$idx2 <- sapply(table_df$index, function(x) {
  rug_obj$tax_idx2[x]
})
table_df$taxonomy_1 <- sapply(table_df$idx1, function(x) {
  tax_levels <- colnames(data$taxonomy)[2:ncol(data$taxonomy)]
  paste(paste(tax_levels,
              data$taxonomy[x,2:ncol(data$taxonomy)]), collapse = " / ")
})
table_df$sequence_1 <- sapply(table_df$idx1, function(x) {
  data$taxonomy[x,1]
})
table_df$taxonomy_2 <- sapply(table_df$idx2, function(x) {
  tax_levels <- colnames(data$taxonomy)[2:ncol(data$taxonomy)]
  paste(paste(tax_levels,
              data$taxonomy[x,2:ncol(data$taxonomy)]), collapse = " / ")
})
table_df$sequence_2 <- sapply(table_df$idx2, function(x) {
  data$taxonomy[x,1]
})

write.table(table_df %>% select(rank, score, net_sign, taxonomy_1, taxonomy_2, sequence_1, sequence_2),
            file = file.path("output", "top100_universal.tsv"),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

table(table_df$net_sign) / 250

