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

ggsave("output/figures/temp.svg",
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
