source("path_fix.R")

library(tidyverse)
library(rulesoflife)
library(driver)
library(cowplot)

# ------------------------------------------------------------------------------
#   PCA plots
# ------------------------------------------------------------------------------

# Reload ASV-level data!
data <- load_data(tax_level = "ASV")
clr_counts <- clr_array(data$counts + 0.5, parts = 1)

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
  dplyr::select(label, year)
labels <- factor(labels_years$label, levels = c(ref_hosts, "other"))
years <- as.numeric(labels_years$year)

plot_df <- data.frame(x = PCA$x[,1],
                      y = PCA$x[,2],
                      label = labels,
                      year = years)
plot_df <- plot_df[plot_df$year %in% year_range,]

plot_df$label <- factor(plot_df$label, levels = c("DUI", "DUX", "LIW", "PEB", "VET", "other"))
levels(plot_df$label) <- c("F09", "F20", "F35", "F27", "F31", "other")

# Label significantly predicts
# summary(aov(x ~ label, plot_df))
# summary(aov(y ~ label, plot_df))

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
  labs(fill = "Host",
       # title = "Selected host samples (2004-2009)",
       x = paste0("PC1 (", round(PoV[1], 3)*100 ,"% variance expl.)"),
       y = paste0("PC2 (", round(PoV[2], 3)*100 ,"% variance expl.)")) +
  theme_bw()

p

ggsave(file.path("output", "figures", "Figure_1_Supplement_2.png"),
       p,
       units = "in",
       dpi = 200,
       height = 6,
       width = 10,
       bg = "white")
