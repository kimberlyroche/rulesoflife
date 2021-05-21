source("path_fix.R")

library(rulesoflife)
library(driver)
library(tidyr)
library(ggplot2)
library(cowplot)

source("ggplot_fix.R")

plot_dir <- check_dir(c("output", "figures"))

# ------------------------------------------------------------------------------
#   Generate PC embedding at ASV level; set host colors
# ------------------------------------------------------------------------------

# PCA plot components (over filtered CLR ASVs)

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

# ------------------------------------------------------------------------------
#   VERSION 1: Duiker's samples x 5 years
# ------------------------------------------------------------------------------

year_range <- seq(from = 2004, to = 2008, by = 1)

labels <- data$metadata %>%
  mutate(year = substr(collection_date, 1, 4)) %>%
  mutate(label = ifelse(sname == "DUI" & year %in% year_range,
                        year,
                        NA)) %>%
  pull(label)
labels <- factor(labels)

plot_df <- data.frame(x = PCA$x[,1],
                      y = PCA$x[,2],
                      label = labels)

p <- ggplot() +
  geom_point(data = plot_df[is.na(plot_df$label),],
             aes(x = x, y = y),
             shape = 21,
             size = 2,
             fill = "#bbbbbb",
             color = "#888888") +
  geom_point(data = plot_df[!is.na(plot_df$label),],
             aes(x = x, y = y, fill = label),
             size = 3,
             shape = 21) +
  scale_fill_manual(values = hex_colors) +
  labs(title = "Duiker's samples (2004-2008)",
       fill = "Year",
       x = paste0("PC1 (", round(PoV[1], 3)*100 ,"% variance expl.)"),
       y = paste0("PC2 (", round(PoV[2], 3)*100 ,"% variance expl.)"))
ggsave(file.path(plot_dir, "PCA_DUI.png"),
       plot = p,
       units = "in",
       dpi = 100,
       height = 6,
       width = 8)
show(p)

# ------------------------------------------------------------------------------
#   VERSION 2: Several hosts, years faceted
# ------------------------------------------------------------------------------

year_range <- seq(from = 2004, to = 2009, by = 1)

labels_years <- data$metadata %>%
  mutate(year = substr(collection_date, 1, 4)) %>%
  mutate(label = ifelse(sname %in% ref_hosts & year %in% year_range,
                        sname,
                        "other")) %>%
  select(label, year)
labels <- factor(labels_years$label, levels = c("DUI", "DUX", "LIW", "PEB", "VET", "other"))
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
       y = paste0("PC2 (", round(PoV[2], 3)*100 ,"% variance expl.)"))
ggsave(file.path(plot_dir, "PCA_5hosts.png"),
       plot = p,
       units = "in",
       dpi = 100,
       height = 7,
       width = 11,
       type = "cairo")
show(p)
