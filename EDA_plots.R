library(rulesoflife)
library(driver)
library(tidyr)
library(ggplot2)

# PCA plot components (over filtered CLR ASVs)

plot_dir <- check_dir(c("output", "figures"))

data_asv <- load_data(tax_level = "ASV")
clr_counts <- clr_array(data_asv$counts + 0.5, parts = 1)
PCA <- prcomp(t(clr_counts))
PoV <- PCA$sdev^2 / sum(PCA$sdev^2)

# Visualize the number of samples per-host, per-year
# ref_hosts <- get_reference_hosts()$sname
# plot_df <- data$metadata %>%
#   filter(sname %in% ref_hosts) %>%
#   select(sname, collection_date) %>%
#   mutate(year = substr(collection_date, 1, 4)) %>%
#   count(sname, year)
# ggplot(plot_df, aes(x = year, y = n, color = factor(sname))) +
#   geom_point(size = 4, alpha = 0.5)

rgb_colors <- matrix(c(240, 7, 7,
                       227, 125, 16,
                       68, 194, 10,
                       90, 141, 242,
                       200, 90, 242), 5, 3, byrow = TRUE)/255
hex_colors <- apply(rgb_colors, 1, function(x) {
  rgb(x[1], x[2], x[3])
})

# ------------------------------------------------------------------------------
#   VERSION 1: One-host (DUI), several years
# ------------------------------------------------------------------------------

year_range <- seq(from = 2004, to = 2008, by = 1)

labels <- data$metadata %>%
  mutate(year = substr(collection_date, 1, 4)) %>%
  mutate(label = ifelse(sname == "DUI" & year %in% year_range,
                        year,
                        NA)) %>%
  pull(label)
labels <- factor(labels)
labels

plot_df <- data.frame(x = PCA$x[,1],
                      y = PCA$x[,2],
                      label = labels)

ggplot() +
  geom_point(data = plot_df[is.na(plot_df$label),],
             aes(x = x, y = y),
             color = "#bbbbbb",
             size = 2) +
  geom_point(data = plot_df[!is.na(plot_df$label),],
             aes(x = x, y = y, color = label),
             size = 3) +
  scale_color_manual(values = hex_colors) +
  labs(title = "Duiker's samples (2004-2008)",
       color = "Year",
       x = paste0("PC1 (", round(PoV[1], 3)*100 ,"% variance expl.)"),
       y = paste0("PC2 (", round(PoV[2], 3)*100 ,"% variance expl.)"))
ggsave(file.path(plot_dir, "PCA_1.png"),
       units = "in",
       dpi = 100,
       height = 6,
       width = 8)

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

ggplot(plot_df) +
  geom_point(data = plot_df[plot_df$label == "other",], mapping = aes(x = x, y = y, color = label),
             size = 2,
             color = "#bbbbbb") +
  geom_point(data = plot_df[plot_df$label != "other",], mapping = aes(x = x, y = y, color = label),
             size = 2) +
  scale_color_manual(values = hex_colors) +
  facet_wrap(. ~ year) +
  labs(title = "Selected host samples (2004-2009)",
       color = "Host",
       x = paste0("PC1 (", round(PoV[1], 3)*100 ,"% variance expl.)"),
       y = paste0("PC2 (", round(PoV[2], 3)*100 ,"% variance expl.)"))
ggsave(file.path(plot_dir, "PCA_2.png"),
       units = "in",
       dpi = 100,
       height = 7,
       width = 11)

# ------------------------------------------------------------------------------
#   Plot relative abundances
# ------------------------------------------------------------------------------

data <- load_data(tax_level = "family")
relative_abundances <- data$counts
for(j in 1:ncol(relative_abundances)) {
  relative_abundances[,j] <- relative_abundances[,j]/sum(relative_abundances[,j])
}

# Collapse rare stuff for visualization
retain_taxa <- which(apply(relative_abundances, 1, max) >= 0.2)
collapse_taxa <- setdiff(1:nrow(relative_abundances), retain_taxa)
trimmed_relab <- relative_abundances[retain_taxa,]
trimmed_relab <- rbind(trimmed_relab,
                       colSums(relative_abundances[collapse_taxa,]))
dim(trimmed_relab)

palette <- generate_highcontrast_palette(nrow(trimmed_relab))

ref_hosts <- get_reference_hosts()$sname
host <- ref_hosts[1]
host_relab <- trimmed_relab[,data$metadata$sname == host]

plot_df <- cbind(1:nrow(host_relab), as.data.frame(host_relab))
colnames(plot_df) <- c("taxon", 1:(ncol(plot_df)-1))
plot_df <- pivot_longer(plot_df, !taxon, names_to = "sample", values_to = "relative_abundance")
plot_df$taxon <- factor(plot_df$taxon)
plot_df$sample <- as.numeric(plot_df$sample)
head(plot_df)

ggplot(plot_df, aes(x = sample, y = relative_abundance, fill = taxon)) +
  geom_area() +
  scale_fill_manual(values = palette) +
  theme(legend.position = "none")
plot_dir <- check_dir(c("output", "figures"))
ggsave(file.path(plot_dir, paste0(host, "_shortseries.png")),
       units = "in",
       dpi = 100,
       height = 3,
       width = 3)
