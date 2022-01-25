source("path_fix.R")

library(tidyverse)
library(rulesoflife)
library(driver)
library(cowplot)

# ------------------------------------------------------------------------------
#
#   Supplemental Figure S2 - relative abundance time courses for all 56 hosts
#                            and PCA plots of samples and selected hosts from
#                            2004-2009
#
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
#   Timecourse plots
# ------------------------------------------------------------------------------

data <- load_data(tax_level = "family")

relative_abundances <- data$counts
for(j in 1:ncol(relative_abundances)) {
  relative_abundances[,j] <- relative_abundances[,j]/sum(relative_abundances[,j])
}

# Collapse rare stuff for visualization
retain_taxa <- which(apply(relative_abundances[1:(nrow(relative_abundances)-1),], 1, mean) >= 0.025)
collapse_taxa <- setdiff(1:nrow(relative_abundances), retain_taxa)
trimmed_relab <- relative_abundances[retain_taxa,]
trimmed_relab <- rbind(trimmed_relab,
                       colSums(relative_abundances[collapse_taxa,]))

labels <- character(length(retain_taxa))
for(i in 1:length(retain_taxa)) {
  taxon_idx <- retain_taxa[i]
  level <- max(which(!is.na(data$taxonomy[taxon_idx,])))
  labels[i] <- paste0(colnames(data$taxonomy)[level],
                      " ",
                      data$taxonomy[taxon_idx,level])
}
labels <- c(labels, "rare phyla")

# palette_fn <- file.path("output", "timecourse_palette.rds")
palette_fn <- file.path("output", "family_palette.rds")
palette <- readRDS(palette_fn)

# Supplement this palette with a few Orders, odds and ends
palette2 <- palette[names(palette) %in% c("Bifidobacteriaceae",
                                          "Clostridiaceae 1",
                                          "Lachnospiraceae",
                                          "Prevotellaceae",
                                          "Ruminococcaceae",
                                          "Unknown",
                                          "Veillonellaceae")]
names(palette2)[2] <- "Clostridiales"
names(palette2)[6] <- "Other taxa"
palette2[["Rikenellaceae"]] <- "#F0B2EB"
palette2[["WCHB1-41"]] <- "#70CDC5"
# palette2[["Other taxa"]] <- "#A68DC8"

palette2 <- palette2[c(1:5,7:9,6)]

# Map taxon names to values in the family-level palette
mapping <- data.frame(taxon = c("domain Bacteria",
                                "family Bifidobacteriaceae",
                                "family Lachnospiraceae",
                                "family Prevotellaceae",
                                "family Rikenellaceae",
                                "family Ruminococcaceae",
                                "family Veillonellaceae",
                                "order Clostridiales",
                                "order WCHB1-41",
                                "rare phyla"),
                      palette_value = c("Unknown",
                                        "Bifidobacteriaceae",
                                        "Lachnospiraceae",
                                        "Prevotellaceae",
                                        "Rikenellaceae",
                                        "Ruminococcaceae",
                                        "Veillonellaceae",
                                        "Clostridiales",
                                        "WCHB1-41",
                                        "Other taxa"))

# Get all hosts
ref_hosts <- unique(data$metadata$sname)

plots <- list()
legend <- NULL
for(host in ref_hosts) {
  host_relab <- trimmed_relab[,data$metadata$sname == host]

  # Downsample
  ds_idx <- round(seq(1, ncol(host_relab), length.out = min(ncol(host_relab), 20)))
  ds_idx[1] <- 1
  ds_idx[length(ds_idx)] <- ncol(host_relab)
  host_ds <- host_relab[,ds_idx]

  plot_df <- cbind(1:nrow(host_ds), as.data.frame(host_ds))
  colnames(plot_df) <- c("taxon", 1:(ncol(plot_df)-1))
  plot_df <- pivot_longer(plot_df, !taxon, names_to = "sample", values_to = "relative_abundance")
  plot_df$taxon <- factor(plot_df$taxon)
  plot_df$sample <- as.numeric(plot_df$sample)

  plot_df <- plot_df %>%
    left_join(data.frame(taxon = levels(plot_df$taxon), name = labels), by = "taxon")
  plot_df$taxon <- plot_df$name

  # Combine 'domain Bacteria' and 'rare phyla' for the sake of labeling
  temp <- plot_df %>%
    filter(!(taxon %in% c("domain Bacteria", "rare phyla")))
  temp2 <- plot_df %>%
    filter(taxon %in% c("domain Bacteria", "rare phyla")) %>%
    group_by(sample) %>%
    mutate(combined_ra = sum(relative_abundance)) %>%
    filter(taxon != "domain Bacteria") %>%
    ungroup() %>%
    select(-c(relative_abundance)) %>%
    mutate(relative_abundance = combined_ra) %>%
    select(-c(combined_ra))
  plot_df <- rbind(temp, temp2)

  plot_df <- plot_df %>%
    left_join(mapping, by = "taxon")

  plot_df$palette_value <- factor(plot_df$palette_value, levels = names(palette2))

  p <- ggplot(plot_df, aes(x = sample, y = relative_abundance, fill = palette_value)) +
    geom_area() +
    scale_fill_manual(values = palette2) +
    theme(legend.position = "bottom") +
    labs(fill = "")
  if(is.null(legend)) {
    legend <- get_legend(p)
  }
  p <- p +
    # geom_area(linetype = 1, size = 0.3, color = "black") +
    theme_nothing() +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme(plot.margin = margin(t = 10, r = , b = 0, l = 1))

  plots[[length(plots) + 1]] <- p
}

anon_labels <- read.delim(file.path("output", "host_labels.tsv"),
                          header = TRUE,
                          sep = "\t")
labels2 <- data.frame(sname = ref_hosts) %>%
  left_join(anon_labels, by = "sname")

# Rearrange by labels
p <- plot_grid(plotlist = plots[order(labels2$host_label)],
               nrow = 4,
               labels = labels2$host_label[order(labels2$host_label)],
               scale = 0.95,
               label_x = -0.1,
               label_y = 1,
               label_size = 8)

p1 <- plot_grid(p, legend, ncol = 1, rel_heights = c(1, 0.15))

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
  select(label, year)
labels <- factor(labels_years$label, levels = c(ref_hosts, "other"))
years <- as.numeric(labels_years$year)

plot_df <- data.frame(x = PCA$x[,1],
                      y = PCA$x[,2],
                      label = labels,
                      year = years)
plot_df <- plot_df[plot_df$year %in% year_range,]

p2 <- ggplot(plot_df) +
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

p1_padded <- plot_grid(p1, scale = 0.9)
p <- plot_grid(p1_padded, NULL, p2,
               ncol = 1,
               rel_heights = c(1, 0.05, 1),
               labels = c("A", "B"),
               label_size = 18,
               label_x = 0,
               label_y = 1,
               scale = 0.95)

ggsave(file.path("output", "figures", "S2.png"),
       p,
       units = "in",
       dpi = 100,
       height = 10,
       width = 8)
