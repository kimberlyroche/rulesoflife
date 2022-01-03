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

palette_fn <- file.path("output", "timecourse_palette.rds")
# palette_fn <- file.path("output", "family_palette.rds")
if(file.exists(palette_fn)) {
  palette <- readRDS(palette_fn)
} else {
  # Random palette
  getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
  palette <- sample(getPalette(nrow(trimmed_relab)))
  saveRDS(palette, palette_fn)
}

# Get all hosts
ref_hosts <- unique(data$metadata$sname)

# Get the top 20 best-sampled hosts
# ref_hosts <- sort(data$metadata %>%
#   group_by(sname) %>%
#   tally() %>%
#   arrange(desc(n)) %>%
#   slice(1:20) %>%
#   pull(sname))

# Get the best represented in the primary social groups
# ref_hosts <- c("DUI", "DUX", "LIW", "PEB", "VET")

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

  p <- ggplot(plot_df, aes(x = sample, y = relative_abundance, fill = taxon)) +
    geom_area() +
    scale_fill_manual(values = palette) +
    theme(legend.position = "bottom")
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

p <- plot_grid(plotlist = plots,
               # ncol = length(plots),
               nrow = 4,
               # labels = letters[1:length(plots)],
               labels = ref_hosts,
               scale = 0.95,
               label_x = -0.1,
               label_y = 1,
               label_size = 8)
pl <- plot_grid(p, legend, ncol = 1, rel_heights = c(1, 0.15))
ggsave(file.path(plot_dir, "F2b.svg"),
       pl,
       units = "in",
       dpi = 100,
       # height = 2,
       height = 6,
       width = 10)

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

ggsave("output/figures/S2.svg",
       plot = p,
       units = "in",
       dpi = 100,
       height = 6,
       width = 10)
