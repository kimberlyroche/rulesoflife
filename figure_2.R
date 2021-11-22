source("path_fix.R")

library(tidyverse)
library(rulesoflife)
library(RColorBrewer)
library(cowplot)

# ------------------------------------------------------------------------------
#
#   Figure 2 - timecourses from all 56 baboon hosts
#
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
ggsave(file.path(plot_dir, "F2.svg"),
       pl,
       units = "in",
       dpi = 100,
       # height = 2,
       height = 6,
       width = 10)
