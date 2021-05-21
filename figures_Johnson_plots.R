source("path_fix.R")

library(rulesoflife)
library(driver)
library(tidyr)
library(ggplot2)
library(cowplot)

source("ggplot_fix.R")

plot_dir <- file.path(check_dir(c("output", "figures", "Johnson_plots")))

# ------------------------------------------------------------------------------
#   Pull and filter data to large-ish relative abundances (>1% avg.)
# ------------------------------------------------------------------------------

data <- load_data(tax_level = "family")
relative_abundances <- data$counts
for(j in 1:ncol(relative_abundances)) {
  relative_abundances[,j] <- relative_abundances[,j]/sum(relative_abundances[,j])
}

# Collapse rare stuff for visualization
retain_taxa <- which(apply(relative_abundances, 1, mean) >= 0.01)
collapse_taxa <- setdiff(1:nrow(relative_abundances), retain_taxa)
trimmed_relab <- relative_abundances[retain_taxa,]
trimmed_relab <- rbind(trimmed_relab,
                       colSums(relative_abundances[collapse_taxa,]))

saved_palette <- file.path(plot_dir, "saved_palette.rds")
if(file.exists(saved_palette)) {
  palette <- readRDS(saved_palette)
} else {
  palette <- generate_highcontrast_palette(nrow(trimmed_relab))
  saveRDS(palette, file = saved_palette)
}

# ------------------------------------------------------------------------------
#   Pull and filter data to large-ish relative abundances (>1% avg.)
# ------------------------------------------------------------------------------

ref_hosts <- unique(data$metadata$sname)

for(host in ref_hosts) {
  host_relab <- trimmed_relab[,data$metadata$sname == host]

  # Downsample
  ds_idx <- round(seq(1, ncol(host_relab), length.out = min(ncol(host_relab), 50)))
  ds_idx[1] <- 1
  ds_idx[length(ds_idx)] <- ncol(host_relab)
  host_ds <- host_relab[,ds_idx]

  plot_df <- cbind(1:nrow(host_ds), as.data.frame(host_ds))
  colnames(plot_df) <- c("taxon", 1:(ncol(plot_df)-1))
  plot_df <- pivot_longer(plot_df, !taxon, names_to = "sample", values_to = "relative_abundance")
  plot_df$taxon <- factor(plot_df$taxon)
  plot_df$sample <- as.numeric(plot_df$sample)
  head(plot_df)

  p <- ggplot(plot_df, aes(x = sample, y = relative_abundance, fill = taxon)) +
    geom_area() +
    scale_fill_manual(values = palette) +
    theme_nothing() +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0))
  ggsave(file.path(plot_dir, paste0(host, "_shortseries.pdf")),
         p,
         units = "in",
         dpi = 100,
         height = 4,
         width = 4)
}

# These figures need to be joined in Illustrator.
