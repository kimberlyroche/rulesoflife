library(rulesoflife)
library(driver)
library(tidyr)
library(ggplot2)

# PCA plot components (over filtered CLR ASVs)

data_asv <- load_data(tax_level = "ASV")
clr_counts <- clr_array(data_asv$counts + 0.5, parts = 1)
PCA <- prcomp(t(clr_counts))
PoV <- PCA$sdev^2 / sum(PCA$sdev^2)

# Plot, labeling for year

labels <- as.numeric(unlist(sapply(data$metadata$collection_date, function(x) {
  substring(x, 1, 4)
})))
labels <- factor(labels)

palette <- generate_palette(hex_min = "#00AFBB",
                            hex_mid = "#BBBBBB",
                            hex_max = "#FC4E07",
                            S = length(levels(labels)))

plot_df <- data.frame(x = PCA$x[,1],
                      y = PCA$x[,2],
                      label = labels)

ggplot(plot_df, aes(x = x, y = y, color = factor(label))) +
  geom_point(size = 2) +
  scale_color_manual(values = palette) +
  labs(title = "PCA BABOONS",
         color = "Year",
         x = paste0("PC1 (", round(PoV[1], 3)*100 ,"% variance expl.)"),
         y = paste0("PC2 (", round(PoV[2], 3)*100 ,"% variance expl.)"))
plot_dir <- check_dir(c("output", "figures"))
ggsave(file.path(plot_dir, "PCA_year.png"),
       units = "in",
       dpi = 100,
       height = 5,
       width = 7)

# Plot, labeling for host

labels <- data_asv$metadata$sname
ref_hosts <- get_reference_hosts()$sname
labels[!(labels %in% ref_hosts)] <- NA
labels <- factor(labels)

plot_df <- data.frame(x = PCA$x[,1],
                      y = PCA$x[,2],
                      label = labels)

ggplot() +
  geom_point(data = plot_df[is.na(plot_df$label),],
             aes(x = x, y = y),
             color = "#BBBBBB",
             size = 2) +
  geom_point(data = plot_df[!is.na(plot_df$label),],
             aes(x = x, y = y, color = label),
             size = 2) +
  labs(title = "PCA BABOONS",
       color = "Host",
       x = paste0("PC1 (", round(PoV[1], 3)*100 ,"% variance expl.)"),
       y = paste0("PC2 (", round(PoV[2], 3)*100 ,"% variance expl.)"))
plot_dir <- check_dir(c("output", "figures"))
ggsave(file.path(plot_dir, "PCA_host.png"),
       units = "in",
       dpi = 100,
       height = 5,
       width = 7)

# Plot relative abundances

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
