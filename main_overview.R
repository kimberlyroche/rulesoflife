source("path_fix.R")

library(tidyverse)
library(rulesoflife)
library(driver)
library(cowplot)
library(RColorBrewer)

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
    dplyr::select(-c(relative_abundance)) %>%
    mutate(relative_abundance = combined_ra) %>%
    dplyr::select(-c(combined_ra))
  plot_df <- rbind(temp, temp2)

  plot_df <- plot_df %>%
    left_join(mapping, by = "taxon")

  plot_df$palette_value <- factor(plot_df$palette_value, levels = names(palette2))

  p <- ggplot(plot_df, aes(x = sample, y = relative_abundance, fill = palette_value)) +
    geom_area() +
    scale_fill_manual(values = palette2) +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 11)) +
    labs(fill = "")
  if(is.null(legend)) {
    legend <- get_legend(p)
  }
  p <- p +
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
               label_y = 1.03,
               label_size = 12)

p1 <- plot_grid(p, legend, ncol = 1, rel_heights = c(1, 0.15))

# ------------------------------------------------------------------------------
#   Sample alignment
# ------------------------------------------------------------------------------

hosts_dates <- data$metadata %>%
  dplyr::select(sname, collection_date)
hosts_dates$sname <- factor(hosts_dates$sname,
                            levels = sort(unique(hosts_dates$sname), decreasing = TRUE))
baseline_time <- min(hosts_dates$collection_date)
hosts_dates$time <- as.numeric(sapply(hosts_dates$collection_date, function(x) difftime(x, baseline_time, units = "days")))

# xticks <- seq(from = 0, to = 5000, length = 10)
# xlabs <- character(length(xticks))
# for(i in 1:length(xticks)) {
#   xlabs[i] <- as.character(as.Date(baseline_time) + xticks[i])
# }
xlabs <- c("2001", "2002", "2003", "2004", "2005",
           "2006", "2007", "2008", "2009", "2010",
           "2011", "2012", "2013")
xticks <- numeric(length(xlabs))
for(i in 1:length(xticks)) {
  xticks[i] <- as.numeric(as.Date(paste0(xlabs[i], "-01-01")) - as.Date(baseline_time))
}

anon_labels <- read.delim(file.path("output", "host_labels.tsv"),
                          header = TRUE,
                          sep = "\t")

hosts_dates <- hosts_dates %>%
  left_join(anon_labels, by = "sname")
hosts_dates$host_label <- factor(hosts_dates$host_label)
levels(hosts_dates$host_label) <- rev(levels(hosts_dates$host_label))

p2 <- ggplot(hosts_dates, aes(x = time, y = host_label)) +
  geom_point(size = 1.25, shape = 21, fill = "#000000") +
  theme_bw() +
  labs(x = "sample collection date",
       y = "host") +
  scale_x_continuous(breaks = xticks, labels = xlabs) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 14),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16))

p1_padded <- plot_grid(NULL, p1, ncol = 2,
                       rel_widths = c(0.04, 1))

prow2 <- plot_grid(p2,
                   p1_padded,
                   ncol = 2,
                   rel_widths = c(0.8, 1.2),
                   labels = c("B", "C"),
                   label_size = 20,
                   label_x = 0,
                   label_y = 1.02,
                   scale = 0.98)

prow1 <- plot_grid(ggdraw() +
                     draw_image(file.path("output", "figures", "Figure_1_include.png")),
                   ncol = 1,
                   labels = c("A"),
                   label_size = 20,
                   label_x = 0,
                   label_y = 1.02,
                   scale = 1.00)

p <- plot_grid(prow1, prow2, ncol = 1,
               rel_heights = c(1, 1))

ggsave(file.path("output", "figures", "Figure_1.png"),
       p,
       units = "in",
       dpi = 200,
       height = 10,
       width = 13,
       bg = "white")
