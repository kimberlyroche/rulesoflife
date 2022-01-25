source("path_fix.R")

library(tidyverse)
library(rulesoflife)
library(RColorBrewer)
library(cowplot)
library(magick)

# ------------------------------------------------------------------------------
#
#   Supplemental Figure 1 - overview of sampling scheme and timecourses for all
#                           baboon hosts
#
# ------------------------------------------------------------------------------

data <- load_data(tax_level = "family")

hosts_dates <- data$metadata %>%
  select(sname, collection_date)
hosts_dates$sname <- factor(hosts_dates$sname,
                            levels = sort(unique(hosts_dates$sname), decreasing = TRUE))
baseline_time <- min(hosts_dates$collection_date)
hosts_dates$time <- as.numeric(sapply(hosts_dates$collection_date, function(x) difftime(x, baseline_time, units = "days")))

xticks <- seq(from = 0, to = 5000, length = 20)
xlabs <- character(length(xticks))
for(i in 1:length(xticks)) {
  xlabs[i] <- as.character(as.Date(baseline_time) + xticks[i])
}

anon_labels <- read.delim(file.path("output", "host_labels.tsv"),
                          header = TRUE,
                          sep = "\t")

hosts_dates <- hosts_dates %>%
  left_join(anon_labels, by = "sname")
hosts_dates$host_label <- factor(hosts_dates$host_label)
levels(hosts_dates$host_label) <- rev(levels(hosts_dates$host_label))

p1 <- ggplot(hosts_dates, aes(x = time, y = host_label)) +
  geom_point(size = 1.5, shape = 21, fill = "#000000") +
  theme_bw() +
  labs(x = "sample collection date",
       y = "host") +
  scale_x_continuous(breaks = xticks, labels = xlabs) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
p1_padded <- plot_grid(NULL, p1, NULL, ncol = 3, rel_widths = c(0.05, 1, 0.2))

p2 <- ggdraw() +
  draw_image("output/figures/social_group_ranges.png")
p2_padded <- plot_grid(NULL, p2, NULL, ncol = 1, rel_heights = c(0.05, 1, 0.05))

p <- plot_grid(p1_padded, p2_padded, ncol = 1,
               labels = c("A", "B"),
               label_size = 20,
               label_x = 0,
               scale = 0.95)

ggsave(file.path("output", "figures", "S1.png"),
       p,
       units = "in",
       dpi = 100,
       height = 12,
       width = 8)
