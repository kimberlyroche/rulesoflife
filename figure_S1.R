source("path_fix.R")

library(tidyverse)
library(rulesoflife)
library(RColorBrewer)
library(cowplot)

# ------------------------------------------------------------------------------
#
#   Supplemental Figure 1a - overview of sampling scheme
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

p <- ggplot(hosts_dates, aes(x = time, y = sname)) +
  geom_point() +
  theme_bw() +
  labs(x = "days from first sample",
       y = "host short name") +
  scale_x_continuous(breaks = xticks, labels = xlabs) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

ggsave(file.path(plot_dir, "S1a.svg"),
       p,
       units = "in",
       dpi = 100,
       height = 7,
       width = 8)
