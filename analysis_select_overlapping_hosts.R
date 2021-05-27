source("path_fix.R")

library(rulesoflife)
library(tidyverse)
library(driver)

source("ggplot_fix.R")

data <- load_data(tax_level = "ASV")
metadata <- data$metadata

# Get max and min day per host, relative to a common baseline
common_baseline <- min(metadata$collection_date)
all_hosts <- unique(metadata$sname)
host_days <- data.frame(min_day = c(),
                        max_day = c(),
                        median_day = c(),
                        host = c())
for(host in all_hosts) {
  days <- sapply(metadata[metadata$sname == host,]$collection_date, function(x) {
    round(difftime(x, common_baseline, units = "days"))
  })
  host_days <- rbind(host_days,
                     data.frame(min_day = min(days),
                                max_day = max(days),
                                median_day = min(days) + (max(days) - min(days))/2,
                                host = host))
}

host_days <- host_days %>%
  arrange(median_day)
host_days$index <- 1:nrow(host_days)

host_days$type <- factor(sapply(1:nrow(host_days), function(i) {
  host_days$min_day[i] <= 750 && host_days$max_day[i] >= 4000
}))
levels(host_days$type) <- c("excluded", "included")

selected_hosts <- host_days[host_days$type == "included",]$host
cat("No. selected hosts:", length(selected_hosts), "\n")

saveRDS(selected_hosts, file.path("output", "overlapped_hosts.rds"))

# ------------------------------------------------------------------------------
#   Visualize overlap in hosts; select a subset with good overlap in time
# ------------------------------------------------------------------------------

if(TRUE) {
  p <- ggplot(host_days, aes(x = min_day,
                             xend = max_day,
                             y = index,
                             yend = index,
                             color = type)) +
    geom_segment(size = 1) +
    scale_color_manual(values = c("#aaaaaa", "#000000")) +
    scale_y_continuous(breaks = host_days$index,
                       labels = host_days$host) +
    theme(axis.text.y = element_text(size = 6)) +
    labs(x = "sampled days",
         y = "host",
         color = "Host status")
  ggsave(file.path(plot_dir, "host_day_overlap.png"),
         plot = p,
         units = "in",
         dpi = 100,
         height = 5,
         width = 8)
  show(p)
}
