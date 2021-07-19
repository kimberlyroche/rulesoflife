source("path_fix.R")

library(rulesoflife)
library(tidyverse)
library(RColorBrewer)

source("ggplot_fix.R")

output_dir <- "asv_days90_diet25_scale1_scrambled"

distro_filename <- file.path("output", "scrambled_within_between_distros_diet-none.rds")
if(!file.exists(distro_filename)) {
  stop(paste0("File not found: ", distro_filename, "\n"))
}
distros <- readRDS(distro_filename)

# rug_filename <- file.path("output", "rug_asv.rds")
# if(file.exists(rug_filename)) {
#   rug_obj <- readRDS(rug_filename)
# } else {
  rug_obj <- summarize_Sigmas(output_dir)
  saveRDS(rug_obj, file = rug_filename)
# }
universalities <- apply(rug_obj$rug, 2, calc_universality_score)
consensus_signs <- apply(rug_obj$rug, 2, calc_consensus_sign)

selected_hosts <- readRDS(file.path("output", "overlapped_hosts.rds"))

distros <- distros %>%
  left_join(data.frame(pair = 1:ncol(rug_obj$rug),
                       score = universalities,
                       sign = consensus_signs), by = "pair")

# Combine between distributions
distros$type <- factor(distros$type)
levels(distros$type) <- c("between", "between", "within")

plot_df <- distros %>%
  group_by(pair, type) %>%
  mutate(mean_corr = mean(correlation)) %>%
  ungroup() %>%
  distinct(pair, type, mean_corr, score, sign) %>%
  pivot_wider(c(pair, score, sign), names_from = type, values_from = mean_corr)

p <- ggplot(plot_df, aes(x = between, y = within, fill = score)) +
  geom_point(size = 3, shape = 21) +
  scale_fill_distiller(palette = "Spectral") +
  labs(fill = "Universality\nscore",
       x = "mean ASV series correlation between hosts",
       y = "mean ASV series correlation within hosts")
# ggsave(file.path(plot_dir, paste0("within-between_all.png")),
#        plot = p,
#        units = "in",
#        dpi = 100,
#        height = 6,
#        width = 8)
show(p)

p <- ggplot(plot_df, aes(x = between, y = score, fill = score)) +
  geom_point(size = 3, shape = 21) +
  scale_fill_distiller(palette = "Spectral") +
  labs(fill = "score",
       x = "mean ASV series correlation between hosts (synchrony)",
       y = "universality score")
ggsave(file.path(plot_dir, paste0("universality-between_all.png")),
       plot = p,
       units = "in",
       dpi = 100,
       height = 6,
       width = 8)
show(p)

# How does the within-host correlation compare to the estimates of Sigma from
# the rug? Ans: really well.
# interp_cor <- distros %>%
#   group_by(pair) %>%
#   filter(type == "within") %>%
#   mutate(mean_within = mean(correlation)) %>%
#   distinct(pair, mean_within) %>%
#   arrange(pair) %>%
#   pull(mean_within)
#
# rug_cor <- colMeans(rug_obj$rug[rug_obj$hosts %in% selected_hosts,])
#
# ggplot(data.frame(x = rug_cor, y = interp_cor), aes(x = x, y = y)) +
#   geom_point(size = 2, shape = 21, fill = "#bbbbbb") +
#   xlim(c(-1, 1)) +
#   ylim(c(-1, 1)) +
#   labs(x = "pairwise correlation from model estimate (rug)",
#        y = "pairwise correlation from interpolated series")
