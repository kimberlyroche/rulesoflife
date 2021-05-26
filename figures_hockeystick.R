source("path_fix.R")

library(rulesoflife)
library(tidyverse)

source("ggplot_fix.R")

# ------------------------------------------------------------------------------
#   Plot violin-style "hockeystick" plots separated for +/- associations
# ------------------------------------------------------------------------------

data <- load_data(tax_level = "ASV")
plot_dir <- check_dir(c("output", "figures"))
rug_obj <- summarize_Sigmas(output_dir = "asv_days90_diet25_scale1")
scores_pieces <- apply(rug_obj$rug, 2, function(x) {
  calc_universality_score(x, return_pieces = TRUE)
})
consensus_signs <- apply(rug_obj$rug, 2, calc_consensus_sign)

p1_factor <- cut(scores_pieces[1,], breaks = c(0, 0.6, 0.7, 0.8, 0.9, Inf))
levels(p1_factor) <- c("<60%", "60-70%", "70-80%", "80-90%", ">90%")
plot_df <- data.frame(p1 = scores_pieces[1,],
                      p1_factor = p1_factor,
                      p2 = scores_pieces[2,],
                      sign = consensus_signs)
plot_df$sign <- factor(plot_df$sign)
levels(plot_df$sign) <- c("consensus negative assoc.",
                          "omit",
                          "consensus positive assoc.")
p <- ggplot(plot_df %>% filter(sign != "omit"), aes(x = p1_factor, y = p2)) +
  geom_violin() +
  facet_wrap(. ~ sign) +
  labs(x = "proportion agreement in sign across hosts",
       y = "average absolute correlation")
ggsave(file.path(plot_dir, "piecewise_universality_scores.png"),
       p,
       units = "in",
       dpi = 100,
       height = 4,
       width = 7)
show(p)

# ------------------------------------------------------------------------------
#   Combine sign of associations and overlay strongly and weakly universal
#   example taxa
# ------------------------------------------------------------------------------

plot_df$tax1 <- rug_obj$tax_idx1
plot_df$tax2 <- rug_obj$tax_idx2

ranked_taxa <- plot_df %>%
  group_by(tax1) %>%
  mutate(mean_score = mean(p1*p2)) %>%
  arrange(desc(mean_score)) %>%
  select(tax1, mean_score) %>%
  distinct()

# Plot a taxon with HIGH universality
p <- ggplot(plot_df, aes(x = p1_factor, y = p2)) +
  geom_violin() +
  geom_jitter(data = plot_df %>% filter(tax1 == ranked_taxa$tax1[1]),
              mapping = aes(x = p1_factor, y = p2),
              width = 0.1,
              size = 3,
              shape = 21,
              fill = "#eb5721") +
  labs(x = "proportion agreement in sign across hosts",
       y = "average absolute correlation",
       title = paste0("Scores associated with ",
                      get_tax_label(data$taxonomy,
                                    ranked_taxa$tax1[1],
                                    "clr"))) +
  theme(plot.title = element_text(size = 9))
ggsave(file.path(plot_dir, "piecewise_universality_scores_strong-ex.png"),
       p,
       units = "in",
       dpi = 100,
       height = 4,
       width = 5)
show(p)

# Plot a taxon with LOW universality
p <- ggplot(plot_df, aes(x = p1_factor, y = p2)) +
  geom_violin() +
  geom_jitter(data = plot_df %>% filter(tax1 == ranked_taxa$tax1[nrow(ranked_taxa)]),
              mapping = aes(x = p1_factor, y = p2),
              width = 0.1,
              size = 3,
              shape = 21,
              fill = "#2cb2f5") +
  labs(x = "proportion agreement in sign across hosts",
       y = "average absolute correlation",
       title = paste0("Scores associated with ",
                      get_tax_label(data$taxonomy,
                                    ranked_taxa$tax1[nrow(ranked_taxa)],
                                    "clr"))) +
  theme(plot.title = element_text(size = 9))
ggsave(file.path(plot_dir, "piecewise_universality_scores_weak-ex.png"),
       p,
       units = "in",
       dpi = 100,
       height = 4,
       width = 5)
show(p)
