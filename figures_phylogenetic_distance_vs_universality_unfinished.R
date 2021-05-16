
universalities <- apply(rug_obj$rug, 2, calc_universality_score)
consensus_sign <- apply(rug_obj$rug, 2, consensus_signs)

plot_df <- data.frame(universality = universalities[canonical_col_order],
                      universality_sign = consensus_sign[canonical_col_order],
                      tax_distance = tax_distances[canonical_col_order])
plot_df$universality_sign <- factor(plot_df$universality_sign)

ggplot(plot_df, aes(x = tax_distance, y = universality, color = universality_sign)) +
  geom_point(alpha = 0.66) +
  scale_color_manual(values = c("#5680e3", "#bbbbbb", "#e35656")) +
  labs(x = "phylogenetic distance",
       color = "Consensus\nsign")
ggsave(file.path(plot_dir, "phylogenetic_distance_vs_universality.png"),
       units = "in",
       dpi = 100,
       height = 5,
       width = 7)
