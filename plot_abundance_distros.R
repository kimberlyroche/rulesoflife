plot_abundance_distros <- function(rug_obj, clr.means, pairs, plot_label) {
  pair1 <- rug_obj$tax_idx1[pairs]
  pair2 <- rug_obj$tax_idx2[pairs]

  clr1 <- clr.means[pair1]
  clr2 <- clr.means[pair2]

  k <- length(pairs)
  max_pair <- sapply(1:k, function(x) {
    max(clr1[x], clr2[x])
  })
  min_pair <- sapply(1:k, function(x) {
    min(clr1[x], clr2[x])
  })

  ggplot() +
    # geom_histogram(data = data.frame(x = clr.means), mapping = aes(x = x), alpha = 0.2) +
    # geom_point(data = data.frame(x = c(min_pair, max_pair),
    #                              y = c(rep(1, k), rep(2, k)),
    #                              type = c(rep("Less abundant", k), rep("More abundant", k))),
    #            mapping = aes(x = x, y = y, color = type),
    #            size = 4,
    #            alpha = 0.5) +
    geom_density(data = data.frame(x = c(clr.means, min_pair, max_pair),
                                   type = c(rep("Background", length(clr.means)),
                                            rep("Less abundant", k),
                                            rep("More abundant", k))),
                 mapping = aes(x = x, fill = type),
                 alpha = 0.33) +
    scale_fill_manual(values = c("black", "red", "blue")) +
    labs(x = "mean CLR abundance",
         fill = "Identity of taxon in\ncorrelated pair",
         title = "Mean CLR abundance of selected taxon pairs again background")

  plot_dir <- check_dir(c("output", "figures"))
  ggsave(file.path(plot_dir, paste0("abundance_pairs_", plot_label, ".png")),
         units = "in",
         dpi = 100,
         height = 5,
         width = 7)
}
