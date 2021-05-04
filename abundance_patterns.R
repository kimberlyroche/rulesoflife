library(rulesoflife)
library(driver)
library(ggplot2)

# Pull data
rug_obj <- summarize_Sigmas(output_dir = "asv_days90_diet25_scale1")

# Get mean CLR abundance of all taxa
data <- load_data(tax_level = "ASV")
clr.counts <- clr_array(data$counts + 0.5, parts = 1)
clr.means <- rowMeans(clr.counts)

col_order <- order(colMeans(rug_obj$rug))
# plot_rug(rug_obj$rug, canonical_col_order = col_order)

k <- 100
pairs <- list(bottom = col_order[1:k],
              middle = col_order[ceiling(ncol(rug_obj$rug)/2-k/2):floor(ncol(rug_obj$rug)/2+k/2)],
              top = col_order[(length(col_order)-k+1):length(col_order)])
# There's a sort of keystone-ness going on in the top and bottom pairs (in terms
# of correlation): the same taxa occur frequently in these pairs.
# Optionally downsample
pairs[["middle"]] <- sample(pairs[["middle"]], size = 50)

for(use in names(pairs)) {
  use_k <- pairs[[use]]
  pair1 <- rug_obj$tax_idx1[use_k]
  pair2 <- rug_obj$tax_idx2[use_k]

  clr1 <- clr.means[pair1]
  clr2 <- clr.means[pair2]

  max_pair <- sapply(1:k, function(x) {
    max(clr1[x], clr2[x])
  })
  min_pair <- sapply(1:k, function(x) {
    min(clr1[x], clr2[x])
  })

  ggplot() +
    geom_histogram(data = data.frame(x = clr.means), mapping = aes(x = x), alpha = 0.2) +
    geom_point(data = data.frame(x = c(min_pair, max_pair),
                                 y = c(rep(1, k), rep(2, k)),
                                 type = c(rep("Less abundant", k), rep("More abundant", k))),
               mapping = aes(x = x, y = y, color = type),
               size = 4,
               alpha = 0.5) +
    labs(x = "mean CLR abundance",
         color = "Identity of taxon in\ncorrelated pair",
         title = "Mean CLR abundance of selected taxon pairs again background")

  plot_dir <- check_dir(c("output", "figures"))
  ggsave(file.path(plot_dir, paste0("abundance_pairs_", use, "_", k, ".png")),
         units = "in",
         dpi = 100,
         height = 5,
         width = 7)
}
