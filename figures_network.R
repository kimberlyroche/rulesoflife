source("path_fix.R")

library(rulesoflife)
library(tidyverse)
library(ggraph)
library(igraph)

source("ggplot_fix.R")

# ------------------------------------------------------------------------------
#   Pull high confidence taxon pairs
# ------------------------------------------------------------------------------

data <- load_data(tax_level = "ASV")
plot_dir <- check_dir(c("output", "figures"))
rug_obj <- summarize_Sigmas(output_dir = "asv_days90_diet25_scale1")
scores <- apply(rug_obj$rug, 2, calc_universality_score)
consensus_signs <- apply(rug_obj$rug, 2, calc_consensus_sign)

palette_file <- file.path("output", "family_palette.rds")
if(file.exists(palette_file)) {
  family_palette <- readRDS(palette_file)
} else {
  unique_families <- unique(data$taxonomy[,6])
  family_palette <- generate_highcontrast_palette(length(unique_families))
  names(family_palette) <- unique_families
  names(family_palette)[is.na(names(family_palette))] <- "Unknown"
  family_palette[names(family_palette) == "Unknown"] <- "#999999"
  saveRDS(family_palette, file = file.path("output", "family_palette.rds"))
}

percents <- c(2.5, 1)
for(percent in percents) {
  k <- round(length(scores)*(percent/100))
  top_pairs <- order(scores, decreasing = TRUE)[1:k]

  pair_idx1 <- rug_obj$tax_idx1[top_pairs]
  pair_idx2 <- rug_obj$tax_idx2[top_pairs]

  map_df <- data.frame(old_idx = unique(c(pair_idx1, pair_idx2)))
  map_df$name <- 1:nrow(map_df)

  # ------------------------------------------------------------------------------
  #   Build node and edge data.frames
  # ------------------------------------------------------------------------------

  node_df <- map_df
  node_df$family <- sapply(node_df$old_idx, function(x) {
    data$taxonomy[x,6]
  })
  node_df <- node_df[,2:3]

  edge1 <- data.frame(old_idx = pair_idx1)
  edge1 <- left_join(edge1, map_df, by = "old_idx")$name
  edge2 <- data.frame(old_idx = pair_idx2)
  edge2 <- left_join(edge2, map_df, by = "old_idx")$name

  edge_df <- data.frame(from = edge1,
                        to = edge2,
                        sign = factor(consensus_signs[top_pairs]),
                        score = scores[top_pairs])
  levels(edge_df$sign) <- c("negative", "positive")

  # ------------------------------------------------------------------------------
  #   Build and plot graph object
  # ------------------------------------------------------------------------------

  na_idx <- which(is.na(node_df$family))
  node_df$family[na_idx] <- "Unknown"

  graph <- graph_from_data_frame(edge_df, node_df, directed = FALSE)

  # Not specifying the layout - defaults to "auto"
  # fr and kk layouts are ok here
  p <- ggraph(graph, layout = "fr") +
    geom_edge_link(aes(color = sign), width = 2, alpha = 1) +
    geom_node_point(aes(color = family), size = 5) +
    geom_node_label(aes(label = name), size = 3, repel = TRUE) +
    scale_colour_manual(values = family_palette) +
    scale_edge_colour_manual(values = c(negative = "gray", positive = "black")) +
    labs(title = paste0("Top ", percent, "% of pairs"),
         x = "dimension 1",
         y = "dimension 2") +
    theme_bw()
  ggsave(file.path(plot_dir, paste0("network_top_", percent, ".png")),
         p,
         units = "in",
         dpi = 100,
         height = 6,
         width = 8)
  show(p)
}

# Improvements:
# 3) Render table of taxon IDs
# 4) Eliminate "islands"; what criteria?

