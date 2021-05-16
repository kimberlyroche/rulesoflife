library(rulesoflife)
library(ggraph)
library(igraph)

source("ggplot_fix.R")

data <- load_data(tax_level = "ASV")
plot_dir <- check_dir(c("output", "figures"))
rug_obj <- summarize_Sigmas(output_dir = "asv_days90_diet25_scale1")
scores <- apply(rug_obj$rug, 2, calc_universality_score)
consensus_signs <- apply(rug_obj$rug, 2, calc_consensus_sign)

# Pull top k features in terms of universality
k <- 100
top_pairs <- order(scores, decreasing = TRUE)[1:k]

pair_idx1 <- rug_obj$tax_idx1[top_pairs]
pair_idx2 <- rug_obj$tax_idx2[top_pairs]

map_df <- data.frame(old_idx = unique(c(pair_idx1, pair_idx2)))
map_df$name <- 1:nrow(map_df)

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

graph <- graph_from_data_frame(edge_df, node_df, directed = FALSE)

# Not specifying the layout - defaults to "auto"
# fr and kk layouts are ok here
p <- ggraph(graph, layout = "fr") +
  geom_edge_link(aes(color = sign, width = score), alpha = 0.5) +
  geom_node_point(size = 5) +
  geom_node_label(aes(label = family), size = 3, repel = TRUE) +
  scale_edge_colour_manual(values = c(negative = "blue", positive = "red")) +
  labs(title = paste0("Top ", k, " pairs"),
       x = "dimension 1",
       y = "dimension 2") +
  theme_bw() +
  theme(legend.position = "none")
ggsave(file.path(plot_dir, paste0("network_top", k, ".png")),
       p,
       units = "in",
       dpi = 100,
       height = 6,
       width = 8)
show(p)

# Improve this by replacing node labels with nodes colored for family (etc.)
