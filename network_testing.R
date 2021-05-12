library(rulesoflife)
library(ggraph)
library(igraph)

head(highschool)

data <- load_data(tax_level = "ASV")

rug_obj <- summarize_Sigmas(output_dir = "asv_days90_diet25_scale1")
scores <- apply(rug_obj$rug, 2, calc_universality_score)
consensus_signs <- apply(rug_obj$rug, 2, calc_consensus_sign)

# Pull top k features in terms of universality
top_pairs <- order(scores, decreasing = TRUE)[1:100]

pair_idx1 <- rug_obj$tax_idx1[top_pairs]
pair_idx2 <- rug_obj$tax_idx2[top_pairs]

map_df <- data.frame(old_idx = unique(c(pair_idx1, pair_idx2)))
map_df$name <- 1:nrow(map_df)

node_df <- map_df
node_df$family <- sapply(node_df$old_idx, function(x) {
  data$taxonomy[x,6]
})
node_df <- node_df[,2:3]
str(node_df)

edge1 <- data.frame(old_idx = pair_idx1)
edge1 <- left_join(edge1, map_df, by = "old_idx")$name
edge2 <- data.frame(old_idx = pair_idx2)
edge2 <- left_join(edge2, map_df, by = "old_idx")$name

edge_df <- data.frame(from = edge1,
                      to = edge2)
str(edge_df)

# graph <- graph_from_data_frame(highschool)
graph <- graph_from_data_frame(edge_df, node_df, directed = FALSE)

# Not specifying the layout - defaults to "auto"
# fr and kk layouts are ok here
ggraph(graph, layout = "kk") +
  geom_edge_link(alpha = 0.5) +
  geom_node_point(size = 5) +
  geom_node_label(aes(label = family), repel = TRUE) +
  theme_bw()




plot_df_links <- data.frame(idx1 = rug_obj$tax_idx1[top_pairs],
                   idx2 = rug_obj$tax_idx2[top_pairs],
                   score = scores[top_pairs],
                   sign = consensus_signs[top_pairs])
plot_df_links$sign <- factor(plot_df_links$sign, levels = c("-1", "1"))
levels(plot_df_links$sign) <- c("negative", "positive")
head(plot_df_links)

plot_df_nodes <- data.frame(idx = unique(c(plot_df_links$idx1, plot_df_links$idx2)))
plot_df_nodes$label <- sapply(plot_df_nodes$idx, function(x) {
  # get_tax_label(data$taxonomy, x, coord_system_label = "clr")
  data$taxonomy[x,6]
})
head(plot_df_nodes)

# graph <- graph_from_data_frame(highschool)
graph <- graph_from_data_frame(plot_df_links, plot_df_nodes, directed = FALSE)

# Not specifying the layout - defaults to "auto"
# fr and kk layouts are ok here
ggraph(graph, layout = "kk") +
  geom_edge_link(aes(width = score, color = sign), alpha = 0.5) +
  geom_node_point(size = 5) +
  geom_node_label(aes(label = label), repel = TRUE) +
  scale_edge_colour_manual(values = c(negative = "blue", positive = "red")) +
  theme_bw()
