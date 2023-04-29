# Universality score thresholds. These are determined by a permutation test.
# thresholds_scores <- data.frame(type = factor(c("Phylum", "Family", "ASV")),
#                                 x0 = c(0.162, 0.140, 0.136))

thresholds_scores <- data.frame(type = factor(c("Phylum", "Family", "ASV")),
                                x0 = c(0.176, 0.154, 0.174))

# Correlation thresholds.
# thresholds <- data.frame(type = factor(c("Phylum", "Family", "ASV")),
#                          lower = c(-0.303, -0.256, -0.263),
#                          upper = c(0.149, 0.207, 0.254))

# Correlation thresholds.
thresholds <- data.frame(type = factor(c("Phylum", "Family", "ASV")),
                         lower = c(-0.370, -0.329, -0.342),
                         upper = c(0.217, 0.279, 0.328))
