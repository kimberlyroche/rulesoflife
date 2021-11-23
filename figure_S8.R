source("path_fix.R")

library(tidyverse)
library(rulesoflife)
# library(shapes)
library(frechet)
library(matrixsampling)
library(ggridges)

# These are calculated by `analysis_compute_Frechets.R`
# We ultimately need to roll that script's content into this one
F1 <- readRDS(file.path("output", "Frechet_1_corr.rds"))
F2 <- readRDS(file.path("output", "Frechet_2_corr.rds"))
F3 <- readRDS(file.path("output", "Frechet_3_corr.rds"))

# This function from the frechet package runs faster and gives a different (and
# frankly more intuitive) answer than the shapes package estimates, where means
# had a very very different scale than the samples they were computed from.
F1$mean <- CovFMean(F1$Sigmas)$Mout[[1]]
F2$mean <- CovFMean(F2$Sigmas)$Mout[[1]]
F3$mean <- CovFMean(F3$Sigmas)$Mout[[1]]

# Compute squared distances in each case as our variance analog
D <- dim(F1$Sigmas)[1]
N <- dim(F1$Sigmas)[3]
D1 <- numeric(N)
D2 <- numeric(N)
D3 <- numeric(N)
for(i in 1:N) {
  D1[i] <- distcov(F1$Sigmas[,,i] + diag(D)*1e-06,
                   F1$mean + diag(D)*1e-06,
                   method = "Riemannian")**2
  D2[i] <- distcov(F2$Sigmas[,,i] + diag(D)*1e-06,
                   F2$mean + diag(D)*1e-06,
                   method = "Riemannian")**2
  D3[i] <- distcov(F3$Sigmas[,,i] + diag(D)*1e-06,
                   F3$mean + diag(D)*1e-06,
                   method = "Riemannian")**2
}

# Plot
plot_df <- data.frame(x = c(D1,
                            D2,
                            D3),
                      label = c(rep("true", N),
                                rep("permuted", N),
                                rep("random", N)))
plot_df$label <- factor(plot_df$label, levels = c("random", "permuted", "true"))
levels(plot_df$label) <- c("Random", "Permutated", "Observed")

p <- ggplot(plot_df,
            aes(x = x, y = label, fill = label)) +
  geom_density_ridges(alpha = 0.8)  +
  theme_bw() +
  labs(x = "distance from mean dynamics",
       y = "",
       fill = "Data source") +
  theme(legend.position = "none")

ggsave(file.path("output", "figures", "S8.svg"),
       dpi = 100,
       units = "in",
       height = 4,
       width = 6)

v1 <- sum(D1) / (N-1)
v2 <- sum(D2) / (N-1)
v3 <- sum(D3) / (N-1)

cat(paste0("Ratio of variance observed to permuted: ", round(v1 / v2, 3), ":1\n"))
cat(paste0("Ratio of variance observed to random: ", round(v1 / v3, 3), ":1\n"))
