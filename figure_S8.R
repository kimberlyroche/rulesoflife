source("path_fix.R")

library(rulesoflife)
library(shapes)
library(matrixsampling)

# These are calculated by `analysis_compute_Frechets.R`
# We ultimately need to roll that script's content into this one
F1 <- readRDS(file.path("output", "Frechet_1.rds"))
F2 <- readRDS(file.path("output", "Frechet_2.rds"))
F3 <- readRDS(file.path("output", "Frechet_3.rds"))

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
plot_df$label <- factor(plot_df$label, levels = c("true", "permuted", "random"))
levels(plot_df$label) <- c("Observed", "Permutated MAP", "Random dynamics")
p <- ggplot(plot_df,
            aes(x = x, fill = label)) +
  geom_density(alpha = 0.5) +
  theme_bw() +
  labs(x = "distance from mean dynamics",
       fill = "Data source") +
  xlim(c(0, 1000))

ggsave(file.path("output", "figures", "S8.svg"),
       dpi = 100,
       units = "in",
       height = 4.5,
       width = 7)

v1 <- sum(D1) / (N-1)
v2 <- sum(D2) / (N-1)
v3 <- sum(D3) / (N-1)

cat(paste0("Ratio of variance observed to permuted: ", round(v1 / v2, 3), ":1\n"))
cat(paste0("Ratio of variance observed to random: ", round(v1 / v3, 3), ":1\n"))
