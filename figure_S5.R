source("path_fix.R")

library(tidyverse)
library(rulesoflife)
library(shapes)
library(frechet)
library(matrixsampling)
library(ggridges)

# ------------------------------------------------------------------------------
#   Supplemental Figure 5 - distributions of distances from mean dynamics for
#                           observed, permuted, and random correlation matrices
#                           plus simulation results for combined global- and
#                           host-level trends in dynamics
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
#   Distributions of distances
# ------------------------------------------------------------------------------

# These are calculated by `analysis_compute_Frechets.R`
# We ultimately need to roll that script's content into this one
F1 <- readRDS(file.path("output", "Frechet_1_corr.rds"))
F2 <- readRDS(file.path("output", "Frechet_2_corr.rds"))
F3 <- readRDS(file.path("output", "Frechet_3_corr.rds"))

# This function from the frechet package runs faster and gives a different (and
# frankly more intuitive) answer than the shapes package estimates, where means
# had a very very different scale than the samples they were computed from.
# F1$mean <- CovFMean(F1$Sigmas)$Mout[[1]]
# F2$mean <- CovFMean(F2$Sigmas)$Mout[[1]]
# F3$mean <- CovFMean(F3$Sigmas)$Mout[[1]]

D <- dim(F1$Sigmas)[1]
N <- dim(F1$Sigmas)[3]

# From: https://stackoverflow.com/questions/14313285/ggplot2-theme-with-no-axes-or-grid
theme_bare <- function(p, no_x = TRUE, no_y = TRUE) {
  p <- p +
    theme(
      axis.line = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "none",
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing = unit(c(0,0,0,0), "lines"),
      plot.margin = unit(c(0,0,0,0), "lines")
    ) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0))
  if(no_x) {
    p <- p +
      theme(axis.title.x = element_blank())
  }
  if(no_y) {
    p <- p +
      theme(axis.title.y = element_blank())
  }
  p
}

# Plot sample correlation matrices
samples <- sample(1:N, size = 2, replace = FALSE)
plots <- list()
legend <- NULL

# Row 1
plots[[1]] <- NULL
plots[[2]] <- plot_kernel_or_cov_matrix(F1$mean) +
  labs(title = "Mean") +
  theme(axis.title.y = element_text(angle = 0, size = 14))
plots[[2]] <- theme_bare(plots[[2]])
temp <- plot_kernel_or_cov_matrix(F1$Sigmas[,,samples[1]]) +
  labs(fill = "CLR correlation")
legend <- get_legend(temp)
plots[[3]] <- temp +
  labs(title = "Sample 1")
plots[[3]] <- theme_bare(plots[[3]])
plots[[4]] <- plot_kernel_or_cov_matrix(F1$Sigmas[,,samples[2]]) +
  labs(title = "Sample 2")
plots[[4]] <- theme_bare(plots[[4]])
plots[[5]] <- NULL

# Row 2
plots[[6]] <- NULL
plots[[7]] <- plot_kernel_or_cov_matrix(F2$mean)
  theme(axis.title.y = element_text(angle = 0, size = 14))
plots[[7]] <- theme_bare(plots[[7]])
plots[[8]] <- plot_kernel_or_cov_matrix(F2$Sigmas[,,samples[1]])
plots[[8]] <- theme_bare(plots[[8]])
plots[[9]] <- plot_kernel_or_cov_matrix(F2$Sigmas[,,samples[2]])
plots[[9]] <- theme_bare(plots[[9]])
plots[[10]] <- legend

# Row 3
plots[[11]] <- NULL
plots[[12]] <- plot_kernel_or_cov_matrix(F3$mean)
  theme(axis.title.y = element_text(angle = 0, size = 14))
plots[[12]] <- theme_bare(plots[[12]])
plots[[13]] <- plot_kernel_or_cov_matrix(F3$Sigmas[,,samples[1]])
plots[[13]] <- theme_bare(plots[[13]])
plots[[14]] <- plot_kernel_or_cov_matrix(F3$Sigmas[,,samples[2]])
plots[[14]] <- theme_bare(plots[[14]])
plots[[15]] <- NULL

r1 <- plot_grid(plotlist = plots[1:5],
                ncol = 5,
                nrow = 1,
                rel_widths = c(1,1,1,1,1))
r2 <- plot_grid(plotlist = plots[6:10],
                ncol = 5,
                nrow = 1,
                rel_widths = c(1,1,1,1,1))
r3 <- plot_grid(plotlist = plots[11:15],
                ncol = 5,
                nrow = 1,
                rel_widths = c(1,1,1,1,1))

plot_grid(r1, r2, r3, nrow = 3, labels = c("Observed", "Permuted", "  Random"),
          label_size = 13, label_fontface = "plain",
          label_y = 0.5,
          hjust = -0.8)

ggsave(file.path("output", "figures", "sample_matrices.png"),
       dpi = 100,
       units = "in",
       height = 5,
       width = 8)

# Compute squared distances in each case as our variance analog
use_Frobenius <- FALSE
D1 <- numeric(N)
D2 <- numeric(N)
D3 <- numeric(N)
for(i in 1:N) {
  if(use_Frobenius) {
    D1[i] <- dist4cov(F1$Sigmas[,,i] + diag(D)*1e-06,
                      F1$mean + diag(D)*1e-06)$dist**2
    D2[i] <- dist4cov(F2$Sigmas[,,i] + diag(D)*1e-06,
                      F2$mean + diag(D)*1e-06)$dist**2
    D3[i] <- dist4cov(F3$Sigmas[,,i] + diag(D)*1e-06,
                      F3$mean + diag(D)*1e-06)$dist**2
  } else {
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
}

v1 <- sum(D1) / (N-1)
v2 <- sum(D2) / (N-1)
v3 <- sum(D3) / (N-1)

cat(paste0("Ratio of variance observed to permuted: ", round(v1 / v2, 3), ":1\n"))
cat(paste0("Ratio of variance observed to random: ", round(v1 / v3, 3), ":1\n"))

# Plot
plot_df <- data.frame(x = c(D1,
                            D2,
                            D3),
                      label = c(rep("true", N),
                                rep("permuted", N),
                                rep("random", N)))
plot_df$label <- factor(plot_df$label, levels = c("random", "permuted", "true"))
levels(plot_df$label) <- c("Random", "Permuted", "Observed")

p <- ggplot(plot_df,
            aes(x = x, y = label, fill = label)) +
  geom_density_ridges(scale = 1)  +
  theme_bw() +
  labs(x = "distance from mean dynamics",
       y = "",
       fill = "Data source") +
  theme(legend.position = "none") +
  scale_fill_brewer(palette = "RdPu") +
  coord_cartesian(clip = "off") +
  scale_y_discrete(expand = expansion(add = c(0.2, 1.15)))

ggsave(file.path("output", "figures", paste0("S5_", ifelse(use_Frobenius, "Fro", ""), ".png")),
       dpi = 100,
       units = "in",
       height = 4,
       width = 6)

# ------------------------------------------------------------------------------
#   Simulation results
# ------------------------------------------------------------------------------

H <- 56
mixing_props <- seq(from = 0.1, to = 0.9, by = 0.1)
T <- 100
D <- 135

all_scores <- NULL
for(mixing_prop in mixing_props) {
  sim_fn <- file.path("output",
                      "simulations_permuted",
                      paste0("simulations_h-", H, "_m-", mixing_prop, ".rds"))
  simulations <- readRDS(sim_fn)

  # Build the "rug" estimated from these H samples
  D2 <- ((D-1)^2 - (D-1))/2
  rug <- matrix(NA, H, D2)
  for(h in 1:H) {
    Sigma <- cov2cor(simulations[[h]]$Sigma_hat)
    Sigma <- Sigma[1:(D-1),1:(D-1)]
    Sigma <- c(Sigma[upper.tri(Sigma)])
    rug[h,] <- Sigma
  }

  # Calculate the universality scores
  all_scores <- rbind(all_scores,
                      data.frame(score = apply(rug, 2, calc_universality_score),
                                 type = "simulated",
                                 mix = mixing_prop))
}

obs_rug <- readRDS("output/rug_asv.rds")
obs_scores <- apply(obs_rug$rug, 2, calc_universality_score)

all_scores <- rbind(all_scores,
                    data.frame(score = apply(obs_rug$rug, 2, calc_universality_score),
                               type = "observed",
                               mix = NA))

all_scores <- all_scores %>%
  mutate(plot_type = paste0(type, "_", mix))
all_scores$plot_type <- factor(all_scores$plot_type,
                               levels = c("simulated_0.9",
                                          "simulated_0.8",
                                          "simulated_0.7",
                                          "simulated_0.6",
                                          "simulated_0.5",
                                          "simulated_0.4",
                                          "simulated_0.3",
                                          "simulated_0.2",
                                          "simulated_0.1",
                                          "observed_NA"))
levels(all_scores$plot_type) <- c("simulated (90% global pattern)",
                                  "simulated (80% global pattern)",
                                  "simulated (70% global pattern)",
                                  "simulated (60% global pattern)",
                                  "simulated (50% global pattern)",
                                  "simulated (40% global pattern)",
                                  "simulated (30% global pattern)",
                                  "simulated (20% global pattern)",
                                  "simulated (10% global pattern)",
                                  "observed")

p <- ggplot(all_scores, aes(x = score, y = plot_type, fill = plot_type)) +
  geom_density_ridges(scale = 1) +
  # stat_binline(bins = 30, scale = 1, draw_baseline = FALSE) +
  theme_bw() +
  theme(legend.position = "none") +
  scale_fill_brewer(palette = "RdPu") +
  labs(x = "universality score",
       y = "") +
  coord_cartesian(clip = "off") +
  scale_y_discrete(expand = expansion(add = c(0.5, 0.8)))

show(p)

ggsave(file.path("output", "figures", "figure_S8-2.png"),
       p,
       dpi = 100,
       units = "in",
       height = 6,
       width = 5)

# Use some metric (KL divergence?) to evaluate which mixing proportion gives a
# simulated distribution closest to our observed one (reference)
# KLD DOES NOT BEHAVE AS EXPECTED ACROSS RESULTS; TRY SOMETHING ELSE
klds <- c()
for(mixing_prop in mixing_props) {
  klds <- c(klds,
            KLD(all_scores %>% filter(type == "simulated" & mix == mixing_prop) %>% pull(score),
                all_scores %>% filter(type == "observed") %>% pull(score))$sum.KLD.py.px)
}
plot(mixing_props, klds)
