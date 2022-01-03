source("path_fix.R")

library(tidyverse)
library(rulesoflife)
library(shapes)
library(frechet)
library(matrixsampling)
library(cowplot)

# ------------------------------------------------------------------------------
#
#   Supplemental Figure S6 - distributions of distances from mean dynamics for
#                            dynamics estimates
#
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
#   Distributions of distances
# ------------------------------------------------------------------------------

# Pull the already-parsed MAP estimates of dynamics from this object, calculated
# by `analysis_compute_Frechets.R`
F1 <- readRDS(file.path("output", "Frechet_1_corr.rds"))

# Calculate the mean using the `frechet` package
F1$mean <- CovFMean(F1$Sigmas)$Mout[[1]]

D <- dim(F1$Sigmas)[1]
N <- dim(F1$Sigmas)[3]

# Compute distances
D1 <- numeric(N)
for(i in 1:N) {
  D1[i] <- dist4cov(F1$Sigmas[,,i] + diag(D)*1e-06,
                    F1$mean + diag(D)*1e-06)$dist
}

# Densities
p1 <- ggplot(data.frame(x = D1), aes(x = x)) +
  geom_histogram(color = "white", binwidth = 0.8) +
  theme_bw() +
  labs(x = "distance from mean dynamics") +
  theme(legend.position = "none")

# ggsave(file.path("output", "figures", "S6a.png"),
#        dpi = 100,
#        units = "in",
#        height = 4,
#        width = 4)

# bounds <- quantile(D1, probs = c(0.025, 0.975))

# ------------------------------------------------------------------------------
#
#   Simulation results (VERSION 1): mapping distances to proportions of global-
#   and host-level signals
#
# ------------------------------------------------------------------------------

D <- 135
# A <- matrixsampling::rinvwishart(1, D + 2, diag(D))[,,1]
A <- F1$mean
addend <- diag(D)*1e-06
mix <- seq(from = 0, to = 1, by = 0.05)

plot_df <- NULL
for(j in 1:1000) {
  cat(paste0("Iteration ", j, "\n"))
  # Include a random host contribution
  B <- cov2cor(matrixsampling::rinvwishart(1, D + 2, diag(D))[,,1])
  # Alternatively, use a random host contribution with strictly the same scale
  # as our observed data (i.e. from permuted estimates)
  # B <- F2$Sigmas[,,sample(1:56, size = 1)]

  dists <- numeric(length(mix))
  for(i in 1:length(mix)) {
    dists[i] <- dist4cov(A + addend,
                         ((1-mix[i])*A + mix[i]*B) + addend)$dist
  }

  plot_df <- rbind(plot_df,
                   data.frame(x = mix, y = dists, iteration = j))
}

p2 <- ggplot(plot_df, aes(x = x, y = y, group = x)) +
  geom_boxplot(outlier.shape = NA) +
  theme_bw() +
  labs(x = "proportion host signal",
       y = "distance") +
  theme(legend.position = "none")

min_obs <- min(D1)
max_obs <- max(D1)
plot_df2 <- data.frame(x = plot_df %>%
                         filter(y > min_obs & y < max_obs) %>%
                         pull(x))
p3 <- ggplot(plot_df2,
       aes(x = x)) +
  geom_histogram(color = "white", binwidth = 0.05) +
  theme_bw() +
  labs(x = "proportion host signal")

p <- plot_grid(p1, p2, p3, ncol = 3,
          scale = 0.95,
          labels = c("a", "b", "c"),
          label_size = 20,
          label_x = -0.02,
          label_y = 1.02)

ggsave(file.path("output", "figures", "S6.png"),
       p,
       dpi = 100,
       units = "in",
       height = 4,
       width = 12)

# ------------------------------------------------------------------------------
#
#   Simulation results (VERSION 2): add in some proportion of the residual until
#   we get data that "looks like" the observed data
#
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
#
#   Cartoon explanatory figure
#
# ------------------------------------------------------------------------------

y <- matrix(c(1, 0.8, -0.5,
              0.8, 1, -0.3,
              -0.5, -0.3, 1), 3, 3, byrow = TRUE)
m <- matrix(c(1, 0.3, -0.4,
              0.3, 1, 0.1,
              -0.4, 0.1, 1), 3, 3, byrow = TRUE)
e <- y - m
diag(e) <- 1

y_plot <- plot_kernel_or_cov_matrix(y) +
  labs(title = "observed host dynamics\n(y)") +
  geom_text(aes(label = round(covariance, 1))) +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        plot.title = element_text(size = 10, hjust = 0.5)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))
m_plot <- plot_kernel_or_cov_matrix(m) +
  labs(title = "population mean dynamics\n(m)") +
  geom_text(aes(label = round(covariance, 1))) +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        plot.title = element_text(size = 10, hjust = 0.5)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))
e_plot <- plot_kernel_or_cov_matrix(e) +
  geom_text(aes(label = round(covariance, 1))) +
  labs(title = "residual host dynamics\n(e)") +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        plot.title = element_text(size = 10, hjust = 0.5)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))
c_plot <- plot_kernel_or_cov_matrix(0.5*m + 0.5*e) +
  geom_text(aes(label = round(covariance, 1))) +
  labs(title = "composite dynamics\n(0.5 m + 0.5 e)") +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        plot.title = element_text(size = 10, hjust = 0.5)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))

p <- plot_grid(y_plot, m_plot, e_plot, c_plot,
               ncol = 4)

ggsave(file.path("output", "figures", "S6_cartoon.png"),
       p,
       dpi = 100,
       units = "in",
       height = 2,
       width = 7)

# Find the distance-minimizer for this tiny example
# props <- seq(0, 1, by = 0.1)
# dists <- c()
# for(prop in props) {
#   dists <- c(dists,
#              dist4cov(y, (1-prop)*m + prop*e)$dist)
# }
# plot(props, dists)

global_mean <- F1$mean
mix <- seq(from = 0, to = 1, by = 0.05)
mins <- NULL
p1 <- NULL
legend <- NULL
for(h in 1:N) {
  host_obs <- F1$Sigmas[,,h]
  host_residual <- host_obs - global_mean
  diag(host_residual) <- 1
  if(is.null(p1)) {
    p1a <- plot_kernel_or_cov_matrix(host_obs) +
      theme(legend.position = "bottom") +
      labs(title = "observed host dynamics",
           fill = "correlation ") +
      guides(fill = guide_colourbar(title.vjust = 0.8))
    legend <- get_legend(p1a)
    p1a <- p1a +
      theme(legend.position = "none",
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            plot.title = element_text(size = 10, hjust = 0.5)) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0))
    p1b <- plot_kernel_or_cov_matrix(global_mean) +
      labs(title = "population mean dynamics") +
      theme(legend.position = "none",
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            plot.title = element_text(size = 10, hjust = 0.5)) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0))
    p1c <- plot_kernel_or_cov_matrix(host_residual) +
      labs(title = "host residual dynamics") +
      theme(legend.position = "none",
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            plot.title = element_text(size = 10, hjust = 0.5)) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0))
    p1 <- plot_grid(p1a, p1b, p1c, ncol = 3)
  }

  for(j in 1:length(mix)) {
    combined_dynamics <- (1-mix[j])*global_mean + mix[j]*host_residual
    mins <- rbind(mins,
                  data.frame(host = h,
                             p = mix[j],
                             dist = dist4cov(host_obs, combined_dynamics)$dist))
  }
}

p2 <- ggplot(mins, aes(x = p, y = dist, color = factor(host))) +
  geom_line() +
  theme_bw() +
  scale_color_manual(values = generate_highcontrast_palette(N)) +
  theme(legend.position = "none") +
  labs(x = "proportion host-level effect",
       y = "distance of composite\ndynamics from observed")

p3 <- plot_grid(p1, NULL, legend, ncol = 1, rel_heights = c(1, -0.1, 0.7))
p4 <- plot_grid(p3, NULL, ncol = 1, rel_heights = c(1, 0.1))
p_all <- plot_grid(p4, p2, ncol = 2,
                   rel_widths = c(1, 0.75),
                   labels = c("a", "b"),
                   label_size = 18,
                   scale = 0.90)

ggsave(file.path("output", "figures", "S6.png"),
       p_all,
       dpi = 100,
       units = "in",
       height = 3.5,
       width = 10)

