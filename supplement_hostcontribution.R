source("path_fix.R")

library(tidyverse)
library(rulesoflife)
library(shapes)
library(frechet)
library(matrixsampling)
library(cowplot)

# ------------------------------------------------------------------------------
#
#   (A) Cartoon explanatory figure
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
  labs(title = "observed host\ndynamics\n(y)") +
  geom_text(aes(label = round(covariance, 2))) +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        plot.title = element_text(size = 12, hjust = 0.5)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))
m_plot <- plot_kernel_or_cov_matrix(m) +
  labs(title = "population mean\ndynamics\n(m)") +
  geom_text(aes(label = round(covariance, 2))) +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        plot.title = element_text(size = 12, hjust = 0.5)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))
e_plot <- plot_kernel_or_cov_matrix(e) +
  geom_text(aes(label = round(covariance, 2))) +
  labs(title = "residual host\ndynamics\n(e)") +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        plot.title = element_text(size = 12, hjust = 0.5)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))
c_plot <- plot_kernel_or_cov_matrix(0.5*m + 0.5*e) +
  geom_text(aes(label = round(covariance, 2))) +
  labs(title = "composite\ndynamics\n(0.5 m + 0.5 e)") +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        plot.title = element_text(size = 12, hjust = 0.5)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))

p1 <- plot_grid(y_plot, m_plot, e_plot, c_plot,
                ncol = 4)

# ------------------------------------------------------------------------------
#
#   (B) Distribution of distances for Amboseli data
#
#   Procedure: add in some proportion of the residual until we get data that
#   "looks like" the observed data
#
# ------------------------------------------------------------------------------

# Pull the already-parsed MAP estimates of dynamics from this object, calculated
# by `analysis_compute_Frechets.R`
# F1 <- readRDS(file.path("output", "Frechet_1_corr.rds"))

# Calculate the mean using the `frechet` package
# F1$mean <- CovFMean(F1$Sigmas)$Mout[[1]]
# N <- dim(F1$Sigmas)[3]
# global_mean <- F1$mean

Sigmas <- pull_Sigmas("asv_days90_diet25_scale1")

data <- load_data(tax_level = "ASV")
rug_asv <- summarize_Sigmas(output_dir = "asv_days90_diet25_scale1")
filtered_pairs <- filter_joint_zeros(data$counts, threshold_and = 0.05, threshold_or = 0.5)
filtered_obj <- rug_asv
filtered_obj$tax_idx1 <- filtered_obj$tax_idx1[filtered_pairs$threshold]
filtered_obj$tax_idx2 <- filtered_obj$tax_idx2[filtered_pairs$threshold]
filtered_obj$rug <- filtered_obj$rug[,filtered_pairs$threshold]

# cov_mean <- estcov(Sigmas, method = "Euclidean")
cov_mean <- colMeans(filtered_obj$rug)

N <- dim(Sigmas)[3]
# global_mean <- cov_mean$mean
global_mean <- cov_mean
mix <- seq(from = 0, to = 1, by = 0.05)
mins <- NULL
p1_real <- NULL
legend <- NULL
# p1 <- NULL
for(h in 1:N) {
  # host_obs <- F1$Sigmas[,,h]
  # host_obs <- Sigmas[,,h]
  host_obs <- filtered_obj$rug[h,]
  host_residual <- host_obs - global_mean
  # diag(host_residual) <- 1
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
    p1_real <- plot_grid(p1a, p1b, p1c, ncol = 3)
  }

  for(j in 1:length(mix)) {
    combined_dynamics <- (1-mix[j])*global_mean + mix[j]*host_residual
    mins <- rbind(mins,
                  data.frame(host = h,
                             p = mix[j],
                             # dist = dist4cov(host_obs, combined_dynamics)$dist))
                             dist = as.numeric(dist(matrix(c(host_obs,
                                                             combined_dynamics),
                                                           byrow = T,
                                                           2,
                                                           length(host_obs))))))
  }
}

p2_alt <- ggplot(mins, aes(x = p, y = dist, color = factor(host))) +
  geom_line() +
  theme_bw() +
  scale_color_manual(values = generate_highcontrast_palette(N)) +
  theme(legend.position = "none") +
  labs(x = "proportion host-level effect",
       y = "distance of composite\ndynamics from observed")

p2 <- ggplot(data.frame(x = mins %>%
                          group_by(host) %>%
                          arrange(dist) %>%
                          slice(1) %>%
                          pull(p)), aes(x = 1 - x)) +
  geom_histogram(color = "white", binwidth = 0.05) +
  theme_bw() +
  scale_alpha_continuous(range = c(0.25, 1.0)) +
  coord_cartesian(clip = "off") +
  guides(fill = "none",
         alpha = "none") +
  labs(x = "population mean weight",
       y = "count")

p1_padded <- plot_grid(NULL, p1, NULL, ncol = 1,
                       rel_heights = c(0.1, 1, 0.18))

p2_padded <- plot_grid(p2, scale = 0.9)

p_all <- plot_grid(p1_padded, NULL, p2_padded, ncol = 3,
                   rel_widths = c(1.2, 0.03, 0.5),
                   labels = c("A", "", "B"),
                   label_size = 18,
                   scale = 0.95)

ggsave(file.path("output", "figures", "Figure_3_Supplement_2"),
       p_all,
       dpi = 200,
       units = "in",
       height = 3.5,
       width = 12,
       bg = "white")
