library(fido)
library(driver)
library(rulesoflife)
library(dplyr)
library(ggplot2)
library(cowplot)

# Code from: https://github.com/yuanpeicao/COAT
# Paper: https://arxiv.org/pdf/1601.04397.pdf
source("COAT.R")

data <- load_data()
hosts <- unique(data$metadata$sname)
rug_basset <- matrix(NA, length(hosts), 8911)
rug_coat <- matrix(NA, length(hosts), 8911)
for(h in 1:length(hosts)) {
  sname <- hosts[h]
  cat(paste0("Processing host: ", sname, "\n"))
  # Pull cached basset parameter estimates
  basset_estimate <- readRDS(file.path("output",
                                       "model_fits",
                                       "asv_days90_diet25_scale1",
                                       "MAP",
                                       paste0(sname, ".rds")))
  # Convert to CLR
  basset_estimate <- to_clr(basset_estimate)
  # Pull covariance matrix
  basset_estimate <- basset_estimate$Sigma[,,1]
  # Scale to correlation
  basset_estimate <- cov2cor(basset_estimate)

  # Pull corresponding raw data for COAT
  data <- load_data(tax_level = "ASV")
  md <- data$metadata
  data <- data$counts
  data <- data[,which(md$sname == sname)]

  # Get estimate from COAT
  coat_estimate <- coat(t(data + 0.5), soft = 1)$corr

  plot_kernel_or_cov_matrix(basset_estimate)
  plot_kernel_or_cov_matrix(coat_estimate)

  # Exclude "other"
  basset_estimate <- basset_estimate[1:134,1:134]
  coat_estimate <- coat_estimate[1:134,1:134]

  x <- basset_estimate[upper.tri(basset_estimate, diag = FALSE)]
  y <- coat_estimate[upper.tri(coat_estimate, diag = FALSE)]

  # ggplot(data.frame(x = x, y = y),
  #        aes(x = x, y = y)) +
  #   geom_point(size = 1, shape = 21, fill = "#aaaaaa") +
  #   theme_bw() +
  #   labs(x = "CLR correlation, original model",
  #        y = "CLR correlation, COAT")

  rug_basset[h,] <- x
  rug_coat[h,] <- y
}

# Plot subset of almost 500K points
x <- c(rug_basset)
y <- c(rug_coat)
subset_idx <- sample(1:length(x), size = 1e5)
p1 <- ggplot(data.frame(x = x[subset_idx], y = y[subset_idx]),
       aes(x = x, y = y)) +
  # geom_point(size = 1, shape = 21, fill = "#aaaaaa") +
  geom_point(size = 1, alpha = 0.1) +
  theme_bw() +
  labs(y = "CLR correlation, COAT",
       x = "CLR correlation, basset")

p2 <- ggplot(data.frame(x = c(x, y),
                        Model = c(rep("basset", length(x)),
                                 rep("COAT", length(y)))),
             aes(x = x, fill = Model)) +
  geom_density(alpha = 0.5) +
  theme_bw() +
  scale_fill_manual(values = c("#59AAD7", "#aaaaaa")) +
  labs(x = "estimated taxon-taxon correlation")

p <- plot_grid(p1, p2, ncol = 2, rel_widths = c(1, 1.2))

ggsave(file.path("output", "figures", "basset_v_COAT.png"),
       plot = p,
       dpi = 100,
       units = "in",
       height = 4,
       width = 9)

# Calculate correlation
cat(paste0("Correlation of basset, COAT estimates: ",
           round(cor(x, y), 3),
           "\n"))
