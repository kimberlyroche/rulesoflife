source("path_fix.R")

library(rulesoflife)
library(fido)
library(shapes)
library(tidyverse)
library(optparse)

option_list = list(
  make_option(c("--corr"), type = "logical", default = FALSE,
              help = "convert covariance to correlation", metavar = "logical"),
  make_option(c("--clr"), type = "logical", default = FALSE,
              help = "convert to CLR", metavar = "logical")
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

use_corr <- opt$corr
use_clr <- opt$clr

output_dir <- "asv_days90_diet25_scale1"
output_dir_full <- check_dir(c("output", "model_fits", output_dir, "MAP"))

# Load a few posterior samples for a few hosts
data <- load_data(tax_level = "ASV")
md <- data$metadata
hosts <- sort(unique(md$sname))

# Pull feature number
fit <- readRDS(file.path(output_dir_full, paste0(hosts[1], ".rds")))
D <- fit$D

if(use_clr) {
  S <- array(NA, dim = c(D, D, length(hosts)))
} else {
  S <- array(NA, dim = c(D-1, D-1, length(hosts)))
}

for(i in 1:length(hosts)) {
  host <- hosts[i]
  cat(paste0("Parsing host ", host, "...\n"))
  fit <- readRDS(file.path(output_dir_full, paste0(host, ".rds")))
  if(use_clr) {
    fit <- to_clr(fit)
  }
  Sigma <- fit$Sigma[,,1] + diag(nrow(fit$Sigma[,,1]))*1e-06
  if(use_corr) {
    Sigma <- cov2cor(Sigma)
  }
  S[,,i] <- Sigma
}

Frechet_mean <- readRDS(file.path("output", paste0("Frechet_corr",
                                                   ifelse(use_corr, 1, 0),
                                                   "_",
                                                   ifelse(use_clr, "clr", "alr"),
                                                   ".rds")))

# How similar does a random sample look to the Frechet mean
# plot_kernel_or_cov_matrix(Frechet_mean$mean)
# plot_kernel_or_cov_matrix(S[,,sample(1:length(hosts), size = 1)])

# Sweep sigma scales and calculate the new distance
scales <- seq(from = 0.05, to = 0.25, length.out = 5)
distances <- data.frame(distance = c(), sigma = c(), type = c())
for(i in 1:length(scales)) {
  sigma <- scales[i]
  cat(paste0("Evaluating sigma = ", round(sigma, 2), "\n"))
  for(j in 1:20) {
    cf <- matrix(rnorm(dim(S)[1]**2, 0, sigma), dim(S)[1], dim(S)[1])
    N <- cf%*%t(cf) + diag(dim(S)[1])*1e-06
    distances <- rbind(distances,
                       data.frame(distance = distcov(Frechet_mean$mean,
                                                     Frechet_mean$mean + N,
                                                     method = "Riemannian"),
                                  sigma = sigma,
                                  type = paste0("sigma = ", round(sigma, 2))))
  }
}

# Calculate distances between host samples and this mean
# Sweep sigma scales and calculate the new distance

for(i in 1:length(hosts)) {
  distances <- rbind(distances,
                     data.frame(distance = distcov(Frechet_mean$mean,
                                                   S[,,i],
                                                   method = "Riemannian"),
                                sigma = sigma,
                                type = "host"))
}

ggplot(distances, aes(x = distance, fill = factor(type))) +
  geom_density(color = "white", alpha = 0.5)

# This seems not right. Squaring a noise matrix just produces mostly zero noise
# with a few big perturbations. Most of the increase in distance in an increase
# in the scale of the diagonal, which isn't really what we want.
#
# That would mean that overall magnitude of variation is was drives differences.
# We're interested in a combination of that and scatter.


