library(rulesoflife)
library(fido)
library(shapes)
library(tidyverse)

noise_mat <- function(sigma, D) {
  cf <- matrix(rnorm(D**2, 0, sigma), D, D)
  t(cf)%*%cf
}

output_dir <- "asv_days90_diet25_scale1"
n_iter <- 1

# Load a few posterior samples for a few hosts
data <- load_data(tax_level = "ASV")
md <- data$metadata
hosts <- sort(unique(md$sname))

Sigmas <- NULL
for(host in hosts) {
  cat(paste0("Parsing host ", host, "...\n"))
  fit <- readRDS(file.path("output",
                           "model_fits",
                           output_dir,
                           "full_posterior",
                           paste0(host, ".rds")))
  fit.clr <- to_clr(fit)
  Sigmas[[host]] <- fit.clr$Sigma[,,sample(1:fit$iter, size = n_iter)]
}
D <- dim(Sigmas[[1]])[1]

# Load Frechet mean (need to re-run this on proper MAP estimates)
Frechet_mean <- readRDS(file.path("output", "Frechet_mean.rds"))

S <- noise_mat(1, D)
plot_kernel_or_cov_matrix(S)

# Sweep sigma scales and calculate the new distance
scales <- seq(from = 1, to = 10, length.out = 5)
distances <- data.frame(distance = c(), sigma = c())
for(i in 1:length(scales)) {
  sigma <- scales[i]
  for(j in 1:40) {
    S <- noise_mat(sigma, D)
    S_null <- diag(D)
    diag(S_null) <- diag(S)
    distances <- rbind(distances,
                       data.frame(distance = distcov(Frechet_mean$mean + S_null + diag(D)*1e-06,
                                                     Frechet_mean$mean + S + diag(D)*1e-06,
                                                     method = "Riemannian"),
                                  sigma = sigma))
  }
}
head(distances)

ggplot(distances, aes(x = distance, fill = factor(sigma))) +
  geom_density(color = "white", alpha = 0.5)

# Is this a reasonable kind of noise structure?
# This will get farther from the mean on the basis of the diagonal increasing
# alone, so I'm removing the effect of the diagonal BUT it keeps distances
# really small!


