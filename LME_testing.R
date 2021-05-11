# This script plots the "rug" heatmap over hosts and pairwise correlation over
# taxa using various row and column orderings.

library(rulesoflife)
library(tidyverse)
library(coxme)
# library(GMMAT)
library(MASS)

source("path_fix.R")
source("ggplot_fix.R")

# Pull reference data and generate rug object.
# output_dir = "asv_days90_diet25_scale1"
# rug_obj <- summarize_Sigmas(output_dir = output_dir)

n_pairs <- 4
n_hosts <- 56

# Build a covariance matrix with two groups.
# alpha <- 0.5
# K <- matrix(-alpha, n_hosts, n_hosts)
# K[1:(n_hosts/2),1:(n_hosts/2)] <- alpha
# K[(n_hosts/2+1):n_hosts,(n_hosts/2+1):n_hosts] <- alpha
# diag(K) <- 1
# rownames(K) <- 1:nrow(K)
# colnames(K) <- 1:nrow(K)

alpha <- 0.5
K <- matrix(0, n_hosts*n_pairs, n_hosts*n_pairs)
for(i in 1:n_pairs) {
  offset <- (i-1)*n_hosts
  K[(offset+1):(offset+n_hosts),(offset+1):(offset+n_hosts)] <- alpha
}
diag(K) <- 1

plot_kernel_or_cov_matrix(K)

# Draw samples from a variance component model
# host_offsets <- mvrnorm(n = 1, mu = rep(0, n_hosts), Sigma = K)
# response <- mvrnorm(n = n_pairs, mu = host_offsets, Sigma = diag(n_hosts)*0.5)
# response <- t(response) # host x pair
# image(t(response))

# But I think we want to fit something more like...
weight <- 0.99
response <- mvrnorm(n = 1,
                    mu = rep(0, n_hosts*n_pairs),
                    Sigma = weight*K + (1-weight)*diag(n_hosts*n_pairs))
plot(response)

# Pivot data long
response_df <- data.frame(y = response,
                          host = rep(1:n_hosts, n_pairs),
                          pair = rep(1:n_pairs, each = n_hosts))
response_df <- response_df %>%
  mutate(host_pair = paste0(host, "_", pair))
response_df$host <- factor(response_df$host)
response_df$pair <- factor(response_df$pair)
response_df$host_pair <- factor(response_df$host_pair)
head(response_df)

# Fit model with "kinship" matrix
fit <- lmekin(y ~ (1|host_pair), data = response_df, varlist = list(K))
fit

# fit <- glmmkin(y ~ pair + host, data = response_df, kins = K,
#                   id = "host", family = gaussian(link = "identity"))
# fit

# Compute something like variance explained
unname(VarCorr(fit)[[1]]) / var(response_df$y)

# Should be able to use roptim to optimize over the log likelihood and just
# solve for a weight for each component

