source("path_fix.R")

library(rulesoflife)
library(tidyverse)
library(coxme)
library(MASS)
library(matrixsampling)

source("ggplot_fix.R")

# ------------------------------------------------------------------------------
#   Simulate data to make sure we're thinking of the right model
# ------------------------------------------------------------------------------

n_pairs <- 10
n_hosts <- 6

rho <- 0.8
K <- matrix(-rho, n_hosts, n_hosts)
K[1:(n_hosts/2),1:(n_hosts/2)] <- rho
K[(n_hosts/2+1):n_hosts,(n_hosts/2+1):n_hosts] <- rho
diag(K) <- 1
# plot_kernel_or_cov_matrix(K)

mu_mat <- matrix(rnorm(n_pairs), n_hosts, n_pairs, byrow = TRUE)

alpha_true <- 20
use_K <- alpha_true*K #+ (1 - alpha_true)*diag(n_hosts)

# Extreme example of shared deviations!
y <- rmatrixnormal(1, mu_mat, use_K, diag(n_pairs))[,,1]

K_large <- kronecker(diag(n_pairs), use_K)
data <- data.frame(y = c(y),
                   host = rep(1:n_hosts, n_pairs),
                   pair = rep(1:n_pairs, each = n_hosts))
data <- data %>%
  mutate(host_pair = paste0(host, "_", pair))
data$host <- factor(data$host)
data$pair <- factor(data$pair)
data$host_pair <- factor(data$host_pair, levels = data$host_pair)
head(data)

fit <- lmekin(y ~ pair + (1 | host_pair),
              data = data, varlist = list(K_large))
# fit <- lmekin(y ~ pair + (1 | host_pair),
#               data = data, varlist = list(diag(nrow(K_large))))
fit

mu_pred <- matrix(fit$coefficients$fixed, n_hosts, n_pairs, byrow = TRUE)

re_pred <- matrix(NA, n_hosts, n_pairs)
for(re_idx in 1:length(fit$coefficients$random$host_pair)) {
  re <- fit$coefficients$random$host_pair[re_idx]
  pieces <- as.numeric(strsplit(names(re), "_")[[1]])
  re_pred[pieces[1], pieces[2]] <- unname(re)
}

par(mfrow = c(1,2))
image(t(y))
image(t(mu_pred) + t(re_pred))

# Residual variance without the random effect
var_y <- var(c(y))
resid_var <- var(c(y - mu_pred))
resid_var_re <- var(c(y - (mu_pred + re_pred)))

(var_y - resid_var_re)/var_y

# # Get a tiny optimization working.
# likelihood <- function(x) {
#   x1 <- x[[1]]
#   x1**2 + x1*3
# }
# gradient <- function(x) {
#   x1 <- x[[1]]
#   c(2*x1 + 3)
# }
# res <- optim(list(-100), likelihood, gradient, method = "BFGS",
#              control = list(trace=TRUE), hessian = FALSE)
# res$par
# res$value
#
#
# alpha_true <- 5
# mu_true <- rep(0, 10)
# K <- diag(10)
# y <- mvrnorm(1, mu_true, K*alpha_true)
#
# likelihood <- function(x) {
#   alpha <- x[1]
#   A <- y - mu_true
#   dim(A) <- c(length(A), 1)
#   X <- alpha*K
#   0.5*(log(det(X)) + t(A)%*%solve(X)%*%A)
# }
# gradient <- function(x) {
#   alpha <- x[1]
#   A <- y - mu_true
#   dim(A) <- c(length(A), 1)
#   X <- alpha*K
#   X_inv <- solve(X)
#   tr1 <- 1/alpha * length(y)
#   tr2 <- sum(diag(t(A)%*%X_inv%*%K%*%X_inv%*%A))
#   cat("Tr1:",tr1,"\n")
#   cat("Tr2:",tr2,"\n")
#   cat("Tr1 - Tr2:",tr1 - tr2,"\n")
#   0.5*(tr1 - tr2)
# }
# res <- optim(c(0.2), likelihood, gradient, method = "BFGS",
#              control = list(trace = TRUE), hessian = FALSE)
# res$par


# # Pull reference data and generate rug object.
# # output_dir = "asv_days90_diet25_scale1"
# # rug_obj <- summarize_Sigmas(output_dir = output_dir)
#
# n_pairs <- 4
# n_hosts <- 56
#
# # Build a covariance matrix with two groups.
# alpha <- 0.5
# K <- matrix(-alpha, n_hosts, n_hosts)
# K[1:(n_hosts/2),1:(n_hosts/2)] <- alpha
# K[(n_hosts/2+1):n_hosts,(n_hosts/2+1):n_hosts] <- alpha
# diag(K) <- 1
# rownames(K) <- 1:nrow(K)
# colnames(K) <- 1:nrow(K)
#
# # alpha <- 0.5
# # K <- matrix(0, n_hosts*n_pairs, n_hosts*n_pairs)
# # for(i in 1:n_pairs) {
# #   offset <- (i-1)*n_hosts
# #   K[(offset+1):(offset+n_hosts),(offset+1):(offset+n_hosts)] <- alpha
# # }
# # diag(K) <- 1
#
# plot_kernel_or_cov_matrix(K)
#
# # Draw samples from a variance component model
# # host_offsets <- mvrnorm(n = 1, mu = rep(0, n_hosts), Sigma = K)
# # response <- mvrnorm(n = n_pairs, mu = host_offsets, Sigma = diag(n_hosts)*0.5)
# response <- mvrnorm(n = n_pairs, mu = rep(0, n_hosts), Sigma = 0.9*K + 0.1*diag(n_hosts))
# response <- t(response) # host x pair
# image(response)
#
# # But I think we want to fit something more like...
# weight <- 0.99
# response <- mvrnorm(n = 1,
#                     mu = rep(0, n_hosts*n_pairs),
#                     Sigma = weight*K + (1-weight)*diag(n_hosts*n_pairs))
# plot(response)
#
# # Pivot data long
# response_df <- data.frame(y = response,
#                           host = rep(1:n_hosts, n_pairs),
#                           pair = rep(1:n_pairs, each = n_hosts))
# response_df <- response_df %>%
#   mutate(host_pair = paste0(host, "_", pair))
# response_df$host <- factor(response_df$host)
# response_df$pair <- factor(response_df$pair)
# response_df$host_pair <- factor(response_df$host_pair)
# head(response_df)
#
# # Fit model with "kinship" matrix
# fit <- lmekin(y ~ (1|host_pair), data = response_df, varlist = list(K))
# fit
#
# # fit <- glmmkin(y ~ pair + host, data = response_df, kins = K,
# #                   id = "host", family = gaussian(link = "identity"))
# # fit
#
# # Compute something like variance explained
# unname(VarCorr(fit)[[1]]) / var(response_df$y)
#
# # Should be able to use roptim to optimize over the log likelihood and just
# # solve for a weight for each component
#
