library(matrixsampling)
library(Rcpp)
library(brms)
library(driver)
library(fido)
library(rulesoflife)

H <- 56
T <- 100
D <- 135

global_dynamics <- rinvwishart(1, D + 2, diag(D))[,,1]
global_baseline <- rdirichlet(1, rep(20/D, D))

mixing_prop <- 0.5

simulated <- NULL
for(h in 1:H) {
  host_X <- round(matrix(runif(T, min = 1, max = T*10), 1, T))
  host_dynamics <- rinvwishart(1, D + 2, diag(D))[,,1]
  combined_dynamics <- mixing_prop*global_dynamics + (1 - mixing_prop)*host_dynamics

  # Convert to ALR; this is a fido function
  combined_alr <- clrvar2alrvar(combined_dynamics, D)

  # Simulate from model using these combined dynamics
  rho <- calc_se_decay(min_correlation = 0.1,
                       days_to_min_autocorrelation = 90)
  Gamma <- function(X) {
    jitter <- 1e-08
    SE(X[1,,drop=F], sigma = 1, rho = rho, jitter = jitter)
  }

  # Sample a random composition, transform to ALR
  Theta <- function(X) matrix(alr(global_baseline), D-1, ncol(X))

  Lambda <- rmatrixnormal(1, Theta(host_X), combined_alr, Gamma(host_X))[,,1]
  Eta <- rmatrixnormal(1, Lambda, combined_alr, diag(T))[,,1]
  proportions <- alrInv_array(Eta, d = D, coords = 1)
  # Eta.clr <- clr_array(proportions, parts = 1)
  Y <- apply(proportions, 2, function(p) {
    rmultinom(1, p, size = rpois(1, 10000))
  })

  prior_cov_taxa <- get_Xi(D, total_variance = 1)

  # Fit model; for some reason this is crashing for me locally (Windows)
  # Need to check fido install -- LEFT OFF HERE!
  fit <- fido::basset(Y = Y, X = host_X, upsilon = prior_cov_taxa$upsilon,
                      Xi = prior_cov_taxa$Xi, Theta, Gamma,
                      n_samples = 0, ret_mean = TRUE)
  print("Model fitted")

  simulations[[h]] <- list(Sigma = combined_dynamics,
                           Y = Y,
                           fit = fit)
}

# Build the "rug" estimated from these H samples

# Calculate the universality scores

# Use some metric (KL divergence?) to evaluate which mixing proportion gives a
# simulated distribution closest to our observed one (reference)

sim <- pibble_sim()
fit <- pibble(sim$Y, sim$X)






