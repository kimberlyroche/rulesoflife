source("path_fix.R")

library(matrixsampling)
library(Rcpp)
library(driver)
library(fido)
library(rulesoflife)
library(LaplacesDemon)
library(frechet)
library(optparse)

option_list <- list(
  make_option(c("--hosts"), type = "integer", default = 56,
              help = "Number of hosts to simulate",
              metavar = "number"),
  make_option(c("--mix"), type = "numeric", default = 0.5,
              help = "Mixing parameter (global proportion)",
              metavar = "number"),
  make_option(c("--permuted"), type = "logical", default = TRUE,
              help = "Use permuted versions of observed matrices?",
              metavar = "logical")
)

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(option_list = option_list))

H <- opt$hosts
mixing_prop <- opt$mix
permuted <- opt$permuted
T <- 100
D <- 135

if(permuted) {
  samples <- readRDS(file.path("output", "Frechet_1_corr.rds"))
  global_dynamics <- CovFMean(samples$Sigmas)$Mout[[1]]
  shuffle_idx <- sample(1:D)
  global_dynamics <- global_dynamics[shuffle_idx,shuffle_idx]
} else {
  global_dynamics <- matrixsampling::rinvwishart(1, D + 2, diag(D))[,,1]
  global_baseline <- LaplacesDemon::rdirichlet(1, rep(100/D, D))
}

simulations <- NULL
for(h in 1:H) {
  cat(paste0("Simulating host ", h, " / ", H, "\n"))
  host_X <- round(matrix(runif(T, min = 1, max = T*10), 1, T))
  if(permuted) {
    shuffle_idx <- sample(1:D)
    host_dynamics <- samples$Sigmas[shuffle_idx,shuffle_idx,h]
  } else {
    host_dynamics <- matrixsampling::rinvwishart(1, D + 2, diag(D))[,,1]
  }
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

  fit <- fido::basset(Y = Y, X = host_X, upsilon = prior_cov_taxa$upsilon,
                      Xi = prior_cov_taxa$Xi, Theta, Gamma,
                      n_samples = 0, ret_mean = TRUE)

  fit.clr <- to_clr(fit)

  simulations[[h]] <- list(Sigma = combined_dynamics,
                           Y = Y,
                           Sigma_hat = fit.clr$Sigma[,,1])
}
saveRDS(simulations, paste0("simulations_h-", H, "_m-", mixing_prop, ".rds"))
