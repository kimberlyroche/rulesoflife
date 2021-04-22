#' Simulate count data from the generative model, downsampling it repeatedly,
#' and estimate the CLR correlation across taxa. Score the fidelity of estimates
#' on this downsampled data by calculating the proportion of taxon-taxon
#' correlations matched for sign.
#'
#' @param sample_no vector of sample counts; simulated data will be iteratively
#' downsampled to these sample counts
#' @import driver
#' @import matrixsampling
#' @import fido
#' @import phyloseq
#' @export
downsample_counts <- function(sample_no = c(10, 20, 30, 40, 50, 75, 100, 200)) {
  results <- data.frame(sample_no = sample_no)
  results$matched_CLR_sign = rep(NA, nrow(results))

  # Simulate data
  D <- 100 # number of taxa
  N <- 1000 # samples
  X <- matrix(1:N, 1, N)
  X <- rbind(X, matrix(0, 3, N))
  taxa_prior <- get_Xi(D, total_variance = 1)
  sample_prior <- get_Gamma(kernel_scale = 1,
                            diet_weight = 0,
                            min_correlation = 0.1,
                            days_to_min_autocorrelation = 90)
  Sigma <- rinvwishart(1, taxa_prior$upsilon, taxa_prior$Xi)[,,1]
  # Convert to CLR correlation matrix for later
  Sigma.clr <- alrvar2clrvar(Sigma, d1 = D)
  Sigma.clr.cor <- cov2cor(Sigma.clr)

  # This (logratio) variance gives ~30-50% zero counts
  Theta <- matrix(rnorm(D-1, 0, 3), D-1, N)
  Lambda <- rmatrixnormal(1, Theta, Sigma, sample_prior$Gamma(X))[,,1]
  Eta <- rmatrixnormal(1, Lambda, Sigma, diag(N))[,,1]
  Pi <- alrInv_array(Eta, coords = 1)
  Y <- matrix(NA, D, N)
  for(i in 1:N) {
    Y[,i] <- rmultinom(1, prob = Pi[,i], size = rpois(1, 10000))
  }

  # Downsample, fit model, and compare results
  for(i in 1:nrow(results)) {
    s_no <- results$sample_no[i]
    ds_idx <- sort(sample(1:N, size = s_no))
    obs_Y <- Y[,ds_idx]
    obs_N <- length(ds_idx)
    obs_X <- matrix(ds_idx, 1, obs_N)
    obs_X <- rbind(obs_X,
                   matrix(0, 3, obs_N))
    alr_Y <- alr_array(obs_Y + 0.5, parts = 1)
    alr_means <- rowMeans(alr_Y)
    Theta_prior <- function(X) matrix(alr_means, D-1, ncol(obs_X))

    fit <- fido::basset(Y = obs_Y,
                        X = obs_X,
                        upsilon = taxa_prior$upsilon,
                        Xi = taxa_prior$Xi,
                        Theta = Theta_prior,
                        Gamma = sample_prior$Gamma,
                        n_samples = 0,
                        ret_mean = TRUE)

    # Get model estimate
    Sigma_hat <- fit$Sigma[,,1]
    # Convert to CLR

    Sigma_hat.clr <- alrvar2clrvar(Sigma_hat, d1 = D)
    Sigma_hat.clr.cor <- cov2cor(Sigma_hat.clr)

    sign1 <- sign(Sigma.clr.cor[upper.tri(Sigma.clr.cor, diag = FALSE)])
    sign2 <- sign(Sigma_hat.clr.cor[upper.tri(Sigma_hat.clr.cor, diag = FALSE)])

    results$matched_CLR_sign[i] <- sum(sign1 == sign2) / length(sign1)
  }
  return(results)
}
