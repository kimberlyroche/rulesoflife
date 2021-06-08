#' Simulate series for one host
#'
#' @param D number of features
#' @param N number of samples (time points)
#' @param Sigma optional correlation matrix to simulate from
#' @param Theta optional baseline for logratios
#' @import matrixsampling
#' @import driver
#' @export
simulate_series <- function(D = 100, N = 1000, Sigma = NULL, Theta = NULL) {
  # Simulate data
  X <- matrix(1:N, 1, N)
  X <- rbind(X, matrix(0, 3, N))
  taxa_prior <- get_Xi(D, total_variance = 1)
  sample_prior <- get_Gamma(kernel_scale = 1,
                            diet_weight = 0,
                            min_correlation = 0.1,
                            days_to_min_autocorrelation = 90)
  if(is.null(Sigma)) {
    Sigma <- rinvwishart(1, taxa_prior$upsilon, taxa_prior$Xi)[,,1]
  }

  # This (logratio) variance gives ~30-50% zero counts
  if(is.null(Theta)) {
    Theta <- matrix(rnorm(D-1, 0, 3), D-1, N)
  }
  Lambda <- rmatrixnormal(1, Theta, Sigma, sample_prior$Gamma(X))[,,1]
  Eta <- rmatrixnormal(1, Lambda, Sigma, diag(N))[,,1]
  Pi <- alrInv_array(Eta, coords = 1)
  Y <- matrix(NA, D, N)
  for(i in 1:N) {
    Y[,i] <- rmultinom(1, prob = Pi[,i], size = rpois(1, 10000))
  }
  return(list(D = D, N = N, X = X, Sigma = Sigma, Gamma = sample_prior$Gamma,
              upsilon = taxa_prior$upsilon, Xi = taxa_prior$Xi,
              Theta = Theta, Lambda = Lambda, Eta = Eta, Pi = Pi, Y = Y))
}

#' Simulate count data from the generative model, downsampling it repeatedly,
#' and estimate the CLR correlation across taxa. Score the fidelity of estimates
#' on this downsampled data by calculating the proportion of taxon-taxon
#' correlations matched for sign.
#'
#' @param sample_no vector of sample counts; simulated data will be iteratively
#' downsampled to these sample counts
#' @import driver
#' @import fido
#' @import phyloseq
#' @export
downsample_counts <- function(sample_no = c(10, 20, 30, 40, 50, 75, 100, 200)) {
  results <- data.frame(sample_no = sample_no)
  results$matched_CLR_sign = rep(NA, nrow(results))

  sim_data <- simulate_series()
  D <- sim_data$D
  N <- sim_data$N
  Y <- sim_data$Y
  upsilon <- sim_data$upsilon
  Xi <- sim_data$Xi
  Gamma <- sim_data$Gamma
  Sigma <- sim_data$Sigma

  # Convert to CLR correlation matrix for later
  Sigma_obj <- convert_alr_Sigma_clr(Sigma)
  Sigma.clr <- Sigma_obj$Sigma.clr
  Sigma.clr.cor <- Sigma_obj$Sigma.clr.cor

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
    Theta_prior <- function(X) matrix(alr_means, D-1, ncol(X))

    fit <- fido::basset(Y = obs_Y,
                        X = obs_X,
                        upsilon = upsilon,
                        Xi = Xi,
                        Theta = Theta_prior,
                        Gamma = Gamma,
                        n_samples = 0,
                        ret_mean = TRUE)

    # Get model estimate
    Sigma_hat <- fit$Sigma[,,1]
    # Convert to CLR
    Sigma_hat_obj <- convert_alr_Sigma_clr(Sigma_hat)
    Sigma_hat.clr <- Sigma_hat_obj$Sigma.clr
    Sigma_hat.clr.cor <- Sigma_hat_obj$Sigma.clr.cor

    sign1 <- sign(Sigma.clr.cor[upper.tri(Sigma.clr.cor, diag = FALSE)])
    sign2 <- sign(Sigma_hat.clr.cor[upper.tri(Sigma_hat.clr.cor, diag = FALSE)])

    results$matched_CLR_sign[i] <- sum(sign1 == sign2) / length(sign1)
  }
  return(results)
}

#' Boxplot data simulated from downsample_counts()
#'
#' @param data data.frame containing results from downsample_counts()
#' @import ggplot2
#' @export
plot_downsampling <- function(data) {
  p <- ggplot(data, aes(x = factor(sample_no), y = matched_CLR_sign)) +
    geom_boxplot() +
    xlab("sample number") +
    ylab("% CLR sign match (true, estimated)")
  output_dir <- check_dir(c("output", "figures"))
  ggsave(file.path(output_dir, "downsampling.png"),
         p,
         units = "in",
         height = 3,
         width = 4)
}
