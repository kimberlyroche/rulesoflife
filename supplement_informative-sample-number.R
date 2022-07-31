source("path_fix.R")

library(matrixsampling)
library(Rcpp)
library(driver)
library(fido)
library(rulesoflife)
library(LaplacesDemon)

save_fn <- file.path("output", "downsampling_experiment.rds")
if(!file.exists(save_fn)) {
  T <- 1000
  D <- 126

  prop_agree <- NULL
  for(i in 1:1) {
    global_baseline <- LaplacesDemon::rdirichlet(1, rep(100/D, D))
    taxon_dynamics <- cov2cor(matrixsampling::rinvwishart(1, D + 2, diag(D))[,,1])
    taxon_dynamics_alr <- clrvar2alrvar(taxon_dynamics, D)

    X <- matrix(1:T, 1, T)

    # Simulate from model using these combined dynamics
    rho <- calc_se_decay(min_correlation = 0.1,
                         days_to_min_autocorrelation = 90)
    Gamma <- function(X) {
      jitter <- 1e-08
      SE(X[1,,drop=F], sigma = 1, rho = rho, jitter = jitter)
    }

    # Sample a random composition, transform to ALR
    Theta <- function(X) matrix(alr(global_baseline), D-1, ncol(X))

    Lambda <- rmatrixnormal(1, Theta(X), taxon_dynamics_alr, Gamma(X))[,,1]
    Eta <- rmatrixnormal(1, Lambda, taxon_dynamics_alr, diag(T))[,,1]
    proportions <- alrInv_array(Eta, d = D, coords = 1)
    # Eta.clr <- clr_array(proportions, parts = 1)
    Y <- apply(proportions, 2, function(p) {
      rmultinom(1, p, size = rpois(1, 10000))
    })

    for(j in c(10, 20, 30, 40, 50, 75, 100, 200)) {
      cat(paste0("Downsapling to ", j, " samples\n"))
      # Downsample
      X2 <- matrix(sort(sample(1:ncol(Y), size = j)), 1, j)
      Y2 <- Y[,c(X2)]

      prior_cov_taxa <- get_Xi(D, total_variance = 1)

      fit <- fido::basset(Y = Y2, X = X2, upsilon = prior_cov_taxa$upsilon,
                          Xi = prior_cov_taxa$Xi, Theta, Gamma,
                          n_samples = 0, ret_mean = TRUE)

      fit.clr <- to_clr(fit)
      est_Sigma <- cov2cor(fit.clr$Sigma[,,1])

      Z1 <- sign(taxon_dynamics[upper.tri(taxon_dynamics, diag = FALSE)])
      Z2 <- sign(est_Sigma[upper.tri(est_Sigma, diag = FALSE)])

      prop_agree <- rbind(NULL,
                          data.frame(prop = sum(Z1 == Z2) / length(Z1),
                                     n = T,
                                     n_subset = j))
    }
  }
  saveRDS(prop_agree, save_fn)
}

prop_agree <- readRDS(save_fn)

p <- ggplot(prop_agree, aes(x = factor(n_subset), y = prop)) +
  geom_boxplot() +
  theme_bw() +
  labs(x = "sample number",
       y = "proportion matched CLR correlation sign")

ggsave(file.path("output", "figures", "informative_sample_no.svg"),
       plot = p,
       dpi = 100,
       units = "in",
       height = 4,
       width = 5)
