#' Fit MLN Gaussian process model to a given host's series using fido
#'
#' @param sname host short name indicating which baboon's series to fit
#' @param counts filtered 16S count table (taxa x samples)
#' @param metadata annotations data.frame
#' @param point_est flag indicating whether or not to return (single) MAP
#' estimate
#' @return fidofit object
#' @import fido
#' @import driver
#' @export
fit_GP <- function(sname, counts, metadata, point_est = TRUE) {
  sname_idx <- which(metadata$sname == sname)
  sub_md <- metadata[sname_idx,]
  sub_counts <- counts[,sname_idx]

  # Response
  Y <- sub_counts
  N <- ncol(Y)
  D <- nrow(Y)

  # Design matrix
  baseline_date <- sub_md$collection_date[1]
  X <- sapply(sub_md$collection_date, function(cd) {
    round(difftime(cd, baseline_date, units = "days"))
  })
  X <- unname(X + 1)
  dim(X) <- c(1, length(X)) # row matrix

  # Prior/hyperparameters for taxonomic covariance and sample covariance
  var_scale_taxa <- 1
  var_scale_samples <- 1
  K_proportions <- c(0.5, 0.5)
  min_correlation <- 0.1
  days_to_baseline <- 90
  cov_taxa <- get_Xi(D, total_variance = var_scale_taxa)
  cov_sample <- get_Gamma(kernel_scale = var_scale_samples,
                          proportions = K_proportions,
                          min_correlation = min_correlation,
                          days_to_baseline = days_to_baseline)

  # Prior over mean
  alr_ys <- driver::alr((t(Y) + 0.5))
  alr_means <- colMeans(alr_ys)
  Theta <- function(X) matrix(alr_means, D-1, ncol(X))

  if(point_est) {
    n_samples <- 0
    ret_mean <- TRUE
  } else {
    n_samples <- 1000
    ret_mean <- FALSE
  }

  cat(paste0(c(paste0("Fitting fido::basset to ",sname,"'s series with the ",
             "following hyperparams:"),
             paste0("Taxa variance scale: ", var_scale_taxa),
             paste0("Sample variance scale: ", var_scale_samples),
             paste0("Minimum autocorrelation: ", min_correlation),
             paste0("Days to min. autocorrelation: ", days_to_baseline),
             paste0("SE kernel proportion: ", K_proportions[1]),
             paste0("PER kernel proportion: ", K_proportions[2])),
             collapse = "\n\t"))
  fit <- fido::basset(Y = Y, X = X, upsilon = cov_taxa$upsilon, Xi = cov_taxa$Xi,
                      Theta, cov_sample$Gamma, n_samples = n_samples,
                      ret_mean = ret_mean)
  fit$sname <- sname
  if(point_est) {
    output_dir <- file.path("output", "GP_fits", "MAP")
  } else {
    output_dir <- file.path("output", "GP_fits")
  }
  check_dir(output_dir) # create if not exist
  output_file <- paste0(sname, ".rds")
  saveRDS(fit, file = file.path(output_dir, output_file))
  return(fit)
}

#' Fit MLN dynamic linear model to a given host's series using fido
#'
#' @param sname host short name indicating which baboon's series to fit
#' @param counts filtered 16S count table (taxa x samples)
#' @param metadata annotations data.frame
#' @param point_est flag indicating whether or not to return (single) MAP
#' estimate
#' @param smoothed flag indicating whether or not to apply the simulation
#' smoother (necessary for trajectory visualization only)
#' @return fidofit object
#' @import fido
#' @import driver
#' @export
fit_DLM <- function(sname, counts, metadata, point_est = TRUE, smoothed = FALSE) {
  sname_idx <- which(metadata$sname == sname)
  sub_md <- metadata[sname_idx,]
  sub_counts <- counts[,sname_idx]

  # Response
  Y <- sub_counts
  N <- ncol(Y)
  D <- nrow(Y)

  # Design matrix
  # Design matrix
  baseline_date <- sub_md$collection_date[1]
  X <- sapply(sub_md$collection_date, function(cd) {
    round(difftime(cd, baseline_date, units = "days"))
  })
  X <- unname(X + 1)
  T <- max(X)
  dim(X) <- c(1, length(X)) # row matrix

  Q <- 1 # number of covariates: intercept only

  # Transition matrices
  F <- matrix(1, 1, T)
  G <- diag(Q)

  W <- diag(Q)
  # Scale the covariate-inclusive and covariate-exclusive models to have
  # similar total variance
  W <- W/nrow(W)

  # Define the prior over baselines
  C0 <- W
  alr_ys <- driver::alr((t(Y) + 0.5))
  alr_means <- colMeans(alr_ys)
  M0 <- matrix(0, Q, D-1)
  M0[1,] <- alr_means

  # # Prior/hyperparameters for taxonomic covariance and sample covariance
  var_scale <- 1
  cov_taxa <- get_Xi(D, total_variance = 1)

  if(point_est) {
    n_samples <- 0
    ret_mean <- TRUE
  } else {
    n_samples <- 1000
    ret_mean <- FALSE
  }

  cat(paste0(c(paste0("Fitting fido::basset to ",sname,"'s series with the ",
             "following hyperparams:"),
             paste0("Overall variance scale: ", var_scale),
             paste0("Gamma-to-W ratio: 2/3 to 1/3"),
             paste0("Smoothed: ", smoothed)),
             collapse = "\n\t"))

  # I'm giving W about 1/2 the scale of gamma
  fit <- labraduck(Y = Y, observations = X, upsilon = cov_taxa$upsilon, Xi = cov_taxa$Xi,
                   F = F, G = G, W = W, M0 = M0, C0 = C0,
                   gamma_scale = (var_scale * 2/3), W_scale = (var_scale * 1/3),
                   apply_smoother = smoothed, n_samples = n_samples, ret_mean = ret_mean)

  fit$sname <- sname
  if(point_est) {
    output_pieces <- c("output", "DLM_fits", "MAP")
  } else {
    output_pieces <- c("output", "DLM_fits")
  }
  output_dir <- check_dir(output_pieces) # create if not exist
  output_file <- paste0(sname, ".rds")
  saveRDS(fit, file = file.path(output_dir, output_file))
  return(fit)
}

#' Set up a basic ALR prior
#'
#' @param D number of features including reference (where the ALR will represent
#'  D-1)
#' @param total_variance scale of the log variance
#' @return list containing inverse Wishart parameters degrees of freedom and
#' scale matrix
#' @export
get_Xi <- function(D, total_variance = 1) {
  upsilon <- D - 1 + 10 # specify low certainty/concentration
  GG <- cbind(diag(D-1), -1) # log contrast for ALR with last taxon as reference
  Xi <- GG%*%(diag(D)*total_variance)%*%t(GG) # take diag as covariance over log
                                              # abundances
  Xi <- Xi*(upsilon-D-1)
  return(list(upsilon = upsilon, Xi = Xi))
}

#' Periodic kernel
#'
#' @param X covariate (dimension Q x N; i.e., covariates x samples)
#' @param sigma scalar parameter
#' @param rho scalar bandwidth parameter
#' @param period period length (in days)
#' @param jitter small scalar to add to off-diagonal of gram matrix
#'   (for numerical underflow issues)
#' @return Gram Matrix (N x N) (e.g., the Kernel evaluated at
#' each pair of points)
#' @export
PER <- function(X, sigma=1, rho=1, period=24, jitter=0){
  dist <- as.matrix(dist(t(X)))
  G <- sigma^2 * exp(-2*(sin(pi*dist/period)^2)/(rho^2)) +
    jitter*diag(ncol(dist))
}

#' Define bandwidth of squared exponential kernel
#'
#' @param min_correlation minimum correlation to assume between (within-host)
#' samples
#' @param days_to_baseline days at which squared exponential kernel decays to
#' baseline correlation of ~0.1
#' @return bandwidth parameter for squared exponential kernel
#' @export
calc_se_decay <- function(min_correlation = 0.1, days_to_baseline = 90) {
  # Back-calculate the squared exponential bandwidth parameter by finding a
  # bandwidth that gives
  # A desired minimum correlation at the number of days specified by
  # days_to_baseline
  rho <- sqrt(-days_to_baseline^2/(2*log(min_correlation))) # back calculate the
                                                            # decay
  return(rho)
}

#' Define a kernel (function) over samples
#'
#' @param kernel_scale total variance for the composite kernel
#' @param proportions proportion variance to attribute to each of 2 kernels (see
#'  details)
#' @param rho bandwidth for SE kernel
#' @details Composite kernel is built from (1) squared exponential kernel (base
#' autocorrelation component) and (2) seasonal kernel (periodic)
#' @return list containing kernel function and bandwidth parameter
#' @import fido
#' @export
get_Gamma <- function(kernel_scale, proportions, min_correlation = 0.1,
                      days_to_baseline = 90) {
  rho <- calc_se_decay(min_correlation = min_correlation,
                       days_to_baseline = days_to_baseline)
  # Back-calculate the squared exponential bandwidth parameter by finding a
  # bandwidth that gives
  # A desired minimum correlation at the number of days specified by
  # SE_days_to_baseline
  Gamma <- function(X) {
    jitter <- 1e-08
    proportions <- abs(proportions)/sum(abs(proportions)) # just in case we pass
                                                          # something that isn't
                                                          # a composition
    part.1 <- kernel_scale * proportions[1]
    part.2 <- kernel_scale * proportions[2]
    SE(X[1,,drop=F], sigma = sqrt(part.1), rho = rho, jitter = jitter) +
      PER(X[1,,drop=F], sigma = sqrt(part.2), rho = 1, period = 365, jitter = jitter)
  }
  return(list(rho = rho, Gamma = Gamma))
}
