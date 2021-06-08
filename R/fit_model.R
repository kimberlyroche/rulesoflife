#' Fit MLN Gaussian process model to a given host's series using fido
#'
#' @param sname host short name indicating which baboon's series to fit
#' @param counts filtered 16S count table (taxa x samples)
#' @param metadata annotations data.frame
#' @param output_dir specifies the subdirectory of the model fits directory in
#' which to save the fitted model output
#' @param MAP flag indicating whether or not to return (single) MAP estimate
#' @param days_to_min_autocorrelation used in the calculation of bandwidth for
#' the kernel over samples that models autocorrelation; indicates the number of
#' days at which the squared exponential correlation decays to a minimum (here,
#' 0.1)
#' @param diet_weight relative contribution of first three diet PCs to
#' covariance across samples
#' @param var_scale_taxa scale of the prior associated with the covariance over
#' taxa
#' @param var_scale_samples scale of the hyperparameter associated with the
#' covariance over samples
#' @param use_adam optimize with Adam (occasionally this converges more reliably
#' that default L-BFGS)
#' @param scramble_time optional flag to scramble a host's observeations in
#' time; we'll use this to generate a "null" version of synchrony across hosts
#' @return fidofit object
#' @import fido
#' @import driver
#' @import dplyr
#' @export
fit_GP <- function(sname, counts, metadata, output_dir, MAP = TRUE,
                   days_to_min_autocorrelation = 90, diet_weight = 0,
                   var_scale_taxa = 1, var_scale_samples = 1, use_adam = FALSE,
                   scramble_time = FALSE) {
  if(diet_weight > 1 | diet_weight < 0) {
    stop("Invalid weight assigned to diet components of kernel!")
  }
  sname_idx <- which(metadata$sname == sname)
  sub_md <- metadata[sname_idx,]
  sub_counts <- counts[,sname_idx]

  # Response
  Y <- sub_counts
  N <- ncol(Y)
  D <- nrow(Y)

  if(scramble_time) {
    Y <- Y[,sample(1:ncol(Y))]
  }

  # Design matrix
  baseline_date <- sub_md$collection_date[1]
  X <- sapply(sub_md$collection_date, function(cd) {
    round(difftime(cd, baseline_date, units = "days"))
  })
  X <- unname(X + 1)
  dim(X) <- c(1, length(X)) # row matrix

  # Diet PCs
  # We need to do something with NAs here
  # For now I'm going to interpolate the mean locally in this function
  sub_md$diet_PC1[which(is.na(sub_md$diet_PC1))] <- mean(sub_md$diet_PC1, na.rm = TRUE)
  sub_md$diet_PC2[which(is.na(sub_md$diet_PC2))] <- mean(sub_md$diet_PC2, na.rm = TRUE)
  sub_md$diet_PC3[which(is.na(sub_md$diet_PC3))] <- mean(sub_md$diet_PC3, na.rm = TRUE)

  X <- rbind(X, sub_md$diet_PC1)
  X <- rbind(X, sub_md$diet_PC2)
  X <- rbind(X, sub_md$diet_PC3)

  # Prior/hyperparameters for taxonomic covariance and sample covariance
  min_correlation <- 0.1
  cov_taxa <- get_Xi(D, total_variance = var_scale_taxa)
  cov_sample <- get_Gamma(kernel_scale = var_scale_samples,
                          diet_weight = diet_weight,
                          min_correlation = min_correlation,
                          days_to_min_autocorrelation = days_to_min_autocorrelation)

  # Prior over mean
  alr_ys <- driver::alr((t(Y) + 0.5))
  alr_means <- colMeans(alr_ys)
  Theta <- function(X) matrix(alr_means, D-1, ncol(X))

  if(MAP) {
    n_samples <- 0
    ret_mean <- TRUE
  } else {
    n_samples <- 100
    ret_mean <- FALSE
  }

  cat(paste0(c(paste0("Fitting fido::basset to ",sname,"'s series with the ",
             "following hyperparams:"),
             paste0("Taxa variance scale: ", var_scale_taxa),
             paste0("Sample variance scale: ", var_scale_samples),
             paste0("Minimum autocorrelation: ", min_correlation),
             paste0("Days to min. autocorrelation: ", days_to_min_autocorrelation),
             paste0("Diet kernel proportion: ", diet_weight),
             paste0("Taxa cov scale: ", var_scale_taxa),
             paste0("Sample cov scale: ", var_scale_samples)),
             collapse = "\n\t"), "\n")
  if(use_adam) {
    fit <- fido::basset(Y = Y, X = X, upsilon = cov_taxa$upsilon, Xi = cov_taxa$Xi,
                        Theta, cov_sample$Gamma, n_samples = n_samples,
                        ret_mean = ret_mean, b2 = 0.98, step_size = 0.003,
                        eps_f = 1e-12, eps_g = 1e-05, max_iter = 10000L,
                        optim_method = "adam")
  } else {
    fit <- fido::basset(Y = Y, X = X, upsilon = cov_taxa$upsilon, Xi = cov_taxa$Xi,
                        Theta, cov_sample$Gamma, n_samples = n_samples,
                        ret_mean = ret_mean)
  }
  fit$sname <- sname
  if(MAP) {
    output_dir <- check_dir(c("output", "model_fits", output_dir, "MAP"))
  } else {
    output_dir <- check_dir(c("output", "model_fits", output_dir, "full_posterior"))
  }
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

#' Define bandwidth of squared exponential kernel
#'
#' @param min_correlation minimum correlation to assume between (within-host)
#' samples
#' @param days_to_min_autocorrelation days at which squared exponential kernel decays to
#' baseline correlation of ~0.1
#' @return bandwidth parameter for squared exponential kernel
#' @export
calc_se_decay <- function(min_correlation = 0.1, days_to_min_autocorrelation = 90) {
  # Back-calculate the squared exponential bandwidth parameter by finding a
  # bandwidth that gives
  # A desired minimum correlation at the number of days specified by
  # days_to_min_autocorrelation
  rho <- sqrt(-days_to_min_autocorrelation^2/(2*log(min_correlation))) # back calculate the
                                                            # decay
  return(rho)
}

#' Define a kernel (function) over samples
#'
#' @param kernel_scale total variance for the composite kernel
#' @param diet_weight relative contribution of first three diet PCs to
#' covariance across samples
#' @param rho bandwidth for SE kernel
#' @details Composite kernel is built from (1) squared exponential kernel (base
#' autocorrelation component) and (2) seasonal kernel (periodic)
#' @return list containing kernel function and bandwidth parameter
#' @import fido
#' @export
get_Gamma <- function(kernel_scale, diet_weight, min_correlation = 0.1,
                      days_to_min_autocorrelation = 90) {
  rho <- calc_se_decay(min_correlation = min_correlation,
                       days_to_min_autocorrelation = days_to_min_autocorrelation)
  # Back-calculate the squared exponential bandwidth parameter by finding a
  # bandwidth that gives
  # A desired minimum correlation at the number of days specified by
  # SE_days_to_min_autocorrelation
  Gamma <- function(X) {
    jitter <- 1e-08
    # These are the relative variances explained by diet PCs 1, 2, 3
    # We'll scale their contribution to the sample-sample covariance
    #   proportional to these.
    var_prop <- c(0.676, 0.117, 0.061)
    var_prop <- var_prop / sum(var_prop)
    AC_component <- kernel_scale * (1 - diet_weight)
    diet_component1 <- kernel_scale * diet_weight * var_prop[1]
    diet_component2 <- kernel_scale * diet_weight * var_prop[2]
    diet_component3 <- kernel_scale * diet_weight * var_prop[3]
    K0 <- SE(X[1,,drop=F], sigma = sqrt(AC_component), rho = rho, jitter = jitter)
    K1 <- SE(X[2,,drop=F], sigma = sqrt(diet_component1), rho = 1, jitter = jitter)
    K2 <- SE(X[3,,drop=F], sigma = sqrt(diet_component2), rho = 1, jitter = jitter)
    K3 <- SE(X[4,,drop=F], sigma = sqrt(diet_component3), rho = 1, jitter = jitter)
    return(K0 + K1 + K2 + K3)
  }
  return(list(rho = rho, Gamma = Gamma))
}
