source("path_fix.R")

library(rulesoflife)
library(tidyverse)
library(coxme) # for lmekin
library(MASS)
library(matrixsampling)
library(driver)

source("ggplot_fix.R")

# ------------------------------------------------------------------------------
#   Simulate data and fit
# ------------------------------------------------------------------------------

if(FALSE) {
  # Simulate data to test the usage/performance of lmekin
  n_pairs <- 30
  n_hosts <- 10

  rho <- 0.5 # within-group correlation / between group anticorrelation
  K <- matrix(-rho, n_hosts, n_hosts)
  K[1:(n_hosts/2),1:(n_hosts/2)] <- rho
  K[(n_hosts/2+1):n_hosts,(n_hosts/2+1):n_hosts] <- rho
  diag(K) <- 1
  K_scale <- 1
  K <- K_scale * K
  # plot_kernel_or_cov_matrix(K)

  mu_mat <- matrix(rnorm(n_pairs), n_hosts, n_pairs, byrow = TRUE)
  y <- rmatrixnormal(1, mu_mat, K, diag(n_pairs))[,,1]
  # y <- rmatrixnormal(1, mu_mat, diag(n_hosts), diag(n_pairs))[,,1]

  data <- data.frame(y = c(y),
                     host = rep(1:n_hosts, n_pairs),
                     pair = rep(1:n_pairs, each = n_hosts))
  data <- data %>%
    mutate(host_pair = paste0(host, "_", pair))
  data$host <- factor(data$host)
  data$pair <- factor(data$pair)
  data$host_pair <- factor(data$host_pair, levels = data$host_pair)

  # Linearize all...
  K_large <- kronecker(diag(n_pairs), K)
  fit <- lmekin(y ~ pair + (1 | host_pair), data = data, varlist = list(K_large))

  # Get fixed and random effect estimates
  mu_pred <- matrix(fit$coefficients$fixed, n_hosts, n_pairs, byrow = TRUE)
  re_pred <- matrix(NA, n_hosts, n_pairs)
  for(re_idx in 1:length(fit$coefficients$random$host_pair)) {
    re <- fit$coefficients$random$host_pair[re_idx]
    pieces <- as.numeric(strsplit(names(re), "_")[[1]])
    re_pred[pieces[1], pieces[2]] <- unname(re)
  }

  # Visualize true (L), fixed effect prediction (C), fixed + RE prediction (R)
  orderings <- plot_rug(y, save_name = "LMM_simulated")
  plot_rug(mu_pred,
           canonical_col_order = orderings$col_order,
           canonical_row_order = orderings$row_order,
           save_name = "LMM_fe_estimate")
  plot_rug(mu_pred + re_pred,
           canonical_col_order = orderings$col_order,
           canonical_row_order = orderings$row_order,
           save_name = "LMM_me_estimate")

  # Total variation in the data
  var_y <- var(c(y))
  # Residual variance (fixed + RE)
  var_model_fe <- var(c(y - mu_pred))
  var_model_re <- var(c(y - (mu_pred + re_pred)))

  (var_y - var_model_fe)/var_y
  (var_y - var_model_re)/var_y
}

# ------------------------------------------------------------------------------
#   Real data
# ------------------------------------------------------------------------------

subset_pairs <- TRUE

cat("Loading 'rug' data...\n")
rug_obj <- summarize_Sigmas(output_dir = "asv_days90_diet25_scale1")

# Code pulled from R/utility.R
# PEDIGREE
pedigree <- readRDS(file.path("input", "pedigree_56hosts.RDS"))
host_order_rug <- data.frame(host = rug_obj$hosts,
                             index_rug = 1:length(rug_obj$hosts))
host_order_ped <- data.frame(host = rownames(pedigree),
                             index_ped = 1:nrow(pedigree))
host_reordering <- left_join(host_order_rug,
                             host_order_ped,
                             by = "host")$index_ped
pedigree2 <- pedigree[host_reordering,host_reordering]
K <- pedigree2
K <- cov2cor(K)

# CLR BASELINE
# data <- load_data(tax_level = "ASV")
# clr.counts <- clr_array(data$counts + 0.5, parts = 1)
# # Get each host's "baseline" composition (avg. CLR composition)
# hosts <- rug_obj$hosts
# host_baselines <- matrix(NA, length(hosts), nrow(clr.counts))
# for(i in 1:length(hosts)) {
#   host <- hosts[i]
#   host_data <- clr.counts[, data$metadata$sname == host]
#   host_baselines[i,] <- rowMeans(host_data)
# }
# K <- as.matrix(dist(host_baselines))
# K <- max(K) - K
# K <- cov2cor(K)

# Subset for testing
# host_idx <- 31:51
# n_hosts <- length(host_idx)
# n_pairs <- 20
# K <- K[host_idx,host_idx]

if(subset_pairs) {
  n_pairs <- 1000
  y <- rug_obj$rug[,sample(1:ncol(rug_obj$rug), size = n_pairs)]
  n_hosts <- length(rug_obj$hosts)
  save_filename <- file.path(check_dir(c("output")),
                             "pedigree_LMM_predicted_subset.rds")
} else {
  n_pairs <- ncol(rug_obj$rug)
  y <- rug_obj$rug
  n_hosts <- length(rug_obj$hosts)
  save_filename <- file.path(check_dir(c("output")),
                             "pedigree_LMM_predicted.rds")
}

data <- data.frame(y = c(y),
                   host = rep(1:n_hosts, n_pairs),
                   pair = rep(1:n_pairs, each = n_hosts))
data <- data %>%
  mutate(host_pair = paste0(host, "_", pair))
data$host <- factor(data$host)
data$pair <- factor(data$pair)
data$host_pair <- factor(data$host_pair, levels = data$host_pair)

# Linearize all...
K_large <- kronecker(diag(n_pairs), K)
start <- Sys.time()
cat("Starting model fit...\n")
fit <- lmekin(y ~ pair + (1 | host_pair), data = data, varlist = list(K_large))
cat("Fit time:", Sys.time() - start,"\n")

# Get fixed and random effect estimates
mu_pred <- matrix(fit$coefficients$fixed, n_hosts, n_pairs, byrow = TRUE)
re_pred <- matrix(NA, n_hosts, n_pairs)
for(re_idx in 1:length(fit$coefficients$random$host_pair)) {
  re <- fit$coefficients$random$host_pair[re_idx]
  pieces <- as.numeric(strsplit(names(re), "_")[[1]])
  re_pred[pieces[1], pieces[2]] <- unname(re)
}

saveRDS(list(fit = fit,
             y = y,
             fe_pred = mu_pred,
             me_pred = mu_pred + re_pred),
        file = save_filename)

# Total variation in the data
var_y <- var(c(y))
# Residual variance (fixed + RE)
var_model_fe <- var(c(y - mu_pred))
var_model_re <- var(c(y - (mu_pred + re_pred)))

var_expl_baseline <- (var_y - var_model_fe)/var_y
var_expl_full <- (var_y - var_model_re)/var_y

cat("Additional percent variance explained:", round((var_expl_full - var_expl_baseline)*100, 3), "\n")
