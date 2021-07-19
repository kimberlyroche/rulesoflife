source("path_fix.R")

library(rulesoflife)
library(fido)
library(shapes)
library(matrixsampling)
library(tidyverse)
library(lme4)

testing <- FALSE
use_corr <- FALSE

# ------------------------------------------------------------------------------
#   Get mean covariance for all hosts
# ------------------------------------------------------------------------------

output_dir <- "asv_days90_diet25_scale1"
output_dir_full <- check_dir(c("output", "model_fits", output_dir, "full_posterior"))
file_list <- list.files(path = output_dir_full, pattern = "*.rds")

# Get taxa number and posterior sample number
fit <- readRDS(file.path(output_dir_full, file_list[1]))
D <- fit$D
iter <- fit$iter
limit <- length(file_list)

mean_matrices <- list()
if(use_corr) {
  save_path <- file.path("output", "mean_matrices_corr.rds")
} else {
  save_path <- file.path("output", "mean_matrices.rds")
}
if(!file.exists(save_path)) {
  # Full posteriors
  for(f in 1:limit) {
    file <- file_list[f]
    fit <- readRDS(file.path(output_dir_full, file))
    # Convert to CLR
    if(fit$coord_system != "clr") {
      fit <- to_clr(fit)
    }
    cat(paste0("Saving covariance for ", fit$sname, "...\n"))
    mean_matrices[[fit$sname]] <- apply(fit$Sigma, c(1,2), mean)
  }
  saveRDS(mean_matrices, file = save_path)
} else {
  mean_matrices <- readRDS(save_path)
}

use_indices <- 1:D

# Add jitter, otherwise we'll get errors about tiny negative eigenvalues
mm2 <- lapply(mean_matrices, function(x) {
  if(use_corr) {
    x <- cov2cor(x)
  }
  x <- x + diag(D)*1e-06
  if(testing) {
    x[use_indices,use_indices]
  } else {
    x
  }
})

# Convert to array
arr_ext <- D
if(testing) {
  arr_ext <- nrow(mm2[[1]])
}
S <- array(as.numeric(unlist(mm2)), dim = c(arr_ext, arr_ext, limit))

if(testing) {
  # Just sub element-wise mean
  Frechet_mean <- list(mean = apply(S, c(1,2), mean))
} else {
  if(use_corr) {
    save_path <- file.path("output", "Frechet_mean_corr.rds")
  } else {
    save_path <- file.path("output", "Frechet_mean.rds")
  }
  if(!file.exists(save_path)) {
    cat("Estimating Frechet mean...\n")
    start <- Sys.time()
    Frechet_mean <- estcov(S, method = "Riemannian", weights = 1)
    diff <- Sys.time() - start
    cat(paste0("Estimate took ",
               round(as.numeric(diff), 3),
               " ",
               attr(diff, "units"),
               "\n"))
    cat("Saving mean estimate...\n")
    saveRDS(Frechet_mean, file = save_path)
  } else {
    Frechet_mean <- readRDS(save_path)
  }
}

plot_kernel_or_cov_matrix(Frechet_mean$mean)
ggsave(file.path("output", "figures",
                 paste0("Frechet_", ifelse(use_corr, "cor", "cov"),
                        ".png")),
       dpi = 100,
       units = "in",
       height = 4,
       width = 5)

plot_kernel_or_cov_matrix(S[,,sample(1:limit, size = 1)])
ggsave(file.path("output", "figures",
                 paste0("Frechet_", ifelse(use_corr, "cor", "cov"),
                        "_sample.png")),
       dpi = 100,
       units = "in",
       height = 4,
       width = 5)

# ------------------------------------------------------------------------------
#   Calculate distribution of distances to mean
# ------------------------------------------------------------------------------

mean_dist <- c()
for(i in 1:limit) {
  mean_dist <- c(mean_dist, distcov(Frechet_mean$mean, S[,,i], method = "Riemannian")^2)
}

plot_df <- data.frame(x = mean_dist)
ggplot(plot_df, aes(x = x)) +
  geom_histogram(color = "#ffffff") +
  labs(x = paste0("Riemannian distance to Frechet mean (",
                  ifelse(use_corr, "correlation", "covariance"),
                  ")"))
ggsave(file.path("output", "figures", paste0("distances_", ifelse(use_corr, "cor", "cov"), ".png")),
       dpi = 100,
       units = "in",
       height = 4,
       width = 5)

# ------------------------------------------------------------------------------
#   Distances across individuals + random effects model
# ------------------------------------------------------------------------------

pedigree <- readRDS(file.path("input", "pedigree_56hosts.RDS"))

# Make into data.frame
pedigree <- cbind(host1 = rownames(pedigree), as.data.frame(pedigree))
pedigree <- pedigree %>%
  pivot_longer(!host1, names_to = "host2", values_to = "relatedness") %>%
  distinct() %>%
  filter(host1 != host2)

# Apply distcov to all combos to generate new column
# This takes ~1-2 min.
host_list <- names(mean_matrices)
pedigree$dist <- sapply(1:nrow(pedigree), function(i) {
  distcov(S[,,which(host_list == pedigree$host1[i])],
          S[,,which(host_list == pedigree$host2[i])],
          method = "Riemannian")
})

# Plot Riemannian distance vs. pedigree (R^2)
ggplot(pedigree, aes(x = relatedness, y = dist)) +
  geom_point(size = 2, shape = 21, fill = "#bbbbbb") +
  labs(x = "host-pair relatedness",
       y = "distance between host-pair covariance matrices")

cat(paste0("R^2: ", round(cor(pedigree$relatedness, pedigree$dist), 3), "\n"))

# Apply lmer
res <- lm(dist ~ 1 + relatedness, data = pedigree)
pred <- predict(res, newdata = pedigree)

var(pedigree$dist)
var(res$residuals)

summary(res)

# ------------------------------------------------------------------------------
#   Get mean covariance for each social group
# ------------------------------------------------------------------------------

# ID social groups
social_groups <- get_host_social_groups(names(mean_matrices))
unique_social_groups <- as.character(unique(social_groups$grp))

group_Frechet_means <- list()
for(grp in unique_social_groups) {
  if(testing) {
    # Just sub element-wise mean
    Frechet_mean <- list(mean = apply(S[,,social_groups$grp == grp], c(1,2), mean))
    group_Frechet_means[[grp]] <- Frechet_mean
  } else {
    save_path <- file.path("output", paste0("Frechet_mean_", grp, ".rds"))
    if(!file.exists(save_path)) {
      cat(paste0("Estimating Frechet mean of group ", grp, "...\n"))
      start <- Sys.time()
      Frechet_mean <- estcov(S[,,social_groups$grp == grp], method = "Riemannian", weights = 1)
      diff <- Sys.time() - start
      cat(paste0("Estimate took ",
                 round(as.numeric(diff), 3),
                 " ",
                 attr(diff, "units"),
                 "\n"))
      cat("Saving mean estimate...\n")
      saveRDS(Frechet_mean, file = save_path)
      group_Frechet_means[[grp]] <- Frechet_mean
    } else {
      Frechet_mean <- readRDS(save_path)
      group_Frechet_means[[grp]] <- Frechet_mean
    }
  }
}

# ------------------------------------------------------------------------------
#   Calculate F-statistic
# ------------------------------------------------------------------------------

K <- length(unique_social_groups)
N <- length(social_groups$grp)
between_sum <- 0
for(i in 1:K) {
  grp <- unique_social_groups[i]
  n_i <- sum(social_groups == grp)
  d <- distcov(group_Frechet_means[[grp]]$mean, Frechet_mean$mean, method = "Riemannian")^2
  between_sum <- between_sum + (n_i * d)
}
numerator <- between_sum / (K-1)

within_sum <- 0
for(i in 1:K) {
  grp <- unique_social_groups[i]
  idx_i <- which(social_groups$grp == grp)
  for(j in idx_i) {
    d <- distcov(S[,,j], group_Frechet_means[[grp]]$mean, method = "Riemannian")^2
    within_sum <- within_sum + d
  }
}
denominator <- within_sum / (N-K)

ratio <- numerator / denominator
ratio

# plot(density(rf(1000, K-1, N-K)))
1 - pf(ratio, K-1, N-K)

length(group_Frechet_means)

plot_kernel_or_cov_matrix(Frechet_mean$mean)

plot_kernel_or_cov_matrix(group_Frechet_means[[unique_social_groups[1]]]$mean)
plot_kernel_or_cov_matrix(group_Frechet_means[[unique_social_groups[2]]]$mean)
plot_kernel_or_cov_matrix(group_Frechet_means[[unique_social_groups[3]]]$mean)
plot_kernel_or_cov_matrix(group_Frechet_means[[unique_social_groups[4]]]$mean)
plot_kernel_or_cov_matrix(group_Frechet_means[[unique_social_groups[5]]]$mean)

# ------------------------------------------------------------------------------
#   Repeat for sex
# ------------------------------------------------------------------------------

sexes <- as.data.frame(load_data()$metadata %>%
                         select(sname, sex) %>%
                         distinct())$sex
unique_sexes <- unique(sexes)

group_Frechet_means <- list()
for(grp in unique_sexes) {
  if(testing) {
    # Just sub element-wise mean
    Frechet_mean <- list(mean = apply(S[,,sexes == grp], c(1,2), mean))
    group_Frechet_means[[grp]] <- Frechet_mean
  } else {
    stop("Only implemented for approx. case!")
  }
}

K <- 2
N <- length(sexes)
between_sum <- 0
for(i in 1:K) {
  grp <- unique_sexes[i]
  n_i <- sum(sexes == grp)
  d <- distcov(group_Frechet_means[[grp]]$mean, Frechet_mean$mean, method = "Riemannian")^2
  between_sum <- between_sum + (n_i * d)
}
numerator <- between_sum / (K-1)

within_sum <- 0
for(i in 1:K) {
  grp <- unique_sexes[i]
  idx_i <- which(sexes == grp)
  for(j in idx_i) {
    d <- distcov(S[,,j], group_Frechet_means[[grp]]$mean, method = "Riemannian")^2
    within_sum <- within_sum + d
  }
}
denominator <- within_sum / (N-K)

ratio <- numerator / denominator
ratio

# plot(density(rf(1000, K-1, N-K)))
1 - pf(ratio, K-1, N-K)

# ------------------------------------------------------------------------------
#   Permute each individual's covariance matrix and compute distance to the
#   mean (here I'm using element-wise mean for time but it should be close)
# ------------------------------------------------------------------------------

perm_S <- S
for(i in 1:dim(perm_S)[3]) {
  scramble_idx <- sample(1:dim(perm_S)[1])
  perm_S[,,i] <- perm_S[scramble_idx,scramble_idx,i]
}

# Element-wise mean
mean_perm_S <- apply(perm_S, c(1,2), mean)

mean_perm_dist <- c()
for(i in 1:limit) {
  mean_perm_dist <- c(mean_perm_dist, distcov(mean_perm_S, perm_S[,,i], method = "Riemannian")^2)
}

plot_df <- data.frame(x = mean_dist, type = "true")
plot_df <- rbind(plot_df,
                 data.frame(x = mean_perm_dist, type = "permuted"))
plot_df$type <- factor(plot_df$type, levels = c("true", "permuted"))
levels(plot_df$type) <- c("orig. estimates", "permuted")
ggplot(plot_df, aes(x = x, fill = type)) +
  geom_histogram(color = "#ffffff") +
  # scale_fill_manual(values = c("#865EC4", "#7BBF37")) +
  scale_fill_manual(values = c("#BC8999", "#9BB49D")) +
  xlim(c(0, max(plot_df$x))) +
  labs(fill = "Data type",
       x = "distance to mean")

1 - mean(mean_dist)/mean(mean_perm_dist)



# pair_dist <- c()
# pairs <- combn(1:limit, 2)
# for(i in 1:ncol(pairs)) {
#   pair_dist <- c(pair_dist, distcov(S[,,pairs[1,i]], S[,,pairs[2,i]], method = "Riemannian")^2)
# }
#
# plot_df <- data.frame(dist = mean_dist, type = "mean")
# plot_df <- rbind(plot_df, data.frame(dist = pair_dist, type = "pairwise"))
# plot_df$type <- factor(plot_df$type, levels = c("mean", "pairwise"))
# levels(plot_df$type) <- c("distances to mean", "pairwise distances")
# # ggplot(plot_df, aes(x = dist, fill = type)) +
# #   geom_histogram(color = "#ffffff")
# ggplot(plot_df, aes(x = type, y = dist)) +
#   geom_boxplot() +
#   ylim(c(0, max(plot_df$dist))) +
#   labs(x = "",
#        y = "Riemannian distance")
#
# # Calculate and F-like statistic
# mean(mean_dist) / mean(pair_dist)

# # PCA
# plot_df <- data.frame(PC1 = Frechet_mean$pco[,1], PC2 = Frechet_mean$pco[,2])
# ggplot(plot_df, aes(x = PC1, y = PC2)) +
#   geom_point(size = 3, shape = 21, fill = "red")
