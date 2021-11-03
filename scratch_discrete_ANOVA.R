source("path_fix.R")

library(rulesoflife)
library(fido)
library(driver)
library(shapes)
library(dplyr)

output_dir <- "asv_days90_diet25_scale1"
testing <- FALSE # subset hosts and posterior samples

heatmap_cov <- function(K, label) {
  K <- cbind(1:nrow(K), K)
  colnames(K) <- c("sample1", 1:nrow(K))
  K <- pivot_longer(as.data.frame(K), !sample1, names_to = "sample2", values_to = "covariance")
  K$sample2 <- as.numeric(K$sample2)
  p <- ggplot(K, aes(x = sample1, y = sample2)) +
    geom_raster(aes(fill = covariance)) +
    scale_fill_gradient2(low = "navy", mid = "white", high = "red",
                         midpoint = 0) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank()) +
    labs(x = "",
         y = "",
         fill = "Covariance")
  show(p)
  output_dir <- check_dir(c("output", "figures"))
  ggsave(file.path(output_dir, paste0(label, ".png")),
         plot = p,
         units = "in",
         dpi = 100,
         height = 3,
         width = 4.25)
}

# ------------------------------------------------------------------------------
#   Load dynamics estimates and Frechet mean
# ------------------------------------------------------------------------------

# Load a few posterior samples for a few hosts
data <- load_data(tax_level = "ASV")
md <- data$metadata
hosts <- sort(unique(md$sname))
if(testing) {
  hosts <- sample(hosts, size = 10)
}

D <- NULL
Sigmas <- list()
for(host in hosts) {
  cat(paste0("Parsing host ", host, "...\n"))
  fit <- readRDS(file.path("output",
                           "model_fits",
                           output_dir,
                           "MAP",
                           paste0(host, ".rds")))
  # CLR
  fit <- to_clr(fit)
  if(is.null(D)) {
    D <- dim(fit$Sigma)[1]
  }
  Sigmas[[host]] <- fit$Sigma[,,1] + diag(D)*1e-06
}

# ------------------------------------------------------------------------------
#   Pull host social group / sex metadata
# ------------------------------------------------------------------------------

host_groups <- get_host_social_groups(hosts)
groups <- unique(host_groups$grp)

metadata <- load_data()$metadata
host_sex <- metadata %>%
  select(sname, sex) %>%
  distinct() %>%
  filter(sname %in% hosts) %>%
  arrange(sname)
sexes <- unique(host_sex$sex)

Frechet_mean <- readRDS(file.path("output", "Frechet_corr0_clr.rds"))

# Load group means
group_means <- list()
for(g in 1:length(groups)) {
  group_means[[g]] <- readRDS(file.path("output",
                                        paste0("Frechet_corr0_clr_group-", groups[g], ".rds")))
}

# Load group means
sex_means <- list()
for(s in 1:length(sexes)) {
  sex_means[[s]] <- readRDS(file.path("output",
                                        paste0("Frechet_corr0_clr_sex-", sexes[s], ".rds")))
}

# ------------------------------------------------------------------------------
#   ANOVA on social group
# ------------------------------------------------------------------------------

K <- length(groups)
N <- nrow(host_groups)
between_sum <- 0
for(i in 1:K) {
  grp <- groups[i]
  n_i <- sum(host_groups == grp)
  d <- distcov(group_means[[i]]$mean, Frechet_mean$mean, method = "Riemannian")^2
  between_sum <- between_sum + (n_i * d)
}
numerator <- between_sum / (K-1)

within_sum <- 0
for(i in 1:K) {
  grp <- groups[i]
  idx_i <- which(host_groups$grp == grp)
  for(j in idx_i) {
    d <- distcov(Sigmas[[host_groups$sname[j]]], group_means[[i]]$mean, method = "Riemannian")^2
    within_sum <- within_sum + d
  }
}
denominator <- within_sum / (N-K)

ratio <- numerator / denominator
ratio

# plot(density(rf(1000, K-1, N-K)))
1 - pf(ratio, K-1, N-K)

# Visualize the differences in scale!
# heatmap_cov(Frechet_mean$mean, "Frechet_mean")
# heatmap_cov(group_means[[1]]$mean, "Frechet_mean_grp-1.1")
# heatmap_cov(Sigmas[[which(host_groups$sname == "COO")]], "Frechet_mean_grp-1.1-COO")

# ------------------------------------------------------------------------------
#   ANOVA on sex
# ------------------------------------------------------------------------------

K <- length(sexes)
N <- nrow(host_sex)
between_sum <- 0
for(i in 1:K) {
  sex <- sexes[i]
  n_i <- sum(host_sex == sex)
  d <- distcov(sex_means[[i]]$mean, Frechet_mean$mean, method = "Riemannian")^2
  between_sum <- between_sum + (n_i * d)
}
numerator <- between_sum / (K-1)

within_sum <- 0
for(i in 1:K) {
  sex <- sexes[i]
  idx_i <- which(host_sex$sex == sex)
  for(j in idx_i) {
    d <- distcov(Sigmas[[host_sex$sname[j]]], sex_means[[i]]$mean, method = "Riemannian")^2
    within_sum <- within_sum + d
  }
}
denominator <- within_sum / (N-K)

ratio <- numerator / denominator
ratio

# plot(density(rf(1000, K-1, N-K)))
1 - pf(ratio, K-1, N-K)
