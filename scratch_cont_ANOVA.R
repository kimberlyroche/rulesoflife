source("path_fix.R")

library(rulesoflife)
library(fido)
library(driver)
library(shapes)
library(tidyverse)
library(vegan)

output_dir <- "asv_days90_diet25_scale1"
testing <- TRUE # subset hosts and posterior samples

# ------------------------------------------------------------------------------
#   Load baseline and dynamics estimates and calculate sets of distances between
#   them
#
#   I'm using the ALR for now but we probably want to consider which data
#   representation is best.
# ------------------------------------------------------------------------------

# Load a few posterior samples for a few hosts
data <- load_data(tax_level = "ASV")
md <- data$metadata
hosts <- sort(unique(md$sname))
if(testing) {
  hosts <- sample(hosts, size = 20)
}

Sigmas <- list()
n_iter <- NULL
for(host in hosts) {
  cat(paste0("Parsing host ", host, "...\n"))
  fit <- readRDS(file.path("output",
                           "model_fits",
                           output_dir,
                           "full_posterior",
                           paste0(host, ".rds")))
  # CLR
  fit <- to_clr(fit)
  if(is.null(n_iter)) {
    if(testing) {
      n_iter <- 20
    } else {
      n_iter <- dim(fit$Sigma)[3]
    }
  }
  Sigmas[[host]] <- fit$Sigma[,,1:n_iter]
}
D <- dim(Sigmas[[1]])[1]

# Calculate samples of the dynamics distance
host_combos <- combn(length(hosts), m = 2)
dynamics_distances <- matrix(NA, ncol(host_combos), n_iter)
for(i in 1:n_iter) {
  cat(paste0("Calculating distances for sample ", i, "\n"))
  dynamics_distances[,i] <- apply(host_combos, 2, function(pair) {
    distcov(Sigmas[[pair[1]]][,,i] + diag(D)*1e-06,
            Sigmas[[pair[2]]][,,i] + diag(D)*1e-06,
            method = "Riemannian")
  })
}

# ------------------------------------------------------------------------------
#   Load and test PEDIGREE
# ------------------------------------------------------------------------------

# Load pedigree
pedigree <- readRDS(file.path("input", "pedigree_56hosts.RDS"))
mapping <- data.frame(sname = names(Sigmas)) %>%
  left_join(data.frame(sname = rownames(pedigree), order = 1:nrow(pedigree)), by = "sname")
pedigree <- pedigree[mapping$order,mapping$order] # subset and order as Sigmas
ped_dist <- 1 - pedigree

x <- ped_dist[lower.tri(ped_dist, diag = FALSE)]
x <- scale(x)

sampled_corr_ped <- c()
for(i in 1:n_iter) {
  y <- dynamics_distances[,i]
  y <- scale(y)
  sampled_corr_ped <- c(sampled_corr_ped, cor(x, y)**2)
}

# Plot a single posterior sample
pl <- ggplot(data.frame(x = x, y = y),
       aes(x = x, y = y)) +
  geom_point(size = 3, shape = 21, fill = "#dddddd") +
  # geom_smooth(method = "lm", color = "black") +
  theme_bw() +
  labs(x = "pedigree dissimilarity (scaled)",
       y = "dynamics distance (scaled)",
       title = paste0("Dynamics x pedigree: R^2 = ", round(mean(sampled_corr_ped), 3), " +/- ", round(sd(sampled_corr_ped), 3)))
show(pl)
ggsave(file.path("output", "figures", "pedigree_R2.png"),
       plot = pl,
       dpi = 100,
       units = "in",
       height = 4,
       width = 5)

# ------------------------------------------------------------------------------
#   Mantel results from the same
# ------------------------------------------------------------------------------

# Tedious and bad code but... render (by re-calculating) these distances as a
# matrix
dynamics_distances <- matrix(0, length(hosts), length(hosts))
for(i in 1:ncol(host_combos)) {
  cat(paste0("Host combo ", i, " / ", ncol(host_combos), "\n"))
  h1 <- host_combos[1,i]
  h2 <- host_combos[2,i]
  dynamics_distances[h1, h2] <- distcov(Sigmas[[h1]][,,1] + diag(D)*1e-06,
                                        Sigmas[[h2]][,,1] + diag(D)*1e-06,
                                        method = "Riemannian")
  dynamics_distances[h2, h1] <- dynamics_distances[h1, h2]
}

fit <- mantel(ped_dist, dynamics_distances)
cat(paste0("Mantel p-value (dynamics x composition): ", round(fit$signif, 5), "\n"))

# ------------------------------------------------------------------------------
#   Load and test BASELINE COMPOSITION
# ------------------------------------------------------------------------------

# Calculate baseline composition -- mean Eta for each host
baselines <- matrix(NA, length(hosts), D)
for(i in 1:length(hosts)) {
  # Average of ALR samples for this host is their baseline
  host_samples <- clr_array(data$counts[,which(md$sname == hosts[i])] + 0.5,
                            parts = 1)
  baselines[i,] <- apply(host_samples, 1, mean)
}

x <- c(dist(baselines))
x <- scale(x)
sampled_corr_base <- c()
for(i in 1:n_iter) {
  y <- dynamics_distances[,i]
  y <- scale(y)
  sampled_corr_base <- c(sampled_corr_base, cor(x, y)**2)
}

# Plot a single posterior sample
pl <- ggplot(data.frame(x = x, y = y),
             aes(x = x, y = y)) +
  geom_point(size = 3, shape = 21, fill = "#dddddd") +
  # geom_smooth(method = "lm", color = "black") +
  theme_bw() +
  labs(x = "baseline compositional distance (scaled)",
       y = "dynamics distance (scaled)",
       title = paste0("Dynamics x composition: R^2 = ",
                      round(mean(sampled_corr_base), 3),
                      " +/- ",
                      round(sd(sampled_corr_base), 3)))
show(pl)
ggsave(file.path("output", "figures", "composition_R2.png"),
       plot = pl,
       dpi = 100,
       units = "in",
       height = 4,
       width = 5)

# ------------------------------------------------------------------------------
#   Mantel results from the same
# ------------------------------------------------------------------------------

fit <- mantel(dist(baselines), dynamics_distances)
cat(paste0("Mantel p-value (dynamics x composition): ", round(fit$signif, 5), "\n"))

