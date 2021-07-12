library(rulesoflife)
library(fido)
library(shapes)
library(tidyverse)

output_dir <- "asv_days90_diet25_scale1"
n_iter <- 10

# Load a few posterior samples for a few hosts
data <- load_data(tax_level = "ASV")
md <- data$metadata
hosts <- sort(unique(md$sname))

Sigmas <- NULL
for(host in hosts) {
  cat(paste0("Parsing host ", host, "...\n"))
  fit <- readRDS(file.path("output",
                           "model_fits",
                           output_dir,
                           "full_posterior",
                           paste0(host, ".rds")))
  fit.clr <- to_clr(fit)
  Sigmas[[host]] <- fit.clr$Sigma[,,sample(1:fit$iter, size = n_iter)]
}
D <- dim(Sigmas[[1]])[1]

# Load Frechet mean (need to re-run this on proper MAP estimates)
# Nvm; no need here
# Frechet_mean <- readRDS(file.path("output", "Frechet_mean.rds"))

# Load pedigree
pedigree <- readRDS(file.path("input", "pedigree_56hosts.RDS"))
# Reorder
pedigree <- pedigree[order(rownames(pedigree)),order(colnames(pedigree))]

# Calculate dynamics distance
s <- 1 # sample i
limit <- 1000
# limit <- length(host_combos)
host_combos <- combn(length(hosts), m = 2)
dynamics_distances <- numeric(limit)
pedigree_vec <- numeric(limit)
for(i in 1:limit) {
  if(i %% 100 == 0) {
    cat(paste0("Iteration ", i, "\n"))
  }
  host1 <- hosts[host_combos[1,i]]
  host2 <- hosts[host_combos[2,i]]
  dynamics_distances[i] <- distcov(Sigmas[[host1]][,,s] + diag(D)*1e-06,
                                   Sigmas[[host2]][,,s] + diag(D)*1e-06,
                                   method = "Riemannian")
  pedigree_vec[i] <- pedigree[host_combos[1,i],host_combos[2,i]]
}
dynamics_distances <- scale(dynamics_distances)
pedigree_vec <- scale(pedigree_vec)

cor(dynamics_distances, pedigree_vec)**2

ggplot(data.frame(x = dynamics_distances, y = pedigree_vec),
       aes(x = x, y = y)) +
  geom_point(size = 3, shape = 21, fill = "#aaaaaa") +
  theme_bw() +
  labs(x = "dynamics distance",
       y = "pedigree similarity")

# ------------------------------------------------------------------------------
#   Is a Mantel test even possible here?
# ------------------------------------------------------------------------------

library(vegan)

# < 1 min. runtime
dynamics_distances <- matrix(0, length(hosts), length(hosts))
for(i in 1:ncol(host_combos)) {
  h1 <- host_combos[1,i]
  h2 <- host_combos[2,i]
  dynamics_distances[h1, h2] <- distcov(Sigmas[[h1]][,,1] + diag(D)*1e-06,
          Sigmas[[h2]][,,1] + diag(D)*1e-06,
          method = "Riemannian")
  dynamics_distances[h2, h1] <- dynamics_distances[h1, h2]
}

pedigree_dissimilarities <- 1 - pedigree

mantel(pedigree_dissimilarities, dynamics_distances)

# Next: Is using R^2 as a plug-in for variance explained reasonable?
