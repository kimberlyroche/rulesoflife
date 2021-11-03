library(rulesoflife)
library(fido)
library(driver)
library(shapes)
library(tidyverse)
library(vegan) # for Mantel test
library(pedtools)
library(kinship2)
library(cowplot)

output_dir <- "asv_days90_diet25_scale1"
p <- 999 # number of permutations

# ------------------------------------------------------------------------------
#   Load data, compute dynamics distances
#
#   Note: I'm using MAP estimates here
# ------------------------------------------------------------------------------

# Load a few posterior samples for a few hosts
data <- load_data(tax_level = "ASV")
md <- data$metadata
hosts <- sort(unique(md$sname))
N <- length(hosts)

Sigmas <- list()
n_iter <- NULL
for(host in hosts) {
  cat(paste0("Parsing host ", host, "...\n"))
  fit <- readRDS(file.path("output",
                           "model_fits",
                           output_dir,
                           "MAP",
                           paste0(host, ".rds")))
  # CLR
  fit <- to_clr(fit)
  Sigmas[[host]] <- cov2cor(fit$Sigma[,,1])
}
D <- dim(Sigmas[[1]])[1]

# Calculate samples of the dynamics distance
host_combos <- combn(length(hosts), m = 2)
dynamics_distances <- matrix(NA, N, N)
for(i in 1:ncol(host_combos)) {
  if(i %% 100 == 0) {
    cat(paste0("Host combo ", i, " / ", ncol(host_combos), "\n"))
  }
  a <- host_combos[1,i]
  b <- host_combos[2,i]
  d <- distcov(Sigmas[[a]] + diag(D)*1e-06,
               Sigmas[[b]] + diag(D)*1e-06,
               method = "Riemannian")
  dynamics_distances[a,b] <- d
  dynamics_distances[b,a] <- d
}

dynamics_dist_vec <- dynamics_distances[lower.tri(dynamics_distances)]

# ------------------------------------------------------------------------------
#
#   Pedigree associations
#
# ------------------------------------------------------------------------------

# Load pedigree
pedigree <- readRDS(file.path("input", "pedigree_56hosts.RDS"))
mapping <- data.frame(sname = names(Sigmas)) %>%
  left_join(data.frame(sname = rownames(pedigree), order = 1:nrow(pedigree)), by = "sname")
pedigree <- pedigree[mapping$order,mapping$order] # subset and order as Sigmas
ped_dist <- 1 - pedigree
ped_vec <- ped_dist[lower.tri(ped_dist, diag = FALSE)]

obs_Rsq_ped <- c(cor(scale(ped_vec), scale(dynamics_dist_vec))**2)

# R^2 plot
p1_ped <- ggplot(data.frame(x = scale(ped_vec), y = scale(dynamics_dist_vec)),
                 aes(x = x, y = y)) +
  # geom_smooth(color = "gray", alpha = 0.5) +
  geom_point(size = 3, shape = 21, fill = "#999999") +
  theme_bw() +
  labs(x = "kinship dissimilarity (1 - % relatedness)",
       y = "host-host dynamics distance (Riemannian)")

# ------------------------------------------------------------------------------
#   Mantel test
#
#   This answers the question: How much variation is explained by X on distances
#   in dynamics *in a population this related*?
# ------------------------------------------------------------------------------

res1_ped <- mantel(as.dist(ped_dist), as.dist(dynamics_distances))
cat(paste0("P-value from Mantel test (kinship): ", round(res1_ped$signif, 3), "\n"))

# ------------------------------------------------------------------------------
#   Null distribution - kinship/pedigree
#
#   Generalized approximately to populations with any amount of relatedness
# ------------------------------------------------------------------------------

# Note: If the number of founders (which must apparently be << N) is not
# specified, this is drawn from a Poisson distribution. A smaller number of
# founders makes the population more related.

simulate_K <- function(g) {
  # Simulate a random pedigree
  temp <- randomPed(g = g)
  # Combine components
  if(any(str_detect(names(temp), "^_comp"))) {
    temp2 <- NULL
    for(i in 1:length(temp)) {
      temp2 <- rbind(as.data.frame(temp2),
                     as.data.frame(temp[[i]]))
    }
    temp2 <- temp2[order(as.numeric(temp2$id)),]
  } else {
    temp2 <- as.data.frame(temp)
  }
  # Convert it to a kinship matrix
  temp3 <- data.frame(id = temp2$id,
                      mom = temp2$mid,
                      dad = temp2$fid,
                      sex = temp2$sex)
  f <- sum(temp3$dad == 0 | temp3$mom == 0)
  temp3$mom[1:f] <- NA
  temp3$dad[1:f] <- NA
  K <- kinship(with(temp3, pedigree(id, dad, mom, sex)))*2
  K <- K[(f+1):nrow(K),(f+1):ncol(K)]
  K
}

draw_random_Rsq_ped <- function(g) {
  K <- simulate_K(g)
  ped_dist <- 1 - K
  ped_vec <- ped_dist[lower.tri(ped_dist, diag = FALSE)]
  c(cor(scale(ped_vec), scale(dynamics_dist_vec))**2)
}

null_Rsq_ped <- sapply(1:p, function(pp) draw_random_Rsq_ped(g))

p2_ped <- ggplot(data.frame(x = null_Rsq_ped), aes(x = x)) +
  geom_histogram(color = "white") +
  theme_bw() +
  labs(x = "R-squared (random kinship dissimilarity x dynamics distance)")

p3_ped <- plot_grid(p1_ped, p2_ped, ncol = 2,
                rel_widths = c(1,1), labels = c("a", "b"),
                scale = 0.9, label_size = 20)
show(p3_ped)

ggsave(file.path("output", "figures", "ANOVA_pedigree_dynamics.svg"),
       plot = p3,
       dpi = 100,
       units = "in",
       height = 4.5,
       width = 10)

# Calculate a pseudo pvalue from quantiles
res2 <- sum(null_Rsq_ped > obs_Rsq_ped) / p
cat(paste0("P-value from random null test: ", round(res2, 3), "\n"))

# ------------------------------------------------------------------------------
#
#   Baseline composition associations
#
# ------------------------------------------------------------------------------

# Calculate baseline composition as mean CLR-transformed composition for each host
baselines <- matrix(NA, length(hosts), D)
for(i in 1:length(hosts)) {
  # Average of ALR samples for this host is their baseline
  host_samples <- clr_array(data$counts[,which(md$sname == hosts[i])] + 0.5,
                            parts = 1)
  baselines[i,] <- apply(host_samples, 1, mean)
}

baseline_distances <- dist(baselines)
baseline_dist_vec <- c(baseline_distances)

obs_Rsq_comp <- c(cor(baseline_dist_vec, dynamics_dist_vec)**2)

# R^2 plot
p1_comp <- ggplot(data.frame(x = scale(baseline_dist_vec), y = scale(dynamics_dist_vec)),
                  aes(x = x, y = y)) +
  # geom_smooth(color = "gray", alpha = 0.5) +
  geom_point(size = 3, shape = 21, fill = "#999999") +
  theme_bw() +
  labs(x = "Aitchison distance over baselines",
       y = "host-host dynamics distance (Riemannian)")

# ------------------------------------------------------------------------------
#   Mantel test
#
#   This answers the question: How much variation is explained by X on distances
#   in dynamics *in a population this related*?
# ------------------------------------------------------------------------------

res1_comp <- mantel(baseline_distances, as.dist(dynamics_distances))
cat(paste0("P-value from Mantel test (composition): ", round(res1_comp$signif, 3), "\n"))

# ------------------------------------------------------------------------------
#   Null distribution - baseline composition
#
#   Generalized approximately to populations with any amount of relatedness
# ------------------------------------------------------------------------------

draw_random_Rsq_comp <- function() {
  # Calculate baseline composition as mean CLR-transformed composition for each host
  # Shuffle existing baselines by shuffling same taxa along all hosts
  new_baselines <- baselines
  for(i in 1:ncol(baselines)) {
    new_baselines[,i] <- sample(new_baselines[,i])
  }
  new_baseline_dist_vec <- c(dist(new_baselines))
  c(cor(new_baseline_dist_vec, dynamics_dist_vec)**2)
}

null_Rsq_comp <- numeric(p)
for(pp in 1:p) {
  null_Rsq_comp[pp] <- draw_random_Rsq_comp()
}

p2_comp <- ggplot(data.frame(x = null_Rsq_comp), aes(x = x)) +
  geom_histogram(color = "white") +
  theme_bw() +
  labs(x = "R-squared (random baseline x dynamics distance)")

p3_comp <- plot_grid(p1_comp, p2_comp, ncol = 2,
                     rel_widths = c(1,1), labels = c("a", "b"),
                     scale = 0.9, label_size = 20)

ggsave(file.path("output", "figures", "ANOVA_baseline_dynamics.svg"),
       plot = p3_comp,
       dpi = 100,
       units = "in",
       height = 4.5,
       width = 10)

# Calculate a pseudo pvalue from quantiles
res2 <- sum(null_Rsq_comp > obs_Rsq_comp) / p
cat(paste0("P-value from random null test: ", round(res2, 3), "\n"))
