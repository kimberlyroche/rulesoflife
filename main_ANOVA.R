source("path_fix.R")

library(tidyverse)
library(rulesoflife)
library(fido)
library(driver)
library(shapes)
library(frechet)
library(vegan) # for Mantel test
library(pedtools)
library(kinship2)
library(RColorBrewer)
library(cowplot)

# ------------------------------------------------------------------------------
#
#   R^2 association testing
#
# ------------------------------------------------------------------------------

# Specify model output directory
output_dir <- "asv_days90_diet25_scale1"
p <- 999 # number of permutations

# ------------------------------------------------------------------------------
#   Load data (MAP estimates), compute dynamics distances
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
  d <- dist4cov(Sigmas[[a]] + diag(D)*1e-06,
                Sigmas[[b]] + diag(D)*1e-06)$dist
  dynamics_distances[a,b] <- d
  dynamics_distances[b,a] <- d
}

dynamics_dist_vec <- dynamics_distances[lower.tri(dynamics_distances)]

# ------------------------------------------------------------------------------
#
#   R^2: Pedigree associations
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
# p1_ped <- ggplot(data.frame(x = scale(ped_vec), y = scale(dynamics_dist_vec)),
p1_ped <- ggplot(data.frame(x = ped_vec, y = dynamics_dist_vec),
                 aes(x = x, y = y)) +
  # geom_smooth(color = "gray", alpha = 0.5) +
  geom_point(size = 3, shape = 21, fill = "#999999") +
  theme_bw() +
  labs(x = "kinship dissimilarity (1 - % relatedness)",
       y = "host-host dynamics distance")

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

null_Rsq_ped <- sapply(1:p, function(pp) draw_random_Rsq_ped(56))

# Calculate a pseudo pvalue from quantiles
res2 <- sum(null_Rsq_ped > obs_Rsq_ped) / p
cat(paste0("Observed R^2: ", round(obs_Rsq_ped, 3), "\n"))
cat(paste0("P-value from random null test: ", round(res2, 3), "\n"))

p2_ped <- ggplot() +
  geom_density(data = data.frame(x = null_Rsq_ped),
               mapping = aes(x = x)) +
  geom_point(data = data.frame(x = obs_Rsq_ped, y = 1),
             mapping = aes(x = x, y = y),
             size = 4,
             shape = 21,
             fill = "#d99e57") +
  theme_bw() +
  labs(x = "R-squared (random kinship dissimilarity x dynamics distance)",
       y = "density")

p3_ped <- plot_grid(p1_ped, p2_ped, ncol = 2,
                    rel_widths = c(1,1), labels = c("A", "B"),
                    scale = 0.9,
                    label_size = 18)

# ------------------------------------------------------------------------------
#
#   R^2: Baseline composition associations
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
p1_comp <- ggplot(data.frame(x = baseline_dist_vec, y = dynamics_dist_vec),
                  aes(x = x, y = y)) +
  # geom_smooth(color = "gray", alpha = 0.5) +
  geom_point(size = 3, shape = 21, fill = "#999999") +
  theme_bw() +
  labs(x = "Aitchison distance over baselines",
       y = "host-host dynamics distance")

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

# Calculate a pseudo pvalue from quantiles
res2 <- sum(null_Rsq_comp > obs_Rsq_comp) / p
cat(paste0("Observed R^2: ", round(obs_Rsq_comp, 3), "\n"))
cat(paste0("P-value from random null test: ", round(res2, 3), "\n"))

p2_comp <- ggplot() +
  geom_density(data = data.frame(x = null_Rsq_comp),
               mapping = aes(x = x)) +
  geom_point(data = data.frame(x = obs_Rsq_comp, y = 1),
             mapping = aes(x = x, y = y),
             size = 4,
             shape = 21,
             fill = "#d99e57") +
  theme_bw() +
  labs(x = "R-squared (random baseline x dynamics distance)",
       y = "density")

p3_comp <- plot_grid(p1_comp, p2_comp, ncol = 2,
                    rel_widths = c(1,1), labels = c("C", "D"),
                    scale = 0.9,
                    label_size = 18)

# p_all <- plot_grid(p3_ped, p3_comp, ncol = 1)

p_all <- plot_grid(p1_ped, p1_comp,
                   ncol = 2,
                   rel_widths = c(1, 1),
                   labels = c("A", "B"),
                   label_size = 18,
                   label_y = 1.00,
                   label_x = -0.01,
                   scale = 0.925)

ggsave(file.path("output", "figures", "anova.svg"),
       plot = p_all,
       dpi = 100,
       units = "in",
       height = 4.5,
       width = 10)

# ------------------------------------------------------------------------------
#
#   R^2: Lifespan
#
# ------------------------------------------------------------------------------

females_only <- FALSE

# Load lifespan data
lifespan <- read.delim(file.path("input", "lifespanForKim_4Dec2021.csv"), sep = ",")
if(females_only) {
  lifespan <- lifespan %>%
    filter(sex == "F")
}
host_list <- lifespan %>% pull(sname)

host_include <- names(Sigmas) %in% host_list
use_Sigmas <- Sigmas[host_include]

# Load pedigree
mapping <- data.frame(sname = names(use_Sigmas)) %>%
  left_join(data.frame(sname = host_list, order = 1:length(host_list)), by = "sname")

lifespan <- lifespan$age_at_death_censor[mapping$order]
lifespan_dist <- as.matrix(dist(lifespan))
lifespan_vec <- lifespan_dist[lower.tri(lifespan_dist)]

# Subset dynamics distances to these hosts
use_dd <- dynamics_distances[host_include,host_include]
use_dd_vec <- use_dd[lower.tri(use_dd)]

obs_Rsq_lifespan <- c(cor(scale(lifespan_vec), scale(use_dd_vec))**2)

# R^2 plot
# p1_ped <- ggplot(data.frame(x = scale(ped_vec), y = scale(dynamics_dist_vec)),
p1_life <- ggplot(data.frame(x = lifespan_vec, y = use_dd_vec),
                  aes(x = x, y = y)) +
  # geom_smooth(color = "gray", alpha = 0.5) +
  geom_point(size = 3, shape = 21, fill = "#999999") +
  theme_bw() +
  labs(x = paste0("lifespance distance", ifelse(females_only, " (females only)", "")),
       y = "host-host dynamics distance")

# ------------------------------------------------------------------------------
#   Mantel test
#
#   This answers the question: How much variation is explained by X on distances
#   in dynamics *in a population this related*?
# ------------------------------------------------------------------------------

res1_life <- mantel(as.dist(lifespan_dist), as.dist(use_dd))
cat(paste0("P-value from Mantel test (lifespan): ", round(res1_life$signif, 3), "\n"))
cat(paste0("Observed R^2: ", round(obs_Rsq_lifespan, 3), "\n"))

# ------------------------------------------------------------------------------
#   Null distribution - lifespan
#
#   Note: I think permutations of observed lifespan will have the same problem
#     as the Mantel test!
#   Essentially we would need an independent model of lifespan here.
#   But luckily, this isn't remotely significant by Mantel.
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
#
#   Pseudo-ANOVA
#
# ------------------------------------------------------------------------------

# Plotting function
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
  ggsave(file.path(output_dir, paste0(label, ".svg")),
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

# Pull the already-parsed MAP estimates of dynamics from this object, calculated
# by `analysis_compute_Frechets.R`
F1 <- readRDS(file.path("output", "Frechet_1_corr.rds"))

# Calculate the mean using the `frechet` package
F1$mean <- CovFMean(F1$Sigmas)$Mout[[1]]

D <- dim(F1$Sigmas)[1]
N <- dim(F1$Sigmas)[3]

# ------------------------------------------------------------------------------
#   Pull host social group / sex metadata
# ------------------------------------------------------------------------------

host_groups <- get_host_social_groups(hosts)
groups <- unique(host_groups$grp)

metadata <- load_data()$metadata
host_sex <- metadata %>%
  dplyr::select(sname, sex) %>%
  distinct() %>%
  filter(sname %in% hosts) %>%
  arrange(sname)
sexes <- unique(host_sex$sex)

# Load group means
group_means <- list()
for(g in 1:length(groups)) {
  group_means[[g]] <- CovFMean(F1$Sigmas[,,host_groups$grp == groups[g]])$Mout[[1]]
}

# Load group means
sex_means <- list()
for(s in 1:length(sexes)) {
  sex_means[[s]] <- CovFMean(F1$Sigmas[,,host_sex$sex == sexes[s]])$Mout[[1]]
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
  d <- dist4cov(group_means[[i]], F1$mean)$dist**2
  between_sum <- between_sum + (n_i * d)
}
numerator <- between_sum / (K-1)

within_sum <- 0
for(i in 1:K) {
  grp <- groups[i]
  idx_i <- which(host_groups$grp == grp)
  for(j in idx_i) {
    d <- dist4cov(F1$Sigmas[,,j], group_means[[i]])$dist**2
    within_sum <- within_sum + d
  }
}
denominator <- within_sum / (N-K)

ratio <- numerator / denominator

cat(paste0("P-value for pseudo-ANOVA on GROUP (F=",
           round(ratio, 3),
           "): ",
           round(1 - pf(ratio, K-1, N-K), 3),
           "\n"))

# ------------------------------------------------------------------------------
#   ANOVA on sex
# ------------------------------------------------------------------------------

K <- length(sexes)
N <- nrow(host_sex)
between_sum <- 0
for(i in 1:K) {
  sex <- sexes[i]
  n_i <- sum(host_sex == sex)
  d <- dist4cov(sex_means[[i]], F1$mean)$dist**2
  between_sum <- between_sum + (n_i * d)
}
numerator <- between_sum / (K-1)

within_sum <- 0
for(i in 1:K) {
  sex <- sexes[i]
  idx_i <- which(host_sex$sex == sex)
  for(j in idx_i) {
    d <- dist4cov(F1$Sigmas[,,j], sex_means[[i]])$dist**2
    within_sum <- within_sum + d
  }
}
denominator <- within_sum / (N-K)

ratio <- numerator / denominator

cat(paste0("P-value for pseudo-ANOVA on SEX (F=",
           round(ratio, 3),
           "): ",
           round(1 - pf(ratio, K-1, N-K), 3),
           "\n"))
