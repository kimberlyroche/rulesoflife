source("path_fix.R")

library(tidyverse)
library(rulesoflife)
library(driver)
library(cowplot)
library(fido)

# ------------------------------------------------------------------------------
#
#   Supplemental Figure S3 - null vs. observed a) CLR correlations and b)
#                            universality scores
#
# ------------------------------------------------------------------------------

psamples <- 10

# ------------------------------------------------------------------------------
#   Phylum
# ------------------------------------------------------------------------------

D_combos <- NULL
permuted_correlation_phy <- NULL
for(i in 1:psamples) {
  cat(paste0("Loading permutation #", i, "...\n"))
  pdir <- paste0("phy_days90_diet25_scale1_scramble-sample-", sprintf("%02d", i))
  fits <- list.files(file.path("output", "model_fits", pdir, "MAP"), full.names = TRUE)
  for(j in 1:length(fits)) {
    Sigma <- cov2cor(to_clr(readRDS(fits[j]))$Sigma[,,1])
    if(is.null(D_combos)) {
      D <- nrow(Sigma)
      D_combos <- (D^2 - D)/2
      permuted_correlation_phy <- matrix(NA, D_combos, length(fits)*psamples)
    }
    offset <- length(fits)*(i-1) + j
    permuted_correlation_phy[,offset] <- Sigma[upper.tri(Sigma)]
  }
}

correlation_phy <- NULL
pdir <- "phy_days90_diet25_scale1"
fits <- list.files(file.path("output", "model_fits", pdir, "MAP"), full.names = TRUE)
for(j in 1:length(fits)) {
  Sigma <- cov2cor(to_clr(readRDS(fits[j]))$Sigma[,,1])
  if(is.null(correlation_phy)) {
    D <- nrow(Sigma)
    D_combos <- (D^2 - D)/2
    correlation_phy <- matrix(NA, D_combos, length(fits))
  }
  correlation_phy[,j] <- Sigma[upper.tri(Sigma)]
}

# corr_distros <- data.frame(x = sample(c(permuted_correlation_phy), size = 10000),
corr_distros <- data.frame(x = c(permuted_correlation_phy),
                           type = "Phylum",
                           scheme = "permuted")
corr_distros <- rbind(corr_distros,
                      data.frame(x = sample(c(correlation_phy)),
                                 type = "Phylum",
                                 scheme = "observed"))

# ------------------------------------------------------------------------------
#   Family
# ------------------------------------------------------------------------------

D_combos <- NULL
permuted_correlation_fam <- NULL
for(i in 1:psamples) {
  cat(paste0("Loading permutation #", i, "...\n"))
  pdir <- paste0("fam_days90_diet25_scale1_scramble-sample-", sprintf("%02d", i))
  fits <- list.files(file.path("output", "model_fits", pdir, "MAP"), full.names = TRUE)
  for(j in 1:length(fits)) {
    Sigma <- cov2cor(to_clr(readRDS(fits[j]))$Sigma[,,1])
    if(is.null(D_combos)) {
      D <- nrow(Sigma)
      D_combos <- (D^2 - D)/2
      permuted_correlation_fam <- matrix(NA, D_combos, length(fits)*psamples)
    }
    offset <- length(fits)*(i-1) + j
    permuted_correlation_fam[,offset] <- Sigma[upper.tri(Sigma)]
  }
}

correlation_fam <- NULL
pdir <- "fam_days90_diet25_scale1"
fits <- list.files(file.path("output", "model_fits", pdir, "MAP"), full.names = TRUE)
for(j in 1:length(fits)) {
  Sigma <- cov2cor(to_clr(readRDS(fits[j]))$Sigma[,,1])
  if(is.null(correlation_fam)) {
    D <- nrow(Sigma)
    D_combos <- (D^2 - D)/2
    correlation_fam <- matrix(NA, D_combos, length(fits))
  }
  correlation_fam[,j] <- Sigma[upper.tri(Sigma)]
}

corr_distros <- rbind(corr_distros,
                      # data.frame(x = sample(c(permuted_correlation_fam), size = 10000),
                      data.frame(x = c(permuted_correlation_fam),
                                 type = "Family",
                                 scheme = "permuted"))
corr_distros <- rbind(corr_distros,
                      data.frame(x = sample(c(correlation_fam)),
                                 type = "Family",
                                 scheme = "observed"))

# ------------------------------------------------------------------------------
#   ASV
# ------------------------------------------------------------------------------

D_combos <- NULL
permuted_correlation_asv <- NULL
for(i in 1:psamples) {
  cat(paste0("Loading permutation #", i, "...\n"))
  pdir <- paste0("asv_days90_diet25_scale1_scramble-sample-", sprintf("%02d", i))
  fits <- list.files(file.path("output", "model_fits", pdir, "MAP"), full.names = TRUE)
  for(j in 1:length(fits)) {
    Sigma <- cov2cor(to_clr(readRDS(fits[j]))$Sigma[,,1])
    if(is.null(D_combos)) {
      D <- nrow(Sigma)
      D_combos <- (D^2 - D)/2
      permuted_correlation_asv <- matrix(NA, D_combos, length(fits)*psamples)
    }
    offset <- length(fits)*(i-1) + j
    permuted_correlation_asv[,offset] <- Sigma[upper.tri(Sigma)]
  }
}

correlation_asv <- NULL
pdir <- "asv_days90_diet25_scale1"
fits <- list.files(file.path("output", "model_fits", pdir, "MAP"), full.names = TRUE)
for(j in 1:length(fits)) {
  Sigma <- cov2cor(to_clr(readRDS(fits[j]))$Sigma[,,1])
  if(is.null(correlation_asv)) {
    D <- nrow(Sigma)
    D_combos <- (D^2 - D)/2
    correlation_asv <- matrix(NA, D_combos, length(fits))
  }
  correlation_asv[,j] <- Sigma[upper.tri(Sigma)]
}

corr_distros <- rbind(corr_distros,
                      # data.frame(x = sample(c(permuted_correlation_asv), size = 10000),
                      data.frame(x = c(permuted_correlation_asv),
                                 type = "ASV",
                                 scheme = "permuted"))
corr_distros <- rbind(corr_distros,
                      data.frame(x = sample(c(correlation_asv)),
                                 type = "ASV",
                                 scheme = "observed"))

corr_distros$type <- factor(corr_distros$type, levels = c("Phylum", "Family", "ASV"))
corr_distros$scheme <- factor(corr_distros$scheme, levels = c("observed", "permuted"))

p1 <- ggplot(corr_distros, aes(x = x, fill = scheme)) +
  geom_density(alpha = 0.5) +
  facet_wrap(. ~ type) +
  scale_fill_manual(values = c("#59AAD7", "#aaaaaa")) +
  theme_bw() +
  labs(fill = "Source",
       x = "correlation")

# ------------------------------------------------------------------------------
#
#   Null vs. observed universality scores
#
# ------------------------------------------------------------------------------

iterations <- 100

# ------------------------------------------------------------------------------
#   Phylum
# ------------------------------------------------------------------------------

rug_phy <- summarize_Sigmas(output_dir = "phy_days90_diet25_scale1")
permuted_scores_phy <- NULL
for(i in 1:iterations) {
  rug_scrambled <- rug_phy$rug
  for(j in 1:nrow(rug_phy$rug)) {
    rug_scrambled[j,] <- sample(rug_scrambled[j,])
  }
  scores <- apply(rug_scrambled, 2, calc_universality_score)
  if(is.null(permuted_scores_phy)) {
    permuted_scores_phy <- matrix(NA, length(scores), iterations)
  }
  permuted_scores_phy[,i] <- scores
}

score_distros <- data.frame(x = c(permuted_scores_phy),
                            type = "Phylum",
                            scheme = "permuted")
score_distros <- rbind(score_distros,
                       data.frame(x = apply(rug_phy$rug, 2, calc_universality_score),
                                  type = "Phylum",
                                  scheme = "observed"))

# ------------------------------------------------------------------------------
#   Family
# ------------------------------------------------------------------------------

rug_fam <- summarize_Sigmas(output_dir = "fam_days90_diet25_scale1")
permuted_scores_fam <- NULL
for(i in 1:iterations) {
  rug_scrambled <- rug_fam$rug
  for(j in 1:nrow(rug_fam$rug)) {
    rug_scrambled[j,] <- sample(rug_scrambled[j,])
  }
  scores <- apply(rug_scrambled, 2, calc_universality_score)
  if(is.null(permuted_scores_fam)) {
    permuted_scores_fam <- matrix(NA, length(scores), iterations)
  }
  permuted_scores_fam[,i] <- scores
}

score_distros <- rbind(score_distros,
                       data.frame(x = c(permuted_scores_fam),
                                  type = "Family",
                                  scheme = "permuted"))
score_distros <- rbind(score_distros,
                       data.frame(x = apply(rug_fam$rug, 2, calc_universality_score),
                                  type = "Family",
                                  scheme = "observed"))

# ------------------------------------------------------------------------------
#   ASV
# ------------------------------------------------------------------------------

rug_asv <- summarize_Sigmas(output_dir = "asv_days90_diet25_scale1")
permuted_scores_asv <- NULL
for(i in 1:iterations) {
  rug_scrambled <- rug_asv$rug
  for(j in 1:nrow(rug_asv$rug)) {
    rug_scrambled[j,] <- sample(rug_scrambled[j,])
  }
  scores <- apply(rug_scrambled, 2, calc_universality_score)
  if(is.null(permuted_scores_asv)) {
    permuted_scores_asv <- matrix(NA, length(scores), iterations)
  }
  permuted_scores_asv[,i] <- scores
}

score_distros <- rbind(score_distros,
                       data.frame(x = c(permuted_scores_asv),
                                  type = "ASV",
                                  scheme = "permuted"))
score_distros <- rbind(score_distros,
                       data.frame(x = apply(rug_asv$rug, 2, calc_universality_score),
                                  type = "ASV",
                                  scheme = "observed"))

score_distros$type <- factor(score_distros$type, levels = c("Phylum", "Family", "ASV"))
score_distros$scheme <- factor(score_distros$scheme, levels = c("observed", "permuted"))

p2 <- ggplot(score_distros, aes(x = x, fill = scheme)) +
  geom_density(alpha = 0.5) +
  facet_wrap(. ~ type) +
  scale_fill_manual(values = c("#59AAD7", "#aaaaaa")) +
  theme_bw() +
  xlim(c(0, 0.7)) +
  labs(fill = "Source",
       x = "universality score")

p <- plot_grid(p1, p2, ncol = 1, labels = c("a", "b"), label_y = 1.05, label_size = 16)

ggsave("output/figures/SF2.svg",
       p,
       dpi = 100,
       units = "in",
       height = 5,
       width = 8)

# ------------------------------------------------------------------------------
#   Report stats: correlation
# ------------------------------------------------------------------------------

# Define cutoffs at <0.025 percentile and >97.5 percentile
phy_thresholds <- quantile(corr_distros %>% filter(type == "Phylum" & scheme == "permuted") %>% pull(x),
                           probs = c(0.025, 0.975))
fam_thresholds <- quantile(corr_distros %>% filter(type == "Family" & scheme == "permuted") %>% pull(x),
                           probs = c(0.025, 0.975))
asv_thresholds <- quantile(corr_distros %>% filter(type == "ASV" & scheme == "permuted") %>% pull(x),
                           probs = c(0.025, 0.975))

cat(paste0("Lower/upper 95% threshold for PHYLUM: ",
           round(phy_thresholds[1], 3),
           " / ",
           round(phy_thresholds[2], 3),
           "\n"))

phy_signif <- rug_phy$rug < phy_thresholds[1] | rug_phy$rug > phy_thresholds[2]
phy_signif_signs <- sign(rug_phy$rug[phy_signif])
cat(paste0("Proportion observed values exceeding these thresholds for PHYLUM (all / positive / negative): ",
           round(sum(phy_signif)/length(c(rug_phy$rug)), 3),
           "\n"))
cat(paste0("\tProportion exceeding threshold and positive: ",
           round(sum(phy_signif_signs > 0)/length(c(phy_signif_signs)), 3),
           "\n"))
cat(paste0("\tProportion exceeding threshold and negative: ",
           round(sum(phy_signif_signs < 0)/length(c(phy_signif_signs)), 3),
           "\n"))

cat(paste0("Mean / median / range of correlations for PHYLUM: ",
           round(mean(c(rug_phy$rug)), 3), " / ",
           round(median(c(rug_phy$rug)), 3), " / ",
           "(", round(min(c(rug_phy$rug)), 3), ", ", round(max(c(rug_phy$rug)), 3), ")",
           "\n"))

cat(paste0("Lower/upper 95% threshold for FAMILY: ",
           round(fam_thresholds[1], 3),
           " / ",
           round(fam_thresholds[2], 3),
           "\n"))

fam_signif <- rug_fam$rug < fam_thresholds[1] | rug_fam$rug > fam_thresholds[2]
fam_signif_signs <- sign(rug_fam$rug[fam_signif])
cat(paste0("Proportion observed values exceeding these thresholds for FAMILY (all / positive / negative): ",
           round(sum(fam_signif)/length(c(rug_fam$rug)), 3),
           "\n"))
cat(paste0("\tProportion exceeding threshold and positive: ",
           round(sum(fam_signif_signs > 0)/length(c(fam_signif_signs)), 3),
           "\n"))
cat(paste0("\tProportion exceeding threshold and negative: ",
           round(sum(fam_signif_signs < 0)/length(c(fam_signif_signs)), 3),
           "\n"))

cat(paste0("Mean / median / range of correlations for FAMILY: ",
           round(mean(c(rug_fam$rug)), 3), " / ",
           round(median(c(rug_fam$rug)), 3), " / ",
           "(", round(min(c(rug_fam$rug)), 3), ", ", round(max(c(rug_fam$rug)), 3), ")",
           "\n"))

cat(paste0("Lower/upper 95% threshold for ASV: ",
           round(asv_thresholds[1], 3),
           " / ",
           round(asv_thresholds[2], 3),
           "\n"))

asv_signif <- rug_asv$rug < asv_thresholds[1] | rug_asv$rug > asv_thresholds[2]
asv_signif_signs <- sign(rug_asv$rug[asv_signif])
cat(paste0("Proportion observed values exceeding these thresholds for ASV (all / positive / negative): ",
           round(sum(asv_signif)/length(c(rug_asv$rug)), 3),
           "\n"))
cat(paste0("\tProportion exceeding threshold and positive: ",
           round(sum(asv_signif_signs > 0)/length(c(asv_signif_signs)), 3),
           "\n"))
cat(paste0("\tProportion exceeding threshold and negative: ",
           round(sum(asv_signif_signs < 0)/length(c(asv_signif_signs)), 3),
           "\n"))

cat(paste0("Mean / median / range of correlations for ASV: ",
           round(mean(c(rug_asv$rug)), 3), " / ",
           round(median(c(rug_asv$rug)), 3), " / ",
           "(", round(min(c(rug_asv$rug)), 3), ", ", round(max(c(rug_asv$rug)), 3), ")",
           "\n"))

# ------------------------------------------------------------------------------
#   Report stats: universality scores
# ------------------------------------------------------------------------------

obs_scores_phy <- apply(rug_phy$rug, 2, calc_universality_score)
obs_scores_fam <- apply(rug_fam$rug, 2, calc_universality_score)
obs_scores_asv <- apply(rug_asv$rug, 2, calc_universality_score)

# Define cutoffs at <0.025 percentile and >97.5 percentile
phy_thresholds <- quantile(score_distros %>% filter(type == "Phylum" & scheme == "permuted") %>% pull(x),
                           probs = c(0.95))
fam_thresholds <- quantile(score_distros %>% filter(type == "Family" & scheme == "permuted") %>% pull(x),
                           probs = c(0.95))
asv_thresholds <- quantile(score_distros %>% filter(type == "ASV" & scheme == "permuted") %>% pull(x),
                           probs = c(0.95))

cat(paste0("Upper 95% threshold for PHYLUM: ",
           round(phy_thresholds[1], 3),
           "\n"))

cat(paste0("Observed values exceeding these thresholds for PHYLUM: ",
           round(sum(obs_scores_phy > phy_thresholds[1])/length(c(obs_scores_phy)), 3),
           "\n"))

cat(paste0("Upper 95% threshold for FAMILY: ",
           round(fam_thresholds[1], 3),
           "\n"))

cat(paste0("Observed values exceeding these thresholds for FAMILY: ",
           round(sum(obs_scores_fam > fam_thresholds[1])/length(c(obs_scores_fam)), 3),
           "\n"))

cat(paste0("Upper 95% threshold for ASV: ",
           round(asv_thresholds[1], 3),
           "\n"))

cat(paste0("Observed values exceeding these thresholds for ASV: ",
           round(sum(obs_scores_asv > asv_thresholds[1])/length(c(obs_scores_asv)), 3),
           "\n"))
