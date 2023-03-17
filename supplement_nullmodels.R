source("path_fix.R")

library(tidyverse)
library(magrittr)
library(rulesoflife)
library(driver)
library(cowplot)
library(fido)

# ------------------------------------------------------------------------------
#
#   Supplemental Figure S4 - null vs. observed a) CLR correlations and b)
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
    D <- nrow(Sigma)
    Sigma <- Sigma[1:(D-1),1:(D-1)]
    if(is.null(D_combos)) {
      D_combos <- ((D-1)^2 - (D-1))/2
      permuted_correlation_phy <- matrix(NA, D_combos, length(fits)*psamples)
    }
    offset <- length(fits)*(i-1) + j
    permuted_correlation_phy[,offset] <- Sigma[upper.tri(Sigma)]
  }
}

# correlation_phy <- NULL
# pdir <- "phy_days90_diet25_scale1"
# fits <- list.files(file.path("output", "model_fits", pdir, "MAP"), full.names = TRUE)
# for(j in 1:length(fits)) {
#   Sigma <- cov2cor(to_clr(readRDS(fits[j]))$Sigma[,,1])
#   D <- nrow(Sigma)
#   Sigma <- Sigma[1:(D-1),1:(D-1)]
#   if(is.null(correlation_phy)) {
#     D_combos <- ((D-1)^2 - (D-1))/2
#     correlation_phy <- matrix(NA, D_combos, length(fits))
#   }
#   correlation_phy[,j] <- Sigma[upper.tri(Sigma)]
# }
# corr_distros <- rbind(corr_distros,
#                       data.frame(x = sample(c(correlation_phy)),
#                                  type = "Phylum",
#                                  scheme = "observed",
#                                  host = unname(sapply(fits, function(x) str_match(x, ".*MAP\\/(.*?)\\.rds")[[2]]))))

data <- load_data(tax_level = "phylum")
rug_phy <- summarize_Sigmas(output_dir = "phy_days90_diet25_scale1")
filtered_pairs <- filter_joint_zeros(data$counts, threshold_and = 0.05, threshold_or = 0.5)
rug <- rug_phy$rug[,filtered_pairs$threshold]

# corr_distros <- data.frame(x = sample(c(permuted_correlation_phy), size = 10000),
corr_distros <- data.frame(x = c(permuted_correlation_phy),
                           type = "Phylum",
                           scheme = "permuted",
                           host = NA)

for(h in 1:nrow(rug)) {
  corr_distros <- rbind(corr_distros,
                        data.frame(x = rug[h,],
                                   type = "Phylum",
                                   scheme = "observed",
                                   host = rug_phy$hosts[h]))
}

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
    D <- nrow(Sigma)
    Sigma <- Sigma[1:(D-1),1:(D-1)]
    if(is.null(D_combos)) {
      D_combos <- ((D-1)^2 - (D-1))/2
      permuted_correlation_fam <- matrix(NA, D_combos, length(fits)*psamples)
    }
    offset <- length(fits)*(i-1) + j
    permuted_correlation_fam[,offset] <- Sigma[upper.tri(Sigma)]
  }
}

# correlation_fam <- NULL
# pdir <- "fam_days90_diet25_scale1"
# fits <- list.files(file.path("output", "model_fits", pdir, "MAP"), full.names = TRUE)
# for(j in 1:length(fits)) {
#   Sigma <- cov2cor(to_clr(readRDS(fits[j]))$Sigma[,,1])
#   D <- nrow(Sigma)
#   Sigma <- Sigma[1:(D-1),1:(D-1)]
#   if(is.null(correlation_fam)) {
#     D_combos <- ((D-1)^2 - (D-1))/2
#     correlation_fam <- matrix(NA, D_combos, length(fits))
#   }
#   correlation_fam[,j] <- Sigma[upper.tri(Sigma)]
# }
# corr_distros <- rbind(corr_distros,
#                       data.frame(x = sample(c(correlation_fam)),
#                                  type = "Family",
#                                  scheme = "observed",
#                                  host = unname(sapply(fits, function(x) str_match(x, ".*MAP\\/(.*?)\\.rds")[[2]]))))

corr_distros <- rbind(corr_distros,
                      # data.frame(x = sample(c(permuted_correlation_fam), size = 10000),
                      data.frame(x = c(permuted_correlation_fam),
                                 type = "Family",
                                 scheme = "permuted",
                                 host = NA))

data <- load_data(tax_level = "family")
rug_fam <- summarize_Sigmas(output_dir = "fam_days90_diet25_scale1")
filtered_pairs <- filter_joint_zeros(data$counts, threshold_and = 0.05, threshold_or = 0.5)
rug <- rug_fam$rug[,filtered_pairs$threshold]

for(h in 1:nrow(rug)) {
  corr_distros <- rbind(corr_distros,
                        data.frame(x = rug[h,],
                                   type = "Family",
                                   scheme = "observed",
                                   host = rug_fam$hosts[h]))
}

# ------------------------------------------------------------------------------
#   ASV
# ------------------------------------------------------------------------------

D_combos <- NULL
permuted_correlation_asv <- NULL
for(i in 1:psamples) {
  cat(paste0("Loading permutation #", i, "...\n"))
  pdir <- paste0("asv_days90_diet25_scale1_scramble-sample-", sprintf("%02d", i))
  fits <- list.files(file.path("output", "model_fits", pdir, "MAP"), full.names = TRUE)
  # full_Y <- matrix(NA, 126, 5534) # hard-code these dimensions for simplicity
  # Y_offset <- 1
  for(j in 1:length(fits)) {
    clr_fit <- to_clr(readRDS(fits[j]))
    # Y <- clr_fit$Y
    # full_Y[,Y_offset:(Y_offset + ncol(Y) - 1)] <- Y
    # Y_offset <- Y_offset + ncol(Y)
    Sigma <- cov2cor(clr_fit$Sigma[,,1])
    D <- nrow(Sigma)
    Sigma <- Sigma[1:(D-1),1:(D-1)]
    if(is.null(D_combos)) {
      D_combos <- ((D-1)^2 - (D-1))/2
      permuted_correlation_asv <- matrix(NA, D_combos, length(fits)*psamples)
    }
    offset <- length(fits)*(i-1) + j
    permuted_correlation_asv[,offset] <- Sigma[upper.tri(Sigma)]
  }
  # filtered_pairs <- filter_joint_zeros(full_Y)
  # Show that permuted pairs don't have the systematic effect whereby increasing
  # frequency of joint zeros is associated with positive correlation. In fact,
  # the ranges for joint zeros and correlation are really small in permuted
  # data!
  # plot(filtered_pairs$frequency_00, rowMeans(permuted_correlation_asv, na.rm = T))
}

# correlation_asv <- NULL
# pdir <- "asv_days90_diet25_scale1"
# fits <- list.files(file.path("output", "model_fits", pdir, "MAP"), full.names = TRUE)
# for(j in 1:length(fits)) {
#   Sigma <- cov2cor(to_clr(readRDS(fits[j]))$Sigma[,,1])
#   D <- nrow(Sigma)
#   Sigma <- Sigma[1:(D-1),1:(D-1)]
#   if(is.null(correlation_asv)) {
#     D_combos <- ((D-1)^2 - (D-1))/2
#     correlation_asv <- matrix(NA, D_combos, length(fits))
#   }
#   correlation_asv[,j] <- Sigma[upper.tri(Sigma)]
# }
# corr_distros <- rbind(corr_distros,
#                       data.frame(x = sample(c(correlation_asv)),
#                                  type = "ASV",
#                                  scheme = "observed",
#                                  host = unname(sapply(fits, function(x) str_match(x, ".*MAP\\/(.*?)\\.rds")[[2]]))))

corr_distros <- rbind(corr_distros,
                      # data.frame(x = sample(c(permuted_correlation_asv), size = 10000),
                      data.frame(x = c(permuted_correlation_asv),
                                 type = "ASV",
                                 scheme = "permuted",
                                 host = NA))

data <- load_data(tax_level = "ASV")
rug_asv <- summarize_Sigmas(output_dir = "asv_days90_diet25_scale1")
filtered_pairs <- filter_joint_zeros(data$counts, threshold_and = 0.05, threshold_or = 0.5)
rug <- rug_asv$rug[,filtered_pairs$threshold]

# This takes ~ 30 sec.
for(h in 1:nrow(rug)) {
  corr_distros <- rbind(corr_distros,
                        data.frame(x = rug[h,],
                                   type = "ASV",
                                   scheme = "observed",
                                   host = rug_asv$hosts[h]))
}

# Get rid of NAs that have resulted from joint zero filtering
corr_distros %<>%
  filter(!is.na(x))

corr_distros$type <- factor(corr_distros$type, levels = c("ASV", "Family", "Phylum"))
levels(corr_distros$type) <- c("ASV", "Family/order/class", "Phylum")
corr_distros$scheme <- factor(corr_distros$scheme, levels = c("observed", "permuted"))

means <- data.frame(type = "Phylum",
                    x = mean(corr_distros %>% filter(scheme == "observed" & type == "Phylum") %>% pull(x)))
means <- rbind(means,
               data.frame(type = "Family/order/class",
                          x = mean(corr_distros %>% filter(scheme == "observed" & type == "Family/order/class") %>% pull(x))))
means <- rbind(means,
               data.frame(type = "ASV",
                          x = mean(corr_distros %>% filter(scheme == "observed" & type == "ASV") %>% pull(x))))

p1 <- ggplot() +
  geom_density(data = corr_distros, mapping = aes(x = x, fill = scheme),
               alpha = 0.5) +
  geom_segment(data = means, mapping = aes(x = x, xend = x, y = 0, yend = -0.1),
               size = 1.5) +
  facet_wrap(. ~ type) +
  scale_fill_manual(values = c("#59AAD7", "#aaaaaa")) +
  theme_bw() +
  labs(fill = "Source",
       x = "correlation",
       y = "density")

# ------------------------------------------------------------------------------
#
#   Null vs. observed universality scores
#
# ------------------------------------------------------------------------------

iterations <- 100

# ------------------------------------------------------------------------------
#   Phylum
# ------------------------------------------------------------------------------

data <- load_data(tax_level = "phylum")
rug_phy <- summarize_Sigmas(output_dir = "phy_days90_diet25_scale1")
filtered_pairs <- filter_joint_zeros(data$counts, threshold_and = 0.05, threshold_or = 0.5)
rug <- rug_phy$rug[,filtered_pairs$threshold]

permuted_scores_phy <- NULL
for(i in 1:iterations) {
  rug_scrambled <- rug
  for(j in 1:nrow(rug)) {
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
                       data.frame(x = apply(rug, 2, calc_universality_score),
                                  type = "Phylum",
                                  scheme = "observed"))

# ------------------------------------------------------------------------------
#   Family
# ------------------------------------------------------------------------------

data <- load_data(tax_level = "family")
rug_fam <- summarize_Sigmas(output_dir = "fam_days90_diet25_scale1")
filtered_pairs <- filter_joint_zeros(data$counts, threshold_and = 0.05, threshold_or = 0.5)
rug <- rug_fam$rug[,filtered_pairs$threshold]

permuted_scores_fam <- NULL
for(i in 1:iterations) {
  rug_scrambled <- rug
  for(j in 1:nrow(rug)) {
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
                       data.frame(x = apply(rug, 2, calc_universality_score),
                                  type = "Family",
                                  scheme = "observed"))

# ------------------------------------------------------------------------------
#   ASV
# ------------------------------------------------------------------------------

data <- load_data(tax_level = "ASV")
rug_asv <- summarize_Sigmas(output_dir = "asv_days90_diet25_scale1")
filtered_pairs <- filter_joint_zeros(data$counts, threshold_and = 0.05, threshold_or = 0.5)
rug <- rug_asv$rug[,filtered_pairs$threshold]

permuted_scores_asv <- NULL
for(i in 1:iterations) {
  if(i %% 10 == 0) {
    cat(paste0("Processing permutation iteration ", i,"\n"))
  }
  rug_scrambled <- rug
  for(j in 1:nrow(rug)) {
    rug_scrambled[j,] <- sample(rug_scrambled[j,])
  }
  scores <- apply(rug_scrambled, 2, calc_universality_score)
  if(is.null(permuted_scores_asv)) {
    permuted_scores_asv <- matrix(NA, length(scores), iterations)
  }
  permuted_scores_asv[,i] <- scores
}

obs_scores <- apply(rug, 2, calc_universality_score)

score_distros <- rbind(score_distros,
                       data.frame(x = c(permuted_scores_asv),
                                  type = "ASV",
                                  scheme = "permuted"))
score_distros <- rbind(score_distros,
                       data.frame(x = obs_scores,
                                  type = "ASV",
                                  scheme = "observed"))

score_distros$type <- factor(score_distros$type, levels = c("ASV", "Family", "Phylum"))
levels(score_distros$type) <- c("ASV", "Family/order/class", "Phylum")
score_distros$scheme <- factor(score_distros$scheme, levels = c("observed", "permuted"))

means <- data.frame(type = "Phylum",
                    x = mean(score_distros %>% filter(scheme == "observed" & type == "Phylum") %>% pull(x)))
means <- rbind(means,
               data.frame(type = "Family/order/class",
                          x = mean(score_distros %>% filter(scheme == "observed" & type == "Family/order/class") %>% pull(x))))
means <- rbind(means,
               data.frame(type = "ASV",
                          x = mean(score_distros %>% filter(scheme == "observed" & type == "ASV") %>% pull(x))))

p2 <- ggplot() +
  geom_density(data = score_distros, mapping = aes(x = x, fill = scheme),
               alpha = 0.5) +
  geom_segment(data = means, mapping = aes(x = x, xend = x, y = 0, yend = -0.6),
               size = 1.5) +
  facet_wrap(. ~ type) +
  scale_fill_manual(values = c("#59AAD7", "#aaaaaa")) +
  theme_bw() +
  xlim(c(-0.05, 0.7)) +
  labs(fill = "Source",
       x = "universality score",
       y = "density")

p <- plot_grid(p1, p2,
               ncol = 1,
               labels = c("A", "B"),
               label_y = 1.02,
               label_size = 16,
               scale = 0.95)

ggsave(file.path("output", "figures", "Figure_2_Supplement_3.png"),
       p,
       dpi = 200,
       units = "in",
       height = 5,
       width = 8,
       bg = "white")

# ------------------------------------------------------------------------------
#   Report stats: How many ASV-level universality scores are lower than
#                 expected by chance?
# ------------------------------------------------------------------------------

if(FALSE) {
  source("thresholds.R")

  bounds <- quantile(score_distros %>% filter(scheme == "permuted") %>% pull(x),
                     probs = c(0.0275, 0.975))

  pw_scores <- apply(rug_asv$rug, 2, function(x) calc_universality_score(x, return_pieces = TRUE))

  # Find pairs will unusually low scores below the 95% interval in the permuted distribution
  lower <- which(apply(pw_scores, 2, function(x) {
    x[1]*x[2] < bounds[1]
  }))

  cat(paste0(round(length(lower) / ncol(rug_asv$rug) * 100, 2),
             "% of pairs have universality scores below the 95% interval of the permuted (null) distribution"))

  subset_correlations <- rug_asv$rug[,lower]

  cat(paste0("The median correlation strength of these pairs is: ",
             round(median(apply(subset_correlations, 2, median)), 3), "\n"))

  cat(paste0("The percent of correlations 'significant' by our correlation threshold is: ",
             round(sum(subset_correlations < thresholds %>% filter(type == "ASV") %>% pull(lower) |
                         subset_correlations > thresholds %>% filter(type == "ASV") %>% pull(upper)) / length(subset_correlations) * 100, 2),
             "\n"))
}

# Optionally visualize the distribution of very-low-universality score pairs
if(FALSE) {
  low_idx <- which(obs_scores < bounds[1])
  subset_rug <- rug_asv$rug[,low_idx]
  scores_pieces <- apply(subset_rug, 2, function(x) {
    calc_universality_score(x, return_pieces = TRUE)
  })

  # Percent agreement across hosts in these very-low scores
  ggplot(data.frame(x = scores_pieces[1,]), aes(x = x)) +
    geom_histogram(color = "white") +
    theme_bw() +
    xlim(c(0.5, 1)) +
    labs(x = "percent agreement in sign across hosts")

  # Percent agreement across hosts in these very-low scores
  subset_rug <- as.data.frame(subset_rug)
  rownames(subset_rug) <- paste0("host_", 1:nrow(subset_rug))
  colnames(subset_rug) <- paste0(1:ncol(subset_rug))
  long_subset <- pivot_longer(subset_rug,
                              everything(),
                              names_to = "pair",
                              values_to = "correlation")
  long_subset$pair <- factor(long_subset$pair)
  # Pull some random pairs to plot
  ggplot(long_subset %>% filter(pair %in% sample(levels(pair), size = 20)),
         aes(x = pair, y = correlation)) +
    geom_violin() +
    theme_bw() +
    ylim(c(-1, 1)) +
    labs(x = "pair index",
         y = "correlation distribution")
}

# ------------------------------------------------------------------------------
#   Report stats: correlation
#
#   Note: these thresholds have been hard-coded into `thresholds.R`
# ------------------------------------------------------------------------------

if(FALSE) {
  corr_distros <- corr_distros[complete.ca]

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
}

# ------------------------------------------------------------------------------
#   Report stats: universality scores
#
#   Note: these thresholds have been hard-coded into `thresholds.R`
# ------------------------------------------------------------------------------

if(FALSE) {
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
}

# ------------------------------------------------------------------------------
#
#   False discovery-adjusted thresholds from BH correction
#
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
#   Correlation
# ------------------------------------------------------------------------------

for(ttype in unique(corr_distros$type)) {
  f <- ecdf(corr_distros %>% filter(type == ttype & scheme == "permuted") %>% pull(x))

  values <- corr_distros %>% filter(type == ttype & scheme == "observed") %>% pull(x)
  pvalues <- sapply(values, function(x) {
    2*min(f(x), 1-f(x))
  })

  pvalues_adj <- p.adjust(pvalues, method = "BH")
  signif_vec <- pvalues_adj < 0.05

  cat(paste0("LEVEL: ", toupper(ttype), "\n"))

  signif_values <- values[signif_vec]
  pos_threshold <- min(signif_values[signif_values > 0])
  neg_threshold <- max(signif_values[signif_values < 0])
  cat(paste0("Negative corr. signif value (cutoff): ", round(neg_threshold, 3), "\n"))
  cat(paste0("Positive corr. signif value (cutoff): ", round(pos_threshold, 3), "\n"))

  cat(paste0("Overall:\n"))
  cat(paste0("\tMedian: ", round(median(values), 3), "\n"))
  cat(paste0("\tPercent positive: ", round(sum(values > 0)/length(values), 3), "\n"))
  cat(paste0("\tPercent negative: ", round(sum(values < 0)/length(values), 3), "\n"))

  prop <- sum(signif_vec) / length(pvalues_adj)
  cat(paste0("Significant after FDR: ", round(prop, 3), " (", sum(signif_vec), " / ", length(pvalues_adj), ")\n"))
  cat(paste0("\tMedian values (signif pairs): ", round(median(abs(values[signif_vec])), 3), "\n"))
  x <- sum(values[signif_vec] > 0)
  n <- sum(signif_vec)
  cat(paste0("\tPercent positive (signif pairs): ", round(x/n, 3), "\n"))
  cat(paste0("\tPercent negative (signif pairs): ", round((n-x)/n, 3), "\n\n"))
  # binom.test(x, n, p = 0.5)
}

# For ASVs, what proportion are significant (and positive) for each host?
temp <- corr_distros %>%
  filter(type == ttype) %>%
  filter(scheme == "observed")
host_list <- sort(unique(temp$host))
prop_signif <- numeric(length(host_list))
prop_positive <- numeric(length(host_list))
for(i in 1:length(host_list)) {
  values <- temp %>%
    filter(host == host_list[i]) %>%
    pull(x)
  pvalues <- sapply(values, function(x) {
    2*min(f(x), 1-f(x))
  })
  pvalues_adj <- p.adjust(pvalues, method = "BH")
  signif_vec <- pvalues_adj < 0.05
  prop_signif[i] <- sum(signif_vec)/length(signif_vec)
  # As a proportion of significant correlations
  # prop_positive[i] <- sum(values[signif_vec] > 0)/sum(signif_vec)
  # As a proportion of all correlations
  prop_positive[i] <- sum(values[signif_vec] > 0)/length(signif_vec)
}

cat(paste0("Min, mean, max proportion significant per host: ",
           round(min(prop_signif), 3),", ",
           round(mean(prop_signif), 3),", ",
           round(max(prop_signif), 3),"\n"))

cat(paste0("Min, mean, max proportion significant and positive per host: ",
           round(min(prop_positive), 3),", ",
           round(mean(prop_positive), 3),", ",
           round(max(prop_positive), 3),"\n"))

# ------------------------------------------------------------------------------
#   "Universality"; I'm using a one-sided test here
# ------------------------------------------------------------------------------

for(ttype in unique(score_distros$type)) {
  f <- ecdf(score_distros %>% filter(type == ttype & scheme == "permuted") %>% pull(x))

  values <- score_distros %>% filter(type == ttype & scheme == "observed") %>% pull(x)
  pvalues <- sapply(values, function(x) {
    1-f(x) # one-sided
    # 2*min(f(x), 1-f(x)) # two-sided
  })

  pvalues_adj <- p.adjust(pvalues, method = "BH")
  signif_vec <- pvalues_adj < 0.05

  min_signif <- min(values[signif_vec])
  cat(paste0("Minimum significant value (cutoff): ", round(min_signif, 5), "\n"))

  prop <- sum(signif_vec) / length(pvalues_adj)
  cat(paste0("Significant after FDR (", ttype, "): ", round(prop, 3), "\n"))
  higher_than_chance <- signif_vec & values > median(values)
  prop <- sum(higher_than_chance) / length(pvalues_adj)
  cat(paste0("\tHIGHER than chance: ", round(prop, 3), "\n"))
  lower_than_chance <- signif_vec & values < median(values)
  prop <- sum(lower_than_chance) / length(pvalues_adj)
  cat(paste0("\tLOWER than chance: ", round(prop, 3), "\n"))

  asv_scores_piecewise <- apply(rug_asv$rug[,lower_than_chance], 2, function(x) calc_universality_score(x, return_pieces = TRUE))

  # hist(rug_asv$rug[,higher_than_chance])
}

