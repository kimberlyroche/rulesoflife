source("path_fix.R")

library(rulesoflife)
library(driver)
library(tidyverse)
library(cowplot)
library(frechet)
library(ggridges)
library(shapes)

# ------------------------------------------------------------------------------
#   Functions
# ------------------------------------------------------------------------------

plot_ac <- function(scores, title = NULL) {
  scores2 <- scores %>%
    filter(complete.cases(.)) %>%
    group_by(diff) %>%
    mutate(bottom90 = quantile(dist, probs = c(0.05)),
           top90 = quantile(dist, probs = c(0.95)),
           bottom50 = quantile(dist, probs = c(0.25)),
           top50 = quantile(dist, probs = c(0.75)),
           mean = mean(dist))
  ggplot() +
    geom_ribbon(data = scores2,
                mapping = aes(x = diff, ymin = bottom90, ymax = top90),
                fill = "#dddddd") +
    geom_ribbon(data = scores2,
                mapping = aes(x = diff, ymin = bottom50, ymax = top50),
                fill = "#aaaaaa") +
    geom_line(data = scores2,
              mapping = aes(x = diff, y = mean),
              color = "#444444",
              size = 1) +
    xlim(c(0, 208)) + # 4 years in weeks
    theme_bw() +
    labs(x = "weeks", y = "correlation", title = title)
}

sin_f <- function(j, x) {
  sin(j*7 + x/(365/(2*pi)))
}

# ------------------------------------------------------------------------------
#   Parse data
# ------------------------------------------------------------------------------

data <- load_data(tax_level = "ASV")
md <- data$metadata
counts <- data$counts
tax <- data$taxonomy

md_full <- md %>%
  left_join(data.frame(readRDS("input/ps_w_covs.RDS")@sam_data) %>%
              dplyr::select(c(sample_id, rain_monthly, tempmax_monthly)), by = "sample_id")

# We won't need these until we need to fit the models from scratch
clr_counts <- clr_array(counts + 0.5, parts = 1)
baseline <- min(md_full$collection_date)
days <- unname(sapply(md_full$collection_date, function(x) difftime(x, baseline, units = "days")))

# ------------------------------------------------------------------------------
#   Fit AR model WITHOUT periodicity
# ------------------------------------------------------------------------------

save_fn <- file.path("output", "AR_fits.rds")
if(file.exists(save_fn)) {
  fit_obj <- readRDS(save_fn)
  offsets <- fit_obj$offsets
  fits_ar <- fit_obj$fits
  clr_ar <- matrix(NA, nrow(clr_counts), ncol(clr_counts))
  for(i in 1:nrow(clr_counts)) {
    clr_ar[i,] <- c(fits_ar[[i]]$residuals)
  }
} else {
  # This takes up to 30 minutes
  fits_ar <- list()
  clr_ar <- matrix(NA, nrow(clr_counts), ncol(clr_counts))
  for(i in 1:nrow(clr_counts)) {
    cat(paste0("Taxon ", i, " / ", nrow(clr_counts), "\n"))
    x.quant <- model.matrix(~ factor(md_full$sname))
    res_aper <- arima(x = clr_counts[i,],
                      xreg = x.quant[,-1], # remove intercept column
                      order = c(1, 0, 0))
    fits_ar[[i]] <- res_aper
    clr_ar[i,] <- unname(unlist(res_aper$residuals))
  }
  saveRDS(list(fits = fits_ar), save_fn)
}

# ------------------------------------------------------------------------------
#   Fit AR model WITH periodicity
# ------------------------------------------------------------------------------

save_fn <- file.path("output", "AR_fits_per.rds")
if(file.exists(save_fn)) {
  fit_obj <- readRDS(save_fn)
  offsets <- fit_obj$offsets
  fits_aper <- fit_obj$fits
  clr_aper <- matrix(NA, nrow(clr_counts), ncol(clr_counts))
  for(i in 1:nrow(clr_counts)) {
    clr_aper[i,] <- c(fits_aper[[i]]$residuals)
  }
} else {
  # Do an initial run to find the best fit to a fixed-periodicity trend in each taxon
  # This takes < 5 min.
  offsets <- numeric(nrow(clr_counts))
  for(i in 1:nrow(clr_counts)) {
    cat(paste0("Taxon ", i, " / ", nrow(clr_counts), "\n"))
    per_bestfit <- NULL
    per_offset <- NULL
    per_sse <- Inf
    for(j in 0:51) {
      res_aper <- lm(clr_counts[i,] ~ md_full$sname + sin_f(j, days))
      if(per_sse > sum(res_aper$residuals**2)) {
        per_bestfit <- res_aper
        per_offset <- j
        per_sse <- sum(res_aper$residuals**2)
      }
    }
    offsets[i] <- per_offset
  }

  # Now fit a fuller model with AR component
  # This takes 30 minutes!
  fits_aper <- list()
  clr_aper <- matrix(NA, nrow(clr_counts), ncol(clr_counts))
  for(i in 1:nrow(clr_counts)) {
    cat(paste0("Taxon ", i, " / ", nrow(clr_counts), "\n"))
    x.quant <- model.matrix(~ factor(md_full$sname) + sin_f(offsets[i], days))
    res_aper <- arima(x = clr_counts[i,],
                      xreg = x.quant[,-1], # remove intercept column
                      order = c(1, 0, 0))
    fits_aper[[i]] <- res_aper
    clr_aper[i,] <- unname(unlist(res_aper$residuals))
  }
  saveRDS(list(offsets = offsets, fits = fits_aper), save_fn)
}

# ------------------------------------------------------------------------------
#   Generate autocorrelation plots
# ------------------------------------------------------------------------------

# Calculate autocorrelation
# This takes ~3 minutes
md_full$count <- 1:nrow(md_full)
scores <- NULL
scores_ar <- NULL
scores_aper <- NULL
hosts <- unique(md_full$sname)
for(host in hosts) {
  cat(paste0("Host is: ", host, "\n"))
  hdates <- md_full %>%
    filter(sname == host) %>%
    dplyr::select(count, collection_date)
  hcombos <- combn(1:length(hdates$collection_date), m = 2)
  hdays <- unname(sapply(1:ncol(hcombos), function(x) difftime(hdates$collection_date[hcombos[2,x]],
                                                               hdates$collection_date[hcombos[1,x]], units = "days")))
  # Round to monthly
  hdays <- round(hdays/7)
  for(i in 0:max(hdays)) {
    d_idx <- which(hdays == i)
    j <- hdates$count[hcombos[1,d_idx[1]]]
    k <- hdates$count[hcombos[2,d_idx[1]]]
    scores <- rbind(scores,
                    data.frame(host = host,
                               diff = i,
                               dist = cor(clr_counts[,j], clr_counts[,k])))
    scores_ar <- rbind(scores_ar,
                       data.frame(host = host,
                                  diff = i,
                                  dist = cor(clr_ar[,j], clr_ar[,k])))
    scores_aper <- rbind(scores_aper,
                         data.frame(host = host,
                                    diff = i,
                                    dist = cor(clr_aper[,j], clr_aper[,k])))
  }
}

# ------------------------------------------------------------------------------
#   Generate autocorrelation plots
# ------------------------------------------------------------------------------

p1 <- plot_ac(scores_ar)
p2 <- plot_ac(scores_aper)

# ------------------------------------------------------------------------------
#   Compare model estimates
# ------------------------------------------------------------------------------

# Pull ASV-ASV consensus correlations derived from fido models
# F1 <- readRDS(file.path("output", "Frechet_1_corr.rds"))
# F1$mean <- CovFMean(F1$Sigmas)$Mout[[1]]

# Sigmas <- pull_Sigmas("asv_days90_diet25_scale1")
rug <- summarize_Sigmas(output_dir = "asv_days90_diet25_scale1")$rug
data <- load_data(tax_level = "ASV")
filtered_pairs <- filter_joint_zeros(data$counts, threshold_and = 0.05, threshold_or = 0.5)
y <- c(colMeans(rug[,filtered_pairs$threshold]))

# Estimate correlation from aperiodic/season-less joint model
# F1 <- estcov(Sigmas, method = "Euclidean")
# D <- dim(F1$mean)[1]
cov_after <- cov(t(clr_aper))
cor_after <- cov2cor(cov_after)
cor_after <- cor_after[1:D,1:D] # omit "other"

x <- cor_after[lower.tri(cor_after, diag = FALSE)]
x <- x[filtered_pairs$threshold]

# x <- c(cor_after[upper.tri(cor_after, diag = FALSE)])
# y <- c(F1$mean[upper.tri(F1$mean, diag = FALSE)])

# Plot correlation of CLR correlation estimates
p3 <- ggplot(data.frame(x = x,
                        y = y),
             aes(x = x, y = y)) +
  geom_point(size = 2, shape = 21, fill = "#bbbbbb") +
  theme_bw() +
  labs(x = "CLR correlation (original)",
       y = "CLR correlation (minus seasonal effect)")

p <- plot_grid(p1, NULL, p2, NULL, p3,
               ncol = 5,
               rel_widths = c(1, 0.05, 1, 0.05, 1),
               labels = c("A", "", "B", "", "C"),
               label_size = 18,
               label_y = 1.01,
               label_x = -0.015)

ggsave(file.path("output", "figures", "lmer_results_alt.png"),
       p,
       dpi = 200,
       units = "in",
       height = 4,
       width = 13)

summary(lm(y ~ x))

# cat(paste0("R of estimates: ", round(cor(c(F1$mean), c(cor_after)), 3), "\n"))
cat(paste0("R of estimates: ", round(cor(x, y), 3), "\n"))

# ------------------------------------------------------------------------------
#   Is R^2 different for Johannes' "seasonal" taxa? (Ans: No)
# ------------------------------------------------------------------------------

if(FALSE) {
  seasonal_families <- list(wet = c("Helicobacteraceae",
                                    "Coriobacteriaceae",
                                    "Burkholderiaceae",
                                    "Eggerthellaceae",
                                    "Atopobiaceae",
                                    "Erysipelotrichaceae",
                                    "Lachnospiraceae"),
                            dry = c("Baceroidales RF16 group",
                                    "vadinBE97",
                                    "Spirochaetaceae",
                                    "Campylobacteraceae",
                                    "Christensenellaceae",
                                    "Syntrophomonadaceae"))

  # What's the R^2 for Johannes' seasonal families?
  pairs <- combn(1:D, m = 2)
  include_pairs <- which(pairs[1,] != D+1 & pairs[2,] != D+1)
  pair1 <- pairs[1,include_pairs]
  pair2 <- pairs[2,include_pairs]
  # MAP estimates
  sigma <- F1$mean[1:D,1:D]
  vsigma <- sigma[lower.tri(sigma,diag = FALSE)]
  vec_sigmas <- data.frame(rho = vsigma,
                           type = "original",
                           tax_idx1 = pair1,
                           tax_idx2 = pair2)
  sigma <- cor_after[1:D,1:D]
  vsigma <- sigma[lower.tri(sigma,diag = FALSE)]
  vec_sigmas <- rbind(vec_sigmas,
                      data.frame(rho = vsigma,
                                 type = "aperiodic",
                                 tax_idx1 = pair1,
                                 tax_idx2 = pair2))
  vec_sigmas <- vec_sigmas %>%
    pivot_wider(id_cols = c(tax_idx1, tax_idx2), names_from = type, values_from = rho)

  vec_sigmas$fam1 <- tax[vec_sigmas$tax_idx1,6]
  vec_sigmas$fam2 <- tax[vec_sigmas$tax_idx2,6]

  vec_sigmas$fam1[!(vec_sigmas$fam1 %in% seasonal_families$wet |
                      vec_sigmas$fam1 %in% seasonal_families$dry)] <- NA
  vec_sigmas$fam2[!(vec_sigmas$fam2 %in% seasonal_families$wet |
                      vec_sigmas$fam2 %in% seasonal_families$dry)] <- NA

  vec_sigmas$n_seasonal <- as.numeric(!sapply(vec_sigmas$fam1, is.na)) + as.numeric(!sapply(vec_sigmas$fam2, is.na))
  vec_sigmas$n_seasonal_factor <- factor(vec_sigmas$n_seasonal)
  levels(vec_sigmas$n_seasonal_factor) <- c("no seasonal partners", "one seasonal partner", "two seasonal partners")
  vec_sigmas$lachno <- vec_sigmas$fam1 == "Lachnospiraceae" & vec_sigmas$fam2 == "Lachnospiraceae"

  # R^2 is very similar for pairs where both partners are one of Johannes'
  # "seasonal" taxa
  n <- 2
  ggplot() +
    geom_point(data = vec_sigmas %>% filter(n_seasonal != n),
               mapping = aes(x = original, y = aperiodic),
               size = 2,
               color = "#bbbbbb") +
    geom_point(data = vec_sigmas %>% filter(n_seasonal == n),
               mapping = aes(x = original, y = aperiodic),
               size = 2,
               shape = 21,
               fill = "red",
               alpha = 0.5) +
    theme_bw()

  # ...or where both partners are lachno-lachno pairs
  ggplot() +
    geom_point(data = vec_sigmas %>% filter(!lachno),
               mapping = aes(x = original, y = aperiodic),
               size = 2,
               color = "#bbbbbb") +
    geom_point(data = vec_sigmas %>% filter(lachno),
               mapping = aes(x = original, y = aperiodic),
               size = 2,
               shape = 21,
               fill = "red",
               alpha = 0.5) +
    theme_bw()
}

# ------------------------------------------------------------------------------
#   Variance explained
# ------------------------------------------------------------------------------

# Get the variance explained from the seasonal component. For simplicity, I'm
# just going to use the difference in variance explained between models with and
# without this component.
ve_ar <- numeric(nrow(clr_counts))
ve_aper <- numeric(nrow(clr_counts))
for(i in 1:nrow(clr_counts)) {
  ve_ar[i] <- 1 - var(fits_ar[[i]]$residuals)/var(clr_counts[i,])
  ve_aper[i] <- 1 - var(fits_aper[[i]]$residuals)/var(clr_counts[i,])
}

ggplot(data.frame(x = ve_aper - ve_ar),
       aes(x = x)) +
  geom_histogram(color = "white", bins = 15) +
  theme_bw() +
  labs(x = "\nchange in percent variance explained",
       y = "count\n")

cat(paste0("Median change in percent variance explained w/ seasonal trend: ", round(median(ve_aper - ve_ar)*100,1 ), "\n"))
round(quantile(ve_aper - ve_ar, probs = c(0, 0.025, 0.5, 0.975, 1)), 3)
