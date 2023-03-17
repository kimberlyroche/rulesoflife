library(rulesoflife)
library(ggplot2)
library(magrittr)
library(dplyr)
library(caret)
library(ggridges)
library(tidyverse)
library(kableExtra)
library(cowplot)
library(conflicted)
library(ggfortify)
library(cowplot)

source("thresholds.R")

conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("mutate", "dplyr")
conflict_prefer("arrange", "dplyr")

data <- load_data(tax_level = "ASV")
filter_obj <- filter_joint_zeros(data$counts)

age_pairs <- list(c(0, 6),
                  c(6, 13),
                  c(13, 30))

# Juvenile vs. prime age with minimum number of samples (n=35)
group1 <- c("ACA", "CAI", "DAS", "EAG", "ELV", "LIW",
            "MBE", "ONY", "OPH", "OXY", "PEB", "VAI", "WYN")

# Prime age vs. older adult with minimum number of samples (n=35)
group2 <- c("COB", "ECH", "GAB", "HOL", "LAD", "NAP",
            "NOB", "TAL", "VEI", "VET", "VIN", "WAD", "WHE")

# Convert these to the labels used in the paper
# map <- read.table("output/host_labels.tsv", header = T)
# map %>% filter(sname %in% group1) %>% pull(host_label) %>% sort() %>% paste0(collapse = ", ")
# map %>% filter(sname %in% group2) %>% pull(host_label) %>% sort() %>% paste0(collapse = ", ")

# We'll split into age ranges **juvenile** (0-6 years), **prime age adult** (6-13
# years), and **older adult** (13+ years).

# Fit models for a handful of baboons for these ranges. We'll compare the CLR
# ASV-ASV correlations that result from the age range-restricted models.

# This only needs to be run once
if(FALSE) {
  # Fit prime age models
  for(sname in c(group1, group2)) {
    odir <- "asv_days90_diet25_scale1_age6-13"
    fit <- fit_GP(sname = sname,
                  counts = data$counts,
                  metadata = data$metadata,
                  output_dir = odir,
                  MAP = TRUE,
                  days_to_min_autocorrelation = 90,
                  diet_weight = 0.25,
                  age_min = 6,
                  age_max = 13,
                  scramble_sample = F,
                  scramble_spacing = F,
                  scramble_order = F,
                  use_adam = F)
  }

  # Fit juvenile models
  for(sname in group1) {
    odir <- "asv_days90_diet25_scale1_age0-6"
    fit <- fit_GP(sname = sname,
                  counts = data$counts,
                  metadata = data$metadata,
                  output_dir = odir,
                  MAP = TRUE,
                  days_to_min_autocorrelation = 90,
                  diet_weight = 0.25,
                  age_min = 0,
                  age_max = 6,
                  scramble_sample = F,
                  scramble_spacing = F,
                  scramble_order = F,
                  use_adam = F)
  }

  # Fit older adult models
  for(sname in group2) {
    odir <- "asv_days90_diet25_scale1_age13-30"
    fit <- fit_GP(sname = sname,
                  counts = data$counts,
                  metadata = data$metadata,
                  output_dir = odir,
                  MAP = TRUE,
                  days_to_min_autocorrelation = 90,
                  diet_weight = 0.25,
                  age_min = 13,
                  age_max = 30,
                  scramble_sample = F,
                  scramble_spacing = F,
                  scramble_order = F,
                  use_adam = F)
  }
}

# ------------------------------------------------------------------------------
# Peek at one host's correlation of CLR ASV-ASV correlation estimates derived
# from juvenile samples (x-axis) and estimates derived from prime age adult
# samples (y-axis)
# ------------------------------------------------------------------------------

if(FALSE) {
  sname <- group1[1]

  est1 <- fido::to_clr(readRDS(paste0("output/model_fits/asv_days90_diet25_scale1_age0-6/MAP/", sname, ".rds")))
  sigma1 <- cov2cor(est1$Sigma[,,1])

  est2 <- fido::to_clr(readRDS(paste0("output/model_fits/asv_days90_diet25_scale1_age6-13/MAP/", sname, ".rds")))
  sigma2 <- cov2cor(est2$Sigma[,,1])

  plot_df <- data.frame(x = c(sigma1[lower.tri(sigma1, diag = F)]),
                        y = c(sigma2[lower.tri(sigma2, diag = F)]))

  ggplot(plot_df, aes(x = x, y = y)) +
    geom_point() +
    xlim(c(-1,1)) +
    ylim(c(-1,1)) +
    theme_bw() +
    labs(x = "correlation (ages 0-6)",
         y = "correlation (ages 6-13)",
         title = paste0("R=", signif(cor(plot_df$x, plot_df$y), 2), " (", sname, ")"))

  # Look at the overall CLR ASV-ASV correlation matrix vs. juvenile vs. prime
  # age adult matrices

  est <- fido::to_clr(readRDS("output/model_fits/asv_days90_diet25_scale1/MAP/ACA.rds"))
  overall <- cov2cor(est$Sigma[,,1])

  p_overall <- plot_kernel_or_cov_matrix(overall) + theme(legend.position = "none")
  p_juvenile <- plot_kernel_or_cov_matrix(sigma1) + theme(legend.position = "none")
  p_prime <- plot_kernel_or_cov_matrix(sigma2) + theme(legend.position = "none")
  cowplot::plot_grid(plotlist = list(p_overall, p_juvenile, p_prime), ncol = 3)
}

# ------------------------------------------------------------------------------
# Peek at one host's correlation of CLR ASV-ASV correlation estimates derived
# from prime age adult samples (x-axis) and estimates derived from older adult
# samples (y-axis)
# ------------------------------------------------------------------------------

if(FALSE) {
  sname <- group2[1]

  est1 <- fido::to_clr(readRDS(paste0("output/model_fits/asv_days90_diet25_scale1_age6-13/MAP/", sname, ".rds")))
  sigma1 <- cov2cor(est1$Sigma[,,1])

  est2 <- fido::to_clr(readRDS(paste0("output/model_fits/asv_days90_diet25_scale1_age13-30/MAP/", sname, ".rds")))
  sigma2 <- cov2cor(est2$Sigma[,,1])

  plot_df <- data.frame(x = c(sigma1[lower.tri(sigma1, diag = F)]),
                        y = c(sigma2[lower.tri(sigma2, diag = F)]))

  ggplot(plot_df, aes(x = x, y = y)) +
    geom_point() +
    xlim(c(-1,1)) +
    ylim(c(-1,1)) +
    theme_bw() +
    labs(x = "correlation (ages 6-13)",
         y = "correlation (ages 13-30)",
         title = paste0("R=", signif(cor(plot_df$x, plot_df$y), 2), " (", sname, ")"))
}

# ------------------------------------------------------------------------------
# Build a "rug" object
# ------------------------------------------------------------------------------

pairs <- combn(1:125, m = 2)
pairs_df <- data.frame(idx1 = pairs[1,], idx2 = pairs[2,]) %>%
  left_join(filter_obj %>%
              select(idx1, idx2, threshold))
rug <- matrix(NA, length(group1)*2 + length(group2)*2, sum(filter_obj$threshold))
host_groups <- data.frame(index = 1:nrow(rug), sname = NA, group = NA)
for(i in 1:length(group1)) {
  sname <- group1[i]
  # Juvenile
  est1 <- fido::to_clr(readRDS(paste0("output/model_fits/asv_days90_diet25_scale1_age0-6/MAP/", sname, ".rds")))
  sigma1 <- cov2cor(est1$Sigma[1:125,1:125,1])
  # Prime age
  est2 <- fido::to_clr(readRDS(paste0("output/model_fits/asv_days90_diet25_scale1_age6-13/MAP/", sname, ".rds")))
  sigma2 <- cov2cor(est2$Sigma[1:125,1:125,1])

  rug[i,] <- c(sigma1[lower.tri(sigma1, diag = F)])[pairs_df$threshold]
  rug[(i+length(group1)),] <- c(sigma2[lower.tri(sigma2, diag = F)])[pairs_df$threshold]
  host_groups[i,]$sname <- sname
  host_groups[i,]$group <- "juvenile"
  host_groups[(i+length(group1)),]$sname <- sname
  host_groups[(i+length(group1)),]$group <- "prime"
}
offset <- length(group1)*2
for(i in 1:length(group2)) {
  sname <- group2[i]
  # Prime age
  est1 <- fido::to_clr(readRDS(paste0("output/model_fits/asv_days90_diet25_scale1_age6-13/MAP/", sname, ".rds")))
  sigma1 <- cov2cor(est1$Sigma[1:125,1:125,1])
  # Older adult
  est2 <- fido::to_clr(readRDS(paste0("output/model_fits/asv_days90_diet25_scale1_age13-30/MAP/", sname, ".rds")))
  sigma2 <- cov2cor(est2$Sigma[1:125,1:125,1])

  rug[(offset+i),] <- c(sigma1[lower.tri(sigma1, diag = F)])[pairs_df$threshold]
  rug[(offset+i+length(group2)),] <- c(sigma2[lower.tri(sigma2, diag = F)])[pairs_df$threshold]
  host_groups[(offset+i),]$sname <- sname
  host_groups[(offset+i),]$group <- "prime"
  host_groups[(offset+i+length(group2)),]$sname <- sname
  host_groups[(offset+i+length(group2)),]$group <- "older"
}

# 1 PCA plot
d <- dist(rug)
d_mat <- as.matrix(d)
coords <- cmdscale(d, k = 2)

# 2 PCA plots
coords1_data <- rug[1:(length(group1)*2),]
colnames(coords1_data) <- pairs_df %>%
  filter(threshold) %>%
  mutate(label = paste0("ASV", idx1, " x ASV", idx2)) %>%
  pull(label)
coords1 <- prcomp(coords1_data)

coords2_data <- rug[(length(group1)*2 + 1):nrow(rug),]
colnames(coords2_data) <- pairs_df %>%
  filter(threshold) %>%
  mutate(label = paste0("ASV", idx1, " x ASV", idx2)) %>%
  pull(label)
coords2 <- prcomp(coords2_data)

host_groups$x <- c(coords1$x[,1], coords2$x[,1])
host_groups$y <- c(coords2$x[,2], coords2$x[,2])

centroid_df <- host_groups %>%
  group_by(group) %>%
  summarize(x_mean = mean(x), y_mean = mean(y))

edge_df <- rbind(host_groups %>%
                   filter(sname %in% group1) %>%
                   group_by(sname) %>%
                   select(-c(index, y)) %>%
                   pivot_wider(id_cols = sname, names_from = group, values_from = x) %>%
                   rename(x = juvenile, xend = prime),
                 host_groups %>%
                   filter(sname %in% group2) %>%
                   group_by(sname) %>%
                   select(-c(index, y)) %>%
                   pivot_wider(id_cols = sname, names_from = group, values_from = x) %>%
                   rename(x = prime, xend = older)) %>%
  left_join(rbind(host_groups %>%
                    filter(sname %in% group1) %>%
                    group_by(sname) %>%
                    select(-c(index, x)) %>%
                    pivot_wider(id_cols = sname, names_from = group, values_from = y) %>%
                    rename(y = juvenile, yend = prime),
                  host_groups %>%
                    filter(sname %in% group2) %>%
                    group_by(sname) %>%
                    select(-c(index, x)) %>%
                    pivot_wider(id_cols = sname, names_from = group, values_from = y) %>%
                    rename(y = prime, yend = older)),
            by = "sname")

p1 <- ggplot() +
  geom_segment(data = edge_df %>%
                 filter(sname %in% group1),
               mapping = aes(x = x, xend = xend, y = y, yend = yend),
               color = "#aaaaaa") +
  geom_point(data = host_groups %>%
               filter(sname %in% group1) %>%
               mutate(group = factor(group, levels = c("juvenile", "prime"))),
             mapping = aes(x = x, y = y, fill = group),
             size = 3, shape = 21) +
  geom_point(data = centroid_df  %>%
               filter(group %in% c("juvenile", "prime")) %>%
               mutate(group = factor(group, levels = c("juvenile", "prime"))),
             mapping = aes(x = x_mean, y = y_mean, fill = group),
             size = 3, shape = 23, stroke = 1.25) +
  theme_bw() +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = c("#FC7300", "#BFDB38")) +
  # guides(color = FALSE) +
  labs(x = "PCA 1", y = "PCA 2", fill = "Age bin", color = "Age bin")

p2 <- ggplot() +
  geom_segment(data = edge_df %>%
                 filter(sname %in% group2),
               mapping = aes(x = x, xend = xend, y = y, yend = yend),
               color = "#aaaaaa") +
  geom_point(data = host_groups %>%
               filter(sname %in% group2) %>%
               mutate(group = factor(group, levels = c("prime", "older"))),
             mapping = aes(x = x, y = y, fill = group),
             size = 3, shape = 21) +
  geom_point(data = centroid_df %>%
               filter(group %in% c("older", "prime")) %>%
               mutate(group = factor(group, levels = c("prime", "older"))),
             mapping = aes(x = x_mean, y = y_mean, fill = group),
             size = 3, shape = 23, stroke = 1.25) +
  theme_bw() +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = c("#BFDB38", "#1F8A70")) +
  # guides(color = FALSE) +
  labs(x = "PCA 1", y = "PCA 2", fill = "Age bin", color = "Age bin")

p <- plot_grid(p1, p2, ncol = 2, labels = c("A", "B"), label_y = 1.015)

ggsave(filename = file.path("output", "figures", "Figure_5_Supplement_3.png"),
       plot = p, bg = "white", device = png(), dpi = 200, width = 8, height = 4.5, units = "in")

# Look at the loadings for the first couple of PCs. Which ASVs most contribute?
# We can visualize the loadings by hand in a biplot but this is messy

representation <- represented_taxa(filter_obj)
renumber_map <- data.frame(idx = 1:length(representation),
                           new_idx = sapply(1:length(representation), function(x) renumber_taxon(representation, x)))
pairs_thresholded <- pairs_df %>%
  filter(threshold) %>%
  left_join(renumber_map, by = c("idx1" = "idx")) %>%
  rename(idx1_renumbered = new_idx) %>%
  left_join(renumber_map, by = c("idx2" = "idx")) %>%
  rename(idx2_renumbered = new_idx)

# Find top three pairs associated with each loading
build_loadings_df <- function(coords_obj) {
  loadings <- cbind(pairs_thresholded,
                    data.frame(pair = 1:nrow(coords_obj$rotation),
                               pc1 = coords_obj$rotation[,1],
                               pc2 = coords_obj$rotation[,2]))

  loadings$taxonomy_1 <- sapply(loadings$idx1, function(x) {
    tax_levels <- colnames(data$taxonomy)[2:ncol(data$taxonomy)]
    paste(paste(tax_levels,
                data$taxonomy[x,2:ncol(data$taxonomy)]), collapse = " / ")
  })
  loadings$taxonomy_2 <- sapply(loadings$idx2, function(x) {
    tax_levels <- colnames(data$taxonomy)[2:ncol(data$taxonomy)]
    paste(paste(tax_levels,
                data$taxonomy[x,2:ncol(data$taxonomy)]), collapse = " / ")
  })

  loadings %<>%
    select(-c(idx1, idx2, threshold, pair)) %>%
    pivot_longer(cols = c(pc1, pc2),
                 names_to = "PC",
                 values_to = "loading") %>%
    mutate(PC = case_when(
      PC == "pc1" ~ 1,
      T ~ 2
    )) %>%
    select(idx1_renumbered, idx2_renumbered, PC, loading, taxonomy_1, taxonomy_2) %>%
    rename(idx1 = idx1_renumbered,
           idx2 = idx2_renumbered)
  loadings$abs_loading <- abs(loadings$loading)
  rownames(loadings) <- NULL
  loadings
}

loadings_juvenile <- cbind(groups = "juvenile vs. prime age adult",
                           build_loadings_df(coords1))
loadings_older <- cbind(groups = "older adult vs. prime age adult",
                        build_loadings_df(coords2))

# Combined table to output
table_out <- loadings_juvenile %>%
  filter(PC == 1) %>%
  arrange(desc(abs_loading)) %>%
  filter(row_number() <= 10) %>%
  rbind(loadings_juvenile %>%
          filter(PC == 2) %>%
          arrange(desc(abs_loading)) %>%
          filter(row_number() <= 10)) %>%
  rbind(loadings_older %>%
          filter(PC == 1) %>%
          arrange(desc(abs_loading)) %>%
          filter(row_number() <= 10)) %>%
  rbind(loadings_older %>%
          filter(PC == 2) %>%
          arrange(desc(abs_loading)) %>%
          filter(row_number() <= 10)) %>%
  select(-abs_loading) %>%
  rename(`Groups` = groups,
         `Identity of ASV 1` = idx1,
         `Identity of ASV 2` = idx2,
         `Principle component` = PC,
         `Loading` = loading)

write.table(table_out,
            file = file.path("output", "Table_S9.tsv"),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

# Biplot
# p0 <- autoplot(coords1, data = coords1_data, loadings = TRUE, loadings.label = 1) # color = 'species'
# p0$layers[[2]]$data <- p0$layers[[2]]$data[labels[top2],]
# p0$layers[[3]]$data <- p0$layers[[3]]$data[labels[top2],]
# p0

# lm_df <- NULL
# for(i in 1:nrow(host_groups)) {
#   if(host_groups$group[i] == "juvenile") {
#     n <- sum(data$metadata$age <= 6 & data$metadata$sname == host_groups$sname[i])
#   } else if(host_groups$group[i] == "prime") {
#     n <- sum(data$metadata$age >= 6 & data$metadata$age <= 13 & data$metadata$sname == host_groups$sname[i])
#   } else {
#     n <- sum(data$metadata$age >= 13 & data$metadata$sname == host_groups$sname[i])
#   }
#   lm_df <- rbind(lm_df,
#                  data.frame(obs_var = var(rug[i,]),
#                             group = host_groups$group[i],
#                             n = n))
# }
# lm_df$group <- factor(lm_df$group,
#                       levels = c("prime", "juvenile", "older"))
#
# summary(lm(obs_var ~ group + n, data = lm_df))

# ANOVA: Does being in the same age bin explain some of the variation we see in
# correlation estimates (for filtered, high-confidence pairs)?

hosts <- combn(1:52, m = 2)
dist_df <- data.frame(idx1 = hosts[1,], idx2 = hosts[2,])
dist_df$d <- d_mat[lower.tri(d_mat, diag = F)]

dist_df %<>%
  left_join(host_groups %>%
              select(index, sname, group),
            by = c("idx1" = "index")) %>%
  rename(sname1 = sname, group1 = group) %>%
  left_join(host_groups %>%
              select(index, sname, group),
            by = c("idx2" = "index")) %>%
  rename(sname2 = sname, group2 = group)

# ANOVA: Juvenile vs. prime

fit <- aov(d ~ same_group, data = dist_df %>%
             mutate(same_group = group1 == group2) %>%
             filter(group1 %in% c("juvenile", "prime"),
                    group2 %in% c("juvenile", "prime")))
s <- summary(fit)
cat(paste0("Juvenile vs. prime age p = ", signif(s[[1]]$`Pr(>F)`[1], 3), "\n"))
ve_age <- s[[1]]$`Sum Sq`[1]
ve_res <- s[[1]]$`Sum Sq`[2]
cat(paste0("\tVariance explained = ", round(ve_age/(ve_age+ve_res)*100, 1), "%\n"))

# ANOVA: Prime vs. older

fit <- aov(d ~ same_group, data = dist_df %>%
             mutate(same_group = group1 == group2) %>%
             filter(group1 %in% c("older", "prime"),
                    group2 %in% c("older", "prime")))
s <- summary(fit)
cat(paste0("Older vs. prime age p = ", signif(s[[1]]$`Pr(>F)`[1], 3), "\n"))
ve_age <- s[[1]]$`Sum Sq`[1]
ve_res <- s[[1]]$`Sum Sq`[2]
cat(paste0("\tVariance explained = ", round(ve_age/(ve_age+ve_res)*100, 1), "%\n"))

# Plot the "rug". Note that rows are juvenile (rows 1-13), prime age (rows 14-39),
# and older adult (rows 40-52).

col_order <- order(apply(rug, 2, median))

# Juvenile vs. prime age adult
use_rows <- host_groups %>%
  filter(group == "juvenile") %>%
  left_join(host_groups %>%
              filter(group == "prime"),
            by = "sname") %>%
  select(index.x, index.y)

p1_plot_rug <- rug[use_rows$index.x, col_order]
plot_rug <- cbind(1:nrow(p1_plot_rug), p1_plot_rug)
colnames(plot_rug) <- c("host", paste0(1:(ncol(plot_rug)-1)))
plot_rug <- pivot_longer(as.data.frame(plot_rug), !host, names_to = "pair", values_to = "correlation")
plot_rug$pair <- as.numeric(plot_rug$pair)

p1_distro <- plot_rug$correlation

p1 <- ggplot(plot_rug, aes(x = pair, y = host)) +
  geom_raster(aes(fill = correlation)) +
  scale_fill_gradient2(low = "navy", mid = "white", high = "red",
                       midpoint = 0) +
  labs(y = "") +
  theme_bw() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.text.x = element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none") +
  labs(fill = "correlation\nmetric")

p2_plot_rug <- rug[use_rows$index.y, col_order]
plot_rug <- cbind(1:nrow(p2_plot_rug), p2_plot_rug)
colnames(plot_rug) <- c("host", paste0(1:(ncol(plot_rug)-1)))
plot_rug <- pivot_longer(as.data.frame(plot_rug), !host, names_to = "pair", values_to = "correlation")
plot_rug$pair <- as.numeric(plot_rug$pair)

p2_distro <- plot_rug$correlation

p2 <- ggplot(plot_rug, aes(x = pair, y = host)) +
  geom_raster(aes(fill = correlation)) +
  scale_fill_gradient2(low = "navy", mid = "white", high = "red",
                       midpoint = 0) +
  labs(y = "") +
  theme_bw() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.text.x = element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none") +
  labs(fill = "correlation\nmetric")

# Prime age adult vs. older adult
use_rows <- host_groups %>%
  filter(group == "older") %>%
  left_join(host_groups %>%
              filter(group == "prime"),
            by = "sname") %>%
  select(index.x, index.y)

p3_plot_rug <- rug[use_rows$index.y, col_order]
plot_rug <- cbind(1:nrow(p3_plot_rug), p3_plot_rug)
colnames(plot_rug) <- c("host", paste0(1:(ncol(plot_rug)-1)))
plot_rug <- pivot_longer(as.data.frame(plot_rug), !host, names_to = "pair", values_to = "correlation")
plot_rug$pair <- as.numeric(plot_rug$pair)

p3_distro <- plot_rug$correlation

p3 <- ggplot(plot_rug, aes(x = pair, y = host)) +
  geom_raster(aes(fill = correlation)) +
  scale_fill_gradient2(low = "navy", mid = "white", high = "red",
                       midpoint = 0) +
  labs(y = "") +
  theme_bw() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.text.x = element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none") +
  labs(fill = "correlation\nmetric")

p4_plot_rug <- rug[use_rows$index.x, col_order]
plot_rug <- cbind(1:nrow(p4_plot_rug), p4_plot_rug)
colnames(plot_rug) <- c("host", paste0(1:(ncol(plot_rug)-1)))
plot_rug <- pivot_longer(as.data.frame(plot_rug), !host, names_to = "pair", values_to = "correlation")
plot_rug$pair <- as.numeric(plot_rug$pair)

p4_distro <- plot_rug$correlation

p4 <- ggplot(plot_rug, aes(x = pair, y = host)) +
  geom_raster(aes(fill = correlation)) +
  scale_fill_gradient2(low = "navy", mid = "white", high = "red",
                       midpoint = 0) +
  labs(y = "") +
  theme_bw() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.text.x = element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "bottom") +
  labs(fill = "correlation\nmetric")

legend <- get_legend(p4)
p4 <- p4 +
  theme(legend.position = "none")

pcol1 <- plot_grid(p1, p2, p3, p4, ncol = 1,
                   labels = c("A", "C", "E", "G"),
                   label_y = 1.04, label_x = -0.01,
                   rel_heights = c(1, 1, 1, 1.15))

pcol1_w_legend <- plot_grid(pcol1, legend, ncol = 1, rel_heights = c(1, 0.2))

plot_density <- function(values, label_x = F) {
  p <- ggplot(data.frame(x = values), aes(x = x, y = 1)) +
    # geom_density() +
    stat_density_ridges(quantile_lines = TRUE) +
    theme_bw() +
    xlim(c(-1, 1)) +
    labs(x = "CLR ASV-ASV correlation") +
    theme(axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          panel.grid = element_blank())
  if(!label_x) {
    p <- p +
      theme(axis.title.x = element_blank())
  }
  p
}

pcol2 <- plot_grid(plot_density(p1_distro),
                   plot_density(p2_distro),
                   plot_density(p3_distro),
                   plot_density(p4_distro, label_x = T),
                   NULL,
                   ncol = 1,
                   labels = c("B", "D", "F", "H", ""),
                   label_y = 1.04, label_x = -0.09,
                   rel_heights = c(1, 1, 1, 1.15, 0.85))

p8panel <- plot_grid(pcol1_w_legend, NULL, pcol2, ncol = 3, rel_widths = c(1, 0.05, 0.5))

ggsave(filename = file.path("output", "figures", "Figure_5_Supplement_1.png"),
       plot = p8panel, bg = "white", device = png(), dpi = 100, width = 8, height = 6, units = "in")

# ------------------------------------------------------------------------------
#
#   Titration-like analysis
#
# ------------------------------------------------------------------------------

minimizer_frequency <- function(rug, bin_size = 0.05, retain_dist = F) {
  N <- nrow(rug)
  global_mean <- colMeans(rug)
  mix <- seq(from = 0, to = 1, by = bin_size)
  mins <- NULL

  for(h in 1:N) {
    host_obs <- rug[h,]
    host_residual <- host_obs - global_mean

    for(j in 1:length(mix)) {
      combined_dynamics <- (1-mix[j])*global_mean + mix[j]*host_residual
      mins <- rbind(mins,
                    data.frame(host = h,
                               p = mix[j],
                               dist = as.numeric(dist(matrix(c(host_obs,
                                                               combined_dynamics),
                                                             byrow = T,
                                                             2,
                                                             length(host_obs))))))
    }
  }

  if(retain_dist) {
    mins
  } else {
    mins %>%
      group_by(host) %>%
      arrange(dist) %>%
      slice_head() %>%
      ungroup() %>%
      select(host, p)
  }
}

# A host component of 45% seems optimal for juveniles
j1 <- minimizer_frequency(p1_plot_rug)

# A host component of 40% seems optimal for prime age set 1
m1 <- minimizer_frequency(p2_plot_rug)

# A host component of 40% seems optimal for prime age set 2
m2 <- minimizer_frequency(p3_plot_rug)

# A host component of 35% seems optimal for older adults
o1 <- minimizer_frequency(p4_plot_rug)

minimizers <- rbind(cbind(bin = "juvenile",
                          j1 %>% select(p)),
                    cbind(bin = "prime age adult (juvenile-matched)",
                          m1 %>% select(p)),
                    cbind(bin = "prime age adult (older adult-matched)",
                          m2 %>% select(p)),
                    cbind(bin = "older adult",
                          o1 %>% select(p)))
minimizers$bin <- factor(minimizers$bin,
                         levels = c("juvenile",
                                    "prime age adult (juvenile-matched)",
                                    "prime age adult (older adult-matched)",
                                    "older adult"))
p <- minimizers %>%
  ggplot(aes(x = 1-p, fill = bin)) +
  geom_histogram(binwidth = 0.05, color = "white") +
  facet_wrap(. ~ bin, ncol = 1) +
  scale_fill_manual(values = c("#987284", "#9dbf9e", "#d0d6b5", "#ee7674")) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "population mean weight",
       y = "host count")

ggsave(filename = file.path("output", "figures", "Figure_5_Supplement_2.png"),
       plot = p,
       bg = "white",
       device = png(),
       dpi = 200,
       width = 4,
       height = 6,
       units = "in")

# Is the difference in proportion with at least 60% global contribution between
# the juvenile group and each of the other groups significant?
# j1_v_m1 <- prop.test(x = c(6, 10), n = c(13, 13), alternative = 'less')$p.value
# j1_v_m2 <- prop.test(x = c(6, 11), n = c(13, 13), alternative = 'less')$p.value
# j1_v_o1 <- prop.test(x = c(6, 9), n = c(13, 13), alternative = 'less')$p.value
# p.adjust(c(j1_v_m1, j1_v_m2, j1_v_o1), method = "BH")

# Are these differences in scatter from the mean significant? (F-test)
# No, not even before multiple test correction
# var.test(j1$p, m1$p) # juvenile vs. paired adults
# var.test(j1$p, o1$p) # juvenile vs. difference set in older age

# Present as a difference in host-level contribution instead
# I think this is too hard to interpret
if(FALSE) {
  minimizers <- rbind(data.frame(bin = "juvenile vs. prime age difference",
                                 diff = j1$p - m1$p),
                      data.frame(bin = "older vs. prime age difference",
                                 diff = o1$p - m2$p))
  minimizers$bin <- factor(minimizers$bin,
                           levels = c("juvenile vs. prime age difference",
                                      "older vs. prime age difference"))
  p <- minimizers %>%
    ggplot(aes(x = diff, fill = bin)) +
    geom_histogram(binwidth = 0.05, color = "white") +
    facet_wrap(. ~ bin, ncol = 1) +
    scale_x_continuous(breaks = seq(-2, 2, by = 0.05)) +
    scale_fill_manual(values = c("#987284", "#9dbf9e", "#d0d6b5", "#ee7674")) +
    theme_bw() +
    theme(legend.position = "none") +
    labs(x = "population mean weight",
         y = "host count")

  ggsave(filename = file.path("output", "figures", "Figure_5_Supplement_2.png"),
         plot = p,
         bg = "white",
         device = png(),
         dpi = 200,
         width = 4,
         height = 6,
         units = "in")
}

# Alternate plot

if(FALSE) {
  plot_minima <- function(df) {
    plot_df <- mins %>%
      left_join(data.frame(host = rep(1:13, 4),
                           sname = host_groups$sname),
                by = "host") %>%
      left_join(read.table("output/host_labels.tsv", sep = "\t", header = T),
                by = "sname")
    plot_df$group <- factor(plot_df$group,
                            levels = c("juvenile",
                                       "prime age adult (juvenile-matched)",
                                       "prime age adult (older-matched)",
                                       "older adult"))
    ggplot(plot_df, aes(x = p, y = dist, color = group, fill = host_label)) +
      geom_line(size = 1, alpha = 1) +
      scale_color_manual(values = c("#987284", "#9dbf9e", "#d0d6b5", "#ee7674")) +
      theme_bw() +
      scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +
      labs(x = "proportion host contribution",
           y = "Frobenius distance")
  }

  mins <- rbind(cbind(group = "juvenile",
                      minimizer_frequency(p1_plot_rug, retain_dist = T)),
                cbind(group = "prime age adult (juvenile-matched)",
                      minimizer_frequency(p2_plot_rug, retain_dist = T)),
                cbind(group = "prime age adult (older-matched)",
                      minimizer_frequency(p3_plot_rug, retain_dist = T)),
                cbind(group = "older adult",
                      minimizer_frequency(p4_plot_rug, retain_dist = T)))

  plot_minima(mins)
}
