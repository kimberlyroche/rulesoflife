source("path_fix.R")

library(rulesoflife)
library(fido)
library(phyloseq)
library(tidyverse) # stringr, dplyr, ggplot2
library(MASS)
library(reshape2)
library(lubridate)
library(ggridges)
library(RColorBrewer)
library(cowplot)
library(frechet)

# ------------------------------------------------------------------------------
#
#   Figure 4 - Comparing results to those from other data sets (Johnson et al.
#              and DIABIMMUNE); the "titration" experiment to estimate population-
#              versus host-level signal is performed at the end of this script
#
#   Supplemental Figure S11 - "rug" heat maps and estimated proportions of host-
#                             vs. population-level signal for the human data sets
#
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
#   Functions
# ------------------------------------------------------------------------------

# Simple version of filtering
# This assumes a matrix of taxa (rows) x samples (columns)
# Identify taxa with > 1/10% average relative abundance as those to retain
filter_taxa <- function(counts) {
  proportions <- as.matrix(counts)
  proportions <- apply(proportions, 2, function(x) x / sum(x))
  mean_rel_abundance <- rowMeans(proportions)
  mean_rel_abundance >= 0.001
}

# This version of filtering matches the Amboseli analyses
# Retain taxa observed at at least a single count in 20% of each subject's
# samples
filter_taxa2 <- function(host_columns, counts, min_relab = 1/10000) {
  relab <- apply(counts, 2, function(x) x/sum(x))
  retain_tax <- rep(TRUE, nrow(counts))
  for(h in 1:length(host_columns)) {
    sub_relab <- relab[,host_columns[[h]]]
    retain_tax <- retain_tax & apply(sub_relab, 1, function(x) {
      sum(x >= min_relab)/length(x) >= 0.2
    })
  }
  agglom_tax <- colSums(as.matrix(counts[!retain_tax,]))
  counts <- rbind(counts[retain_tax,],
                  agglom_tax)
  return(list(counts = counts, filter_vec = retain_tax))
}

# `counts` is taxa x samples
# `host_columns` is a list (length = num. hosts) of host column indices in the
#   `counts` object
# `host_dates` is a list (length = num. hosts) of host sample dates associated
#   with the columns in `counts`
fit_model <- function(counts, host_columns, host_dates, dataset_name) {
  D <- nrow(counts)
  matrix_Sigmas <- array(NA, dim = c(D, D, length(host_columns)))

  pair_indices <- NULL
  for(subject in 1:length(host_columns)) {
    cat(paste0("Evaluating ", dataset_name, " (", subject, " / ", length(host_columns), ")\n"))
    subject_samples <- host_columns[[subject]]
    subject_dates <- host_dates[[subject]] + 1
    if(length(subject_samples) > 0) {
      subject_counts <- counts[,subject_samples] # omit taxonomy

      Y <- as.matrix(subject_counts)
      X <- matrix(subject_dates, 1, length(subject_dates))

      alr_ys <- driver::alr((t(Y) + 0.5))
      alr_means <- colMeans(alr_ys)
      Theta <- function(X) matrix(alr_means, D-1, ncol(X))

      taxa_covariance <- get_Xi(D, total_variance = 1)

      # rho <- calc_se_decay(min_correlation = 0.1, days_to_min_autocorrelation = 7)
      rho <- calc_se_decay(min_correlation = 0.1, days_to_min_autocorrelation = 90)
      Gamma <- function(X) {
        SE(X, sigma = 1, rho = rho, jitter = 1e-08)
      }

      fit <- fido::basset(Y, X, taxa_covariance$upsilon, Theta, Gamma, taxa_covariance$Xi,
                          n_samples = 0, ret_mean = TRUE)

      fit.clr <- to_clr(fit)
      matrix_Sigmas[,,subject] <- cov2cor(fit.clr$Sigma[,,1])
    }
  }
  matrix_Sigmas
}

vectorize_results <- function(matrix_Sigmas) {
  hosts <- dim(matrix_Sigmas)[3]
  D <- dim(matrix_Sigmas)[1]
  n_interactions <- (D^2 - D)/2
  pair_indices <- NULL
  vectorized_Sigmas <- matrix(NA, hosts, n_interactions)
  for(k in 1:hosts) {
    temp <- matrix_Sigmas[,,k]
    if(is.null(pair_indices)) {
      pair_indices <- combn(D:1, m = 2)[2:1,((D^2 - D)/2):1]
    }
    temp <- c(temp[upper.tri(temp, diag = F)])
    vectorized_Sigmas[k,] <- temp
  }
  list(Sigmas = vectorized_Sigmas, pairs = pair_indices)
}

plot_heatmap2 <- function(rug, x_label, title) {
  canonical_col_order <- order(colMeans(rug))
  rug <- rug[,canonical_col_order]

  rug <- cbind(1:nrow(rug), rug)
  colnames(rug) <- c("host", paste0(1:(ncol(rug)-1)))
  rug <- pivot_longer(as.data.frame(rug), !host, names_to = "pair", values_to = "correlation")
  rug$pair <- as.numeric(rug$pair)

  p <- ggplot(rug, aes(x = pair, y = host)) +
    geom_raster(aes(fill = correlation)) +
    scale_fill_gradient2(low = "navy", mid = "white", high = "red", midpoint = 0,
                         guide = guide_colorbar(frame.colour = "black",
                                                ticks.colour = "black")) +
    scale_x_continuous(expand = c(0, 0)) +
    labs(fill = "Correlation",
         x = x_label,
         y = "",
         title = title) +
    scale_y_continuous(expand = c(0, 0)) +
    theme(axis.title = element_text(size = 14, face = "plain"),
          plot.title = element_text(size = 20, hjust = 0.5),
          legend.title = element_text(size = 14),
          plot.margin = margin(t = 20, r = 10, b = 10, l = 0),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank())
  p
}

get_density_obj <- function(scores) {
  scores <- scores + rnorm(nrow(scores)*ncol(scores), 0, 0.01)
  scores.kde <- kde2d(scores[,1], scores[,2], n = 400)
  dx <- diff(scores.kde$x[1:2]) # from emdbook::HPDregionplot()
  dy <- diff(scores.kde$y[1:2])
  sz <- sort(scores.kde$z)
  c1 <- cumsum(sz) * dx * dy
  dimnames(scores.kde$z) <- list(scores.kde$x, scores.kde$y)
  dc <- melt(scores.kde$z)
  dc$prob <- approx(sz, 1 - c1, dc$value)$y
  dc
}

# ------------------------------------------------------------------------------
#   Global data.frames to store: per-study universality scores and statistics
# ------------------------------------------------------------------------------

study_stats <- NULL
all_scores <- NULL

# ------------------------------------------------------------------------------
#   Johnson et al.
# ------------------------------------------------------------------------------

# Read count table
# Note that first column is taxonomy in the form k__XXX;p__YYY;c__ZZZ;...
counts <- read.delim(file.path("input", "johnson2019", "taxonomy_counts_s.txt"),
                     sep = "\t" )

# Read sample ID to subject ID mapping
mapping <- read.delim(file.path("input", "johnson2019", "SampleID_map.txt"),
                      sep = "\t")

# Basic filtering
# filter_vec <- filter_taxa(counts[,2:ncol(counts)])
# tax_johnson <- counts$taxonomy[filter_vec]
# counts <- counts[filter_vec,]
tax_johnson <- counts$taxonomy
counts <- counts[,2:ncol(counts)]

subjects <- unique(mapping$UserName)
host_columns <- list()
host_dates <- list()
for(i in 1:length(subjects)) {
  subject <- subjects[i]
  subject_sample_IDs <- mapping[mapping$UserName == subject,]$SampleID
  subject_days <- mapping[mapping$UserName == subject,]$StudyDayNo
  subject_columns <- which(colnames(counts) %in% subject_sample_IDs)
  subject_days <- which(subject_sample_IDs %in% colnames(counts)[subject_columns])
  if(length(subject_days) > 0) {
    subject_days <- subject_days - min(subject_days)
    idx <- length(host_columns) + 1
    host_columns[[idx]] <- subject_columns
    host_dates[[idx]] <- subject_days
  }
}

# Filter counts
counts <- filter_taxa2(host_columns, counts)$counts

# Stats
study_stats <- rbind(study_stats,
                     data.frame(name = "Johnson et al.",
                                n_subjects = length(subjects),
                                n_samples = ncol(counts),
                                n_taxa = nrow(counts),
                                tax_level = "species"))

saved_fn <- file.path("input", "johnson2019", "fitted_results.rds")
if(file.exists(saved_fn)) {
  johnson <- readRDS(saved_fn)
  johnson_vec <- vectorize_results(johnson)
} else {
  johnson <- fit_model(counts, host_columns, host_dates, dataset_name = "Johnson et al.")
  johnson_vec <- vectorize_results(johnson)
  saveRDS(johnson, saved_fn)
}
pairs_johnson <- johnson_vec$pairs
johnson_vec <- johnson_vec$Sigmas

# Plot rug
rug <- johnson_vec[,order(colMeans(johnson_vec))]
rug <- cbind(1:nrow(rug), rug)
colnames(rug) <- c("host", paste0(1:(ncol(rug)-1)))
rug <- pivot_longer(as.data.frame(rug), !host, names_to = "pair", values_to = "correlation")
rug$pair <- as.numeric(rug$pair)

s1 <-  ggplot(rug, aes(x = pair, y = host)) +
  geom_raster(aes(fill = correlation)) +
  scale_fill_gradient2(low = "navy", mid = "white", high = "red",
                       midpoint = 0) +
  labs(y = "host") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  labs(fill = "Correlation") +
  theme(axis.text = element_text(size = 7),
        axis.title = element_text(size = 12, face = "plain"),
        legend.title = element_text(size = 12))

scores <- apply(johnson_vec, 2, function(x) calc_universality_score(x, return_pieces = TRUE))
consensus_sign <- apply(johnson_vec, 2, calc_consensus_sign)

all_scores <- rbind(all_scores,
                    data.frame(X1 = scores[1,],
                               X2 = scores[2,],
                               idx1 = pairs_johnson[1,],
                               idx2 = pairs_johnson[2,],
                               sign = consensus_sign,
                               dataset = "Johnson et al.",
                               n_subjects = length(subjects)))

# ------------------------------------------------------------------------------
#   DIABIMMUNE
# ------------------------------------------------------------------------------

load(file.path("input", "DIABIMMUNE", "diabimmune_karelia_16s_data.rdata"))
load(file.path("input", "DIABIMMUNE", "DIABIMMUNE_Karelia_metadata.RData"))

counts <- t(data_16s)
rm(data_16s)
tax_diabimmune <- rownames(counts)
rownames(counts) <- NULL # omit taxonomy

# These aren't counts but relative abundances (summing to 5).
# Let's scale them into CPM.
counts <- counts*1e06/5
# filter_vec <- filter_taxa(counts)
# tax_diabimmune <- tax_diabimmune[filter_vec]
# counts <- counts[filter_vec,]

# Use kids with at least 15 samples
use_subjects <- names(which(table(metadata$subjectID) >= 15))

host_columns <- list()
host_dates <- list()
for(i in 1:ncol(counts)) {
  sample_ID <- colnames(counts)[i]
  metadata_idx <- which(metadata$SampleID == sample_ID)
  subject_ID <- metadata$subjectID[metadata_idx]
  if(subject_ID %in% use_subjects) {
    age_at_collection <- metadata$age_at_collection[metadata_idx]
    if(subject_ID %in% names(host_columns)) {
      host_columns[[subject_ID]] <- c(host_columns[[subject_ID]], i)
      host_dates[[subject_ID]] <- c(host_dates[[subject_ID]], age_at_collection)
    } else {
      host_columns[[subject_ID]] <- c(i)
      host_dates[[subject_ID]] <- age_at_collection
    }
  }
}

for(i in 1:length(host_dates)) {
  baseline_day <- min(host_dates[[i]])
  sample_days <- sapply(host_dates[[i]], function(x) x - baseline_day) + 1
  host_dates[[i]] <- sample_days
}

# Filter counts
counts <- filter_taxa2(host_columns, counts)$counts

# Stats
study_stats <- rbind(study_stats,
                     data.frame(name = "DIABIMMUNE",
                                n_subjects = length(use_subjects),
                                n_samples = ncol(counts),
                                n_taxa = nrow(counts),
                                tax_level = "unknown"))

saved_fn <- file.path("input", "DIABIMMUNE", "fitted_results.rds")
if(file.exists(saved_fn)) {
  diabimmune <- readRDS(saved_fn)
  diabimmune_vec <- vectorize_results(diabimmune)
} else {
  diabimmune <- fit_model(counts, host_columns, host_dates, "DIABIMMUNE")
  diabimmune_vec <- vectorize_results(diabimmune)
  saveRDS(diabimmune, saved_fn)
}
pairs_diabimmune <- diabimmune_vec$pairs
diabimmune_vec <- diabimmune_vec$Sigmas

# Plot rug
rug <- diabimmune_vec[,order(colMeans(diabimmune_vec))]
rug <- cbind(1:nrow(rug), rug)
colnames(rug) <- c("host", paste0(1:(ncol(rug)-1)))
rug <- pivot_longer(as.data.frame(rug), !host, names_to = "pair", values_to = "correlation")
rug$pair <- as.numeric(rug$pair)

s2 <- ggplot(rug, aes(x = pair, y = host)) +
  geom_raster(aes(fill = correlation)) +
  scale_fill_gradient2(low = "navy", mid = "white", high = "red",
                       midpoint = 0) +
  labs(y = "host") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  labs(fill = "Correlation") +
  theme(axis.text = element_text(size = 7),
        axis.title = element_text(size = 12, face = "plain"),
        legend.title = element_text(size = 12))

scores <- apply(diabimmune_vec, 2, function(x) calc_universality_score(x, return_pieces = TRUE))
consensus_sign <- apply(diabimmune_vec, 2, calc_consensus_sign)

all_scores <- rbind(all_scores,
                    data.frame(X1 = scores[1,],
                               X2 = scores[2,],
                               idx1 = pairs_diabimmune[1,],
                               idx2 = pairs_diabimmune[2,],
                               sign = consensus_sign,
                               dataset = "DIABIMMUNE",
                               n_subjects = length(use_subjects)))

# ------------------------------------------------------------------------------
#   Amboseli baboons
# ------------------------------------------------------------------------------

amboseli <- summarize_Sigmas("asv_days90_diet25_scale1")
scores <- apply(amboseli$rug, 2, function(x) calc_universality_score(x, return_pieces = TRUE))
consensus_sign <- apply(amboseli$rug, 2, calc_consensus_sign)
tax_abrp <- load_data(tax_level = "ASV")$taxonomy
all_scores <- rbind(all_scores,
                    data.frame(X1 = scores[1,],
                               X2 = scores[2,],
                               idx1 = amboseli$tax_idx1,
                               idx2 = amboseli$tax_idx2,
                               sign = consensus_sign,
                               dataset = "Amboseli",
                               n_subjects = 56))

# ------------------------------------------------------------------------------
#   Figure 4 panels
# ------------------------------------------------------------------------------

dataset_palette <- c(brewer.pal(9, "Spectral")[c(2,8,9)])
names(dataset_palette) <- sort(unique(all_scores$dataset))[c(1,2,3)]

plot_scores <- all_scores %>%
  mutate(prop_subjects = n_subjects / max(n_subjects))

p1 <- ggplot(plot_scores %>% filter(!(dataset %in% c("McMahon et al. (epilimnion)",
                                                    "McMahon et al. (hypolimnion)",
                                                    "Grossart et al."))),
             aes(x = X2, y = reorder(dataset, n_subjects), alpha = prop_subjects, fill = dataset)) +
  geom_density_ridges(stat = "binline", bins = 20, scale = 0.95, draw_baseline = TRUE) +
  theme_bw() +
  scale_alpha_continuous(range = c(0.25, 1.0)) +
  scale_fill_manual(values = dataset_palette) +
  scale_y_discrete(expand = expansion(add = c(0.35, 0.75))) +
  coord_cartesian(clip = "off") +
  guides(fill = "none",
         alpha = "none") +
  labs(x = "universality scores",
       y = "")

# Density plot from:
# https://stackoverflow.com/questions/23437000/how-to-plot-a-contour-line-showing-where-95-of-values-fall-within-in-r-and-in

# Specify desired contour levels
prob <- c(0.95, 0.5)

# I'm just using the studies with 5+ subjects
dc_combined <- rbind(cbind(get_density_obj(all_scores %>% filter(dataset == "DIABIMMUNE") %>% dplyr::select(X1, X2)), type = "DIABIMMUNE"),
                     cbind(get_density_obj(all_scores %>% filter(dataset == "Amboseli") %>% dplyr::select(X1, X2)), type = "Amboseli"),
                     cbind(get_density_obj(all_scores %>% filter(dataset == "Johnson et al.") %>% dplyr::select(X1, X2)), type = "Johnson et al."))

dc_combined$type <- factor(dc_combined$type, levels = c("DIABIMMUNE", "Johnson et al.", "Amboseli"))

p2 <- ggplot(data = dc_combined,
             aes(x = Var1, y = Var2, z = prob, fill = type, alpha = ..level..)) +
  geom_contour_filled(breaks = c(prob, 0)) +
  scale_alpha_discrete(range = c(0.4, 0.7)) +
  # scale_fill_manual(values = dataset_palette[c(2,7,8)]) +
  scale_fill_manual(values = dataset_palette) +
  theme_bw() +
  xlim(c(0.5, 1)) +
  ylim(c(0, 1)) +
  labs(x = "proportion shared sign",
       y = "median association strength",
       fill = "Data set") +
  guides(alpha = "none",
         color = "none")

p <- plot_grid(p1, NULL, p2,
               ncol = 3,
               rel_widths = c(1, 0.1, 1.2),
               scale = 0.95,
               labels = c("A", "B"),
               label_size = 20)

ggsave(file.path("output", "figures", "F4.png"),
       p,
       dpi = 100,
       units = "in",
       height = 4.5,
       width = 12)

# ------------------------------------------------------------------------------
#   Supplemental Figure 11 panels
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
#
#   "Titration"-type analyses
#
# ------------------------------------------------------------------------------

output_dir_full <- check_dir(c("output", "model_fits", "asv_days90_diet25_scale1", "MAP"))
file_list <- list.files(path = output_dir_full, pattern = "*.rds")
if(length(file_list) == 0) {
  output_dir_full <- check_dir(c("output", "model_fits", output_dir, "full_posterior"))
  file_list <- list.files(path = output_dir_full, pattern = "*.rds")
}
# Get taxa number and posterior sample number
fit <- readRDS(file.path(output_dir_full, file_list[1]))
D <- fit$D
amboseli <- array(NA, dim=c(D, D, length(file_list)))
for(f in 1:length(file_list)) {
  file <- file_list[f]
  fit <- readRDS(file.path(output_dir_full, file))
  # Convert to CLR
  if(fit$coord_system != "clr") {
    fit <- to_clr(fit)
  }
  amboseli[,,f] <- cov2cor(fit$Sigma[,,1])
}

datasets <- list("Amboseli" = amboseli,
                 "DIABIMMUNE" = diabimmune,
                 "Johnson et al." = johnson)
mins <- NULL
for(i in 1:length(datasets)) {
  dataset <- datasets[[i]]
  D <- dim(dataset)[1]
  N <- dim(dataset)[3]
  global_mean <- CovFMean(dataset)$Mout[[1]]
  addend <- diag(D)*1e-06
  mix <- seq(from = 0, to = 1, by = 0.05)
  for(h in 1:N) {
    host_obs <- dataset[,,h]
    host_residual <- host_obs - global_mean
    diag(host_residual) <- 1

    for(j in 1:length(mix)) {
      combined_dynamics <- (1-mix[j])*global_mean + mix[j]*host_residual
      mins <- rbind(mins,
                    data.frame(host = h,
                               p = mix[j],
                               dist = dist4cov(host_obs, combined_dynamics)$dist,
                               dataset = names(datasets)[i]))
    }
  }
}

minimizing_proportions <- mins %>%
  group_by(dataset, host) %>%
  arrange(dist) %>%
  slice(1) %>%
  ungroup()

minimizing_proportions$dataset <- factor(minimizing_proportions$dataset, levels = names(datasets)[3:1])

s3 <- ggplot(minimizing_proportions, aes(x = p, y = dataset, fill = dataset)) +
  geom_density_ridges(stat = "binline", binwidth = 0.05, scale = 0.95) +
  theme_bw() +
  scale_alpha_continuous(range = c(0.25, 1.0)) +
  scale_fill_manual(values = dataset_palette) +
  scale_y_discrete(expand = expansion(add = c(0.15, 1.05))) +
  coord_cartesian(clip = "off") +
  guides(fill = "none",
         alpha = "none") +
  labs(x = "host-level proportion",
       y = "")

p <- plot_grid(s1, s2, NULL, s3, ncol = 4,
               rel_widths = c(1, 1, -0.05, 0.8),
               labels = c("A", "B", "", "C"),
               label_size = 18,
               label_y = 1.01,
               label_x = -0.01,
               scale = 0.95)

ggsave(file.path("output", "figures", "S11.png"),
       p,
       dpi = 100,
       units = "in",
       height = 3,
       width = 12)

cat(paste0("Median host-level contribution (Amboseli): ",
           round(median(minimizing_proportions %>% filter(dataset == "Amboseli") %>% pull(p)), 2), "\n"))
cat(paste0("Median host-level contribution (DIABIMMUNE): ",
           round(median(minimizing_proportions %>% filter(dataset == "DIABIMMUNE") %>% pull(p)), 2), "\n"))
cat(paste0("Median host-level contribution (Johnson et al.): ",
           round(median(minimizing_proportions %>% filter(dataset == "Johnson et al.") %>% pull(p)), 2), "\n"))

# ------------------------------------------------------------------------------
#
#   ADDITIONAL ANALYSES
#
#   QUESTION 1: Are strong relationships (large median unsigned association)
#               the most universal? I.e. what is the rank-correlation of the
#               median association and the universality score? Does large median
#               association imply large universality score in this data set?
# ------------------------------------------------------------------------------

for(this_dataset in c("Johnson et al.", "DIABIMMUNE", "Amboseli")) {
  x <- all_scores[all_scores$dataset == this_dataset,]$X2
  y <- all_scores[all_scores$dataset == this_dataset,]$X1

  cat(paste0(toupper(this_dataset),
             " Spearman cor. betw. % agreement and median assoc.: ",
             round(cor(x, y, method = "spearman"), 3),
             "\n"))
}

# ------------------------------------------------------------------------------
#   QUESTION 2: Are the top 1% to 2.5% most universal pairs also skewed towards
#               positive associations? This doesn't make much sense to ask of
#               Johnson et al. because there is so little evidence of
#               "universality"...
# ------------------------------------------------------------------------------

for(this_dataset in c("DIABIMMUNE", "Amboseli")) {
  subset_scores <- all_scores %>%
    filter(dataset == this_dataset) %>%
    mutate(score = X1*X2) %>%
    arrange(desc(score))
  n_pairs <- nrow(subset_scores)
  top_1 <- round(n_pairs * 0.01)
  top_2p5 <- round(n_pairs * 0.025)
  pn_tally_1 <- subset_scores %>%
    slice(1:top_1) %>%
    group_by(sign) %>%
    tally()
  ppos_1 <- pn_tally_1 %>%
    mutate(n_all = sum(n)) %>%
    filter(sign == 1) %>%
    mutate(prop = n/n_all) %>%
    pull(prop)
  pn_tally_2p5 <- subset_scores %>%
    slice(1:top_2p5) %>%
    group_by(sign) %>%
    tally()
  ppos_2p5 <- pn_tally_2p5 %>%
    mutate(n_all = sum(n)) %>%
    filter(sign == 1) %>%
    mutate(prop = n/n_all) %>%
    pull(prop)
  cat(paste0(toupper(this_dataset),
             " proportion positive in top 1%: ",
             round(ppos_1, 2),
             " (of ", sum(pn_tally_1$n), " pairs)",
             "\n"))
  cat(paste0(toupper(this_dataset),
             " proportion positive in top 2.5%: ",
             round(ppos_2p5, 2),
             " (of ", sum(pn_tally_2p5$n), " pairs)",
             "\n"))
}

# ------------------------------------------------------------------------------
#   QUESTION 3: Are the most universal pairs from the same families?
# ------------------------------------------------------------------------------

# tax_map <- list("Johnson et al." = tax_johnson,
#                 "DIABIMMUNE" = tax_diabimmune)
# for(this_dataset in c("Johnson et al.", "DIABIMMUNE")) {
#   subset_scores <- all_scores %>%
#     filter(dataset == this_dataset) %>%
#     mutate(score = X1*X2) %>%
#     arrange(desc(score))
#
#   n_pairs <- nrow(subset_scores)
#   top_2p5 <- round(n_pairs * 0.025)
#   subset_scores$top <- c(rep(TRUE, top_2p5),
#                          rep(FALSE, n_pairs - top_2p5))
#   # Append families
#   if(this_dataset == "Johnson et al.") {
#     families <- unname(sapply(tax_johnson, function(x) {
#       str_replace(str_split(x, pattern = ";")[[1]][5], "f__", "")
#     }))
#   } else {
#     families <- unname(sapply(tax_diabimmune, function(x) {
#       pieces <- str_split(x, pattern = "\\|")[[1]]
#       if(length(pieces) < 5) {
#         ""
#       } else {
#         str_replace(pieces[5], "f__", "")
#       }
#     }))
#   }
#   families[families %in% c("", "NA")] <- ""
#   subset_scores$fam1 <- families[subset_scores$idx1]
#   subset_scores$fam2 <- families[subset_scores$idx2]
#
#   # Remove any with missing families
#   subset_scores <- subset_scores %>%
#     filter(fam1 != "" & fam2 != "")
#
#   subset_scores$taxpair <- paste0(subset_scores$fam1, " - ", subset_scores$fam2)
#
#   # ----------------------------------------------------------------------------
#   #   Enrichment of family pairs
#   # ----------------------------------------------------------------------------
#
#   frequencies <- table(subset_scores$taxpair)
#   frequencies_subset <- table(subset_scores %>% filter(top == TRUE) %>% pull(taxpair))
#
#   signif <- c()
#   for(fam in names(frequencies_subset)) {
#     fam_in_sample <- unname(unlist(frequencies_subset[fam]))
#     sample_size <- unname(unlist(sum(frequencies_subset)))
#     fam_in_bg <- unname(unlist(frequencies[fam]))
#     bg_size <- unname(unlist(sum(frequencies)))
#     ctab <- matrix(c(fam_in_sample,
#                      sample_size - fam_in_sample,
#                      fam_in_bg,
#                      bg_size - fam_in_bg),
#                    2, 2, byrow = TRUE)
#     prob <- fisher.test(ctab, alternative = "greater")$p.value
#     if(prob < 0.05) {
#       signif <- c(signif, fam)
#       cat(paste0("ASV family: ", fam, ", p-value: ", round(prob, 3), "\n"))
#     }
#   }
#
#   temp_palette <- generate_highcontrast_palette(length(signif))
#   names(temp_palette) <- signif
#
#   plot_enrichment(frequencies_subset1 = frequencies_subset,
#                   frequencies = frequencies,
#                   significant_families1 = signif,
#                   plot_height = 6,
#                   plot_width = 6,
#                   legend_topmargin = 100,
#                   use_pairs = TRUE,
#                   rel_widths = c(1, 0.35, 1, 0.25, 3),
#                   labels = c("overall", "top 2.5% pairs"),
#                   palette = temp_palette,
#                   save_name = paste0(this_dataset, "-enrichment-pair.png"))
#
#   # ----------------------------------------------------------------------------
#   #   Enrichment of families
#   # ----------------------------------------------------------------------------
#
#   frequencies <- table(c(subset_scores$fam1, subset_scores$fam2))
#   frequencies_subset <- table(c(subset_scores %>% filter(top == TRUE) %>% pull(fam1),
#                                 subset_scores %>% filter(top == TRUE) %>% pull(fam2)))
#
#   signif <- c()
#   for(fam in names(frequencies_subset)) {
#     fam_in_sample <- unname(unlist(frequencies_subset[fam]))
#     sample_size <- unname(unlist(sum(frequencies_subset)))
#     fam_in_bg <- unname(unlist(frequencies[fam]))
#     bg_size <- unname(unlist(sum(frequencies)))
#     ctab <- matrix(c(fam_in_sample,
#                      sample_size - fam_in_sample,
#                      fam_in_bg,
#                      bg_size - fam_in_bg),
#                    2, 2, byrow = TRUE)
#     prob <- fisher.test(ctab, alternative = "greater")$p.value
#     if(prob < 0.05) {
#       signif <- c(signif, fam)
#       cat(paste0("ASV family: ", fam, ", p-value: ", round(prob, 3), "\n"))
#     }
#   }
#
#   temp_palette <- generate_highcontrast_palette(length(signif))
#   names(temp_palette) <- signif
#
#   plot_enrichment(frequencies_subset1 = frequencies_subset,
#                   frequencies = frequencies,
#                   significant_families1 = signif,
#                   plot_height = 6,
#                   plot_width = 4.5,
#                   legend_topmargin = 100,
#                   use_pairs = FALSE,
#                   rel_widths = c(1, 0.35, 1, 0.1, 1.75),
#                   labels = c("overall", "top 2.5% pairs"),
#                   palette = temp_palette,
#                   save_name = paste0(this_dataset, "-enrichment.png"))
# }
