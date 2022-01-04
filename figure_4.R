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

# ------------------------------------------------------------------------------
#
#   Figure 4 - Comparing results to those from other data sets, both host-
#              associated and free-living bacterial communities
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
filter_taxa2 <- function(host_columns, counts) {
  retain_tax <- rep(TRUE, nrow(counts))
  for(h in 1:length(host_columns)) {
    sub_counts <- counts[,host_columns[[h]]]
    retain_tax <- retain_tax & apply(sub_counts, 1, function(x) {
      sum(x >= 1)/length(x) >= 0.2
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
# `depth` is the number of posterior samples to draw (depth = 1 is MAP
#   estimation)
fit_model <- function(counts, host_columns, host_dates, dataset_name, depth = 1) {
  D <- nrow(counts)
  n_interactions <- (D^2 - D)/2
  vectorized_Sigmas <- array(NA, dim = c(length(host_columns), n_interactions, depth))

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

      if(depth == 1) {
        n_samples <- 0
        ret_mean <- TRUE
      } else {
        n_samples <- 20
        ret_mean <- FALSE
      }

      fit <- fido::basset(Y, X, taxa_covariance$upsilon, Theta, Gamma, taxa_covariance$Xi,
                          n_samples = n_samples, ret_mean = ret_mean)

      fit.clr <- to_clr(fit)
      Sigmas <- fit.clr$Sigma
      for(k in 1:depth) {
        temp <- cov2cor(Sigmas[,,k])
        if(is.null(pair_indices)) {
          pair_indices <- combn(D:1, m = 2)[2:1,((D^2 - D)/2):1]
        }
        temp <- c(temp[upper.tri(temp, diag = F)])
        vectorized_Sigmas[subject,,k] <- temp
      }
    }
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
} else {
  johnson <- fit_model(counts, host_columns, host_dates, dataset_name = "Johnson et al.", depth = 1)
  saveRDS(johnson, saved_fn)
}
pairs_johnson <- johnson$pairs
johnson <- johnson$Sigmas[,,1]

# Plot rug
rug <- johnson[,order(colMeans(johnson))]
rug <- cbind(1:nrow(rug), rug)
colnames(rug) <- c("host", paste0(1:(ncol(rug)-1)))
rug <- pivot_longer(as.data.frame(rug), !host, names_to = "pair", values_to = "correlation")
rug$pair <- as.numeric(rug$pair)

p <- ggplot(rug, aes(x = pair, y = host)) +
  geom_raster(aes(fill = correlation)) +
  scale_fill_gradient2(low = "navy", mid = "white", high = "red",
                       midpoint = 0) +
  labs(y = "") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  labs(fill = "correlation\nmetric") +
  theme(axis.text = element_text(size = 7),
          axis.title = element_text(size = 12, face = "plain"),
          legend.title = element_text(size = 12))

ggsave(file.path("output", "figures", "Johnson_rug.png"),
       plot = p,
       dpi = 100,
       units = "in",
       height = 4,
       width = 6)

scores <- apply(johnson, 2, function(x) calc_universality_score(x, return_pieces = TRUE))
consensus_sign <- apply(johnson, 2, calc_consensus_sign)

all_scores <- rbind(all_scores,
                    data.frame(X1 = scores[1,],
                               X2 = scores[2,],
                               idx1 = pairs_johnson[1,],
                               idx2 = pairs_johnson[2,],
                               sign = consensus_sign,
                               dataset = "Johnson et al.",
                               n_subjects = length(subjects)))

# ------------------------------------------------------------------------------
#   Dethlefsen & Relman
# ------------------------------------------------------------------------------

# Read count table
counts <- read.delim(file.path("input", "dethlefsen_relman2011", "sd01.txt"),
                     sep = "\t")

# Basic filtering
# counts <- counts[filter_taxa(counts[,2:163]),]
counts <- counts[,2:163]

# Read the sampling schedule metadata
metadata <- read.delim(file.path("input", "dethlefsen_relman2011", "sampling_schedule.txt"),
                       header = FALSE,
                       sep = "\t")
colnames(metadata) <- c("sample", "date")
metadata$date <- as.Date(mdy(metadata$date))

host_columns <- list(unname(which(sapply(metadata$sample, function(x) str_detect(x, "^D") ))),
                     unname(which(sapply(metadata$sample, function(x) str_detect(x, "^E") ))),
                     unname(which(sapply(metadata$sample, function(x) str_detect(x, "^F") ))))

host_dates <- list()
for(i in 1:3) {
  subject_dates <- metadata[host_columns[[i]],]$date
  host_dates[[i]] <- sapply(subject_dates, function(x) {
    difftime(x, subject_dates[1], units = "days")
  })
}

# Filter counts
counts <- filter_taxa2(host_columns, counts)$counts

# Stats
study_stats <- rbind(study_stats,
                     data.frame(name = "Dethlefsen & Relman",
                                n_subjects = 3,
                                n_samples = ncol(counts),
                                n_taxa = nrow(counts),
                                tax_level = "unknown"))

saved_fn <- file.path("input", "dethlefsen_relman2011", "fitted_results.rds")
if(file.exists(saved_fn)) {
  dethlefsen_relman <- readRDS(saved_fn)
} else {
  dethlefsen_relman <- fit_model(counts, host_columns, host_dates, "Dethlefsen & Relman", depth = 1)
  saveRDS(dethlefsen_relman, saved_fn)
}
dethlefsen_relman <- dethlefsen_relman$Sigmas[,,1]
scores <- apply(dethlefsen_relman, 2, function(x) calc_universality_score(x, return_pieces = TRUE))
consensus_sign <- apply(dethlefsen_relman, 2, calc_consensus_sign)

all_scores <- rbind(all_scores,
                    data.frame(X1 = scores[1,],
                               X2 = scores[2,],
                               sign = consensus_sign,
                               dataset = "Dethlefsen & Relman",
                               n_subjects = 3))

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
} else {
  diabimmune <- fit_model(counts, host_columns, host_dates, "DIABIMMUNE", depth = 1)
  saveRDS(diabimmune, saved_fn)
}
pairs_diabimmune <- diabimmune$pairs
diabimmune <- diabimmune$Sigmas[,,1]

# Plot rug
rug <- diabimmune[,order(colMeans(diabimmune))]
rug <- cbind(1:nrow(rug), rug)
colnames(rug) <- c("host", paste0(1:(ncol(rug)-1)))
rug <- pivot_longer(as.data.frame(rug), !host, names_to = "pair", values_to = "correlation")
rug$pair <- as.numeric(rug$pair)

p <- ggplot(rug, aes(x = pair, y = host)) +
  geom_raster(aes(fill = correlation)) +
  scale_fill_gradient2(low = "navy", mid = "white", high = "red",
                       midpoint = 0) +
  labs(y = "") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  labs(fill = "correlation\nmetric") +
  theme(axis.text = element_text(size = 7),
        axis.title = element_text(size = 12, face = "plain"),
        legend.title = element_text(size = 12))

ggsave(file.path("output", "figures", "DIABIMMUNE_rug.png"),
       plot = p,
       dpi = 100,
       units = "in",
       height = 4,
       width = 6)

scores <- apply(diabimmune, 2, function(x) calc_universality_score(x, return_pieces = TRUE))
consensus_sign <- apply(diabimmune, 2, calc_consensus_sign)

all_scores <- rbind(all_scores,
                    data.frame(X1 = scores[1,],
                               X2 = scores[2,],
                               idx1 = pairs_diabimmune[1,],
                               idx2 = pairs_diabimmune[2,],
                               sign = consensus_sign,
                               dataset = "DIABIMMUNE",
                               n_subjects = length(use_subjects)))

# ------------------------------------------------------------------------------
#   Grossart et al.
# ------------------------------------------------------------------------------

# Read count table
counts <- readRDS(file.path("input", "grossart_lakes", "data.rds"))
counts <- otu_table(counts)@.Data

# Filter data
# counts <- counts[filter_taxa(counts),]

metadata <- read.delim(file.path("input", "grossart_lakes", "945_20190102-074005.txt"),
                       header = TRUE,
                       sep = "\t")

# Look at the sample numbers per unique combo of depth and filter poresize for
# each lake. We'll hand pick similar conditions with large sample numbers!
conditions <- list(c("freshwater metagenome, Lake Breiter Luzin", "0-5", "0.2micron"),
                   c("freshwater metagenome, Lake Grosse Fuchskuhle", "0-2", "0.2micron"),
                   c("freshwater metagenome, Lake Melzer", "0-1", "0.2micron"),
                   c("freshwater metagenome, Lake Stechlin", "0-10", "0.2micron"),
                   c("freshwater metagenome, Lake Tiefwaren", "0-10", "0.2micron"))

# cycle through and collect the samples for all of these
host_columns <- list()
host_dates <- list()
for(i in 1:length(conditions)) {
  focal_md <- metadata[metadata$description == conditions[[i]][1] &
                         metadata$depth == conditions[[i]][2] &
                         metadata$filter_poresize == conditions[[i]][3],]
  sample_names <- focal_md$sample_name
  timestamps <- focal_md$collection_timestamp
  reorder <- order(timestamps)
  sample_names <- sample_names[reorder]
  timestamps <- timestamps[reorder]

  # it looks like some sample names aren't in the count table?
  matched_sample_names <- sample_names %in% colnames(counts)
  sample_names <- sample_names[matched_sample_names]
  timestamps <- timestamps[matched_sample_names]

  baseline_date <- timestamps[1]
  days_vector <- unname(sapply(timestamps, function(x) difftime(as.Date(x), as.Date(baseline_date), units = "days")) + 1)
  host_columns[[i]] <- which(colnames(counts) %in% sample_names)
  host_dates[[i]] <- days_vector
}

# Filter counts
counts <- filter_taxa2(host_columns, counts)$counts

# Stats
study_stats <- rbind(study_stats,
                     data.frame(name = "Grossart et al.",
                                n_subjects = length(conditions),
                                n_samples = ncol(counts),
                                n_taxa = nrow(counts),
                                tax_level = "unknown"))

saved_fn <- file.path("input", "grossart_lakes", "fitted_results.rds")
if(file.exists(saved_fn)) {
  grossart <- readRDS(saved_fn)
} else {
  grossart <- fit_model(counts, host_columns, host_dates, "Grossart", depth = 1)
  saveRDS(grossart, saved_fn)
}
grossart <- grossart$Sigmas[,,1]
scores <- apply(grossart, 2, function(x) calc_universality_score(x, return_pieces = TRUE))
consensus_sign <- apply(grossart, 2, calc_consensus_sign)

all_scores <- rbind(all_scores,
                    data.frame(X1 = scores[1,],
                               X2 = scores[2,],
                               sign = consensus_sign,
                               dataset = "Grossart et al.",
                               n_subjects = length(conditions)))

# ------------------------------------------------------------------------------
#   McMahon et al.
# ------------------------------------------------------------------------------

# Read count table
counts <- readRDS(file.path("input", "mcmahon_lakes", "data.rds"))
counts <- otu_table(counts)@.Data

# Filter data
# counts <- counts[filter_taxa(counts),]

metadata <- read.delim(file.path("input", "mcmahon_lakes", "1288_20180418-110149.txt"),
                       header = TRUE,
                       sep = "\t")

# lake short name to list index lookup
# lake_names <- c("MA", "TB", "CB", "NS", "WSB", "SSB", "HK", "NSB")
metadata$lake <- NA
metadata$layer <- NA
# append more useful columns to the metadata
for(i in 1:nrow(metadata)) {
  matches <- str_match_all(metadata$description[i], regex("^freshwater metagenome (\\D+)(\\d+)(\\D+)(\\d+)"))
  lake_name <- matches[[1]][1,2]
  if(lake_name != "NSb") {
    layer <- str_match_all(lake_name, regex("[H|E]$"))
    if(nrow(layer[[1]]) > 0) {
      if(layer[[1]][1,1] == "H") {
        # hypolimnion; deep water
        metadata$lake[i] <- substr(lake_name, 1, str_length(lake_name)-1)
        metadata$layer[i] <- "H"
      } else {
        # epilimnion; shallow water
        metadata$lake[i] <- substr(lake_name, 1, str_length(lake_name)-1)
        metadata$layer[i] <- "E"
      }
    } else {
      metadata$lake[i] <- lake_name
    }
  }
}

host_columns_E <- list()
host_columns_H <- list()
host_dates_E <- list()
host_dates_H <- list()
for(i in 1:ncol(counts)) {
  sample_id <- colnames(counts)[i]
  metadata_idx <- which(metadata$sample_name == sample_id)
  lake_name <- metadata$lake[metadata_idx]
  layer <- metadata$layer[metadata_idx]
  if(!is.na(layer) & !is.na(lake_name)) { # "NSb", the singleton
    sample_date <- as.Date(metadata$collection_timestamp[metadata_idx])
    if(layer == "E") {
      if(lake_name %in% names(host_columns_E)) {
        host_columns_E[[lake_name]] <- c(host_columns_E[[lake_name]], i)
        host_dates_E[[lake_name]] <- c(host_dates_E[[lake_name]], sample_date)
      } else {
        host_columns_E[[lake_name]] <- c(i)
        host_dates_E[[lake_name]] <- sample_date
      }
    } else {
      if(lake_name %in% names(host_columns_H)) {
        host_columns_H[[lake_name]] <- c(host_columns_H[[lake_name]], i)
        host_dates_H[[lake_name]] <- c(host_dates_H[[lake_name]], sample_date)
      } else {
        host_columns_H[[lake_name]] <- c(i)
        host_dates_H[[lake_name]] <- sample_date
      }
    }
  }
}

# convert dates to days
for(i in 1:length(host_dates_E)) {
  baseline_date <- min(host_dates_E[[i]])
  sample_days <- sapply(host_dates_E[[i]], function(x) difftime(x, baseline_date, units = "days")) + 1
  host_dates_E[[i]] <- sample_days
}
for(i in 1:length(host_dates_H)) {
  baseline_date <- min(host_dates_H[[i]])
  sample_days <- sapply(host_dates_H[[i]], function(x) difftime(x, baseline_date, units = "days")) + 1
  host_dates_H[[i]] <- sample_days
}

# Filter counts
counts_E <- filter_taxa2(host_columns_E, counts)$counts
counts_H <- filter_taxa2(host_columns_H, counts)$counts

# Stats
study_stats <- rbind(study_stats,
                     data.frame(name = "McMahon et al.",
                                n_subjects = length(unique(metadata$lake)),
                                n_samples = ncol(counts),
                                n_taxa = min(nrow(counts_E), nrow(counts_H)),
                                tax_level = "unknown"))

saved_fn <- file.path("input", "mcmahon_lakes", "fitted_results.rds")
if(file.exists(saved_fn)) {
  mcmahon <- readRDS(saved_fn)
  mc_e <- mcmahon[[1]]
  mc_h <- mcmahon[[2]]
  rm(mcmahon)
} else {
  mc_e <- fit_model(counts_E, host_columns_E, host_dates_E, "McMahon (epilimnion; shallow water)", depth = 1)
  mc_h <- fit_model(counts_H, host_columns_H, host_dates_H, "McMahon (hypolimnion; deep water)", depth = 1)
  saveRDS(list(mc_e, mc_h), saved_fn)
}
mc_e <- mc_e$Sigmas[,,1]
mc_h <- mc_h$Sigmas[,,1]
scores <- apply(mc_e, 2, function(x) calc_universality_score(x, return_pieces = TRUE))
consensus_sign <- apply(mc_e, 2, calc_consensus_sign)
all_scores <- rbind(all_scores,
                    data.frame(X1 = scores[1,],
                               X2 = scores[2,],
                               sign = consensus_sign,
                               dataset = "McMahon et al. (epilimnion)",
                               n_subjects = length(unique(metadata$lake))))
scores <- apply(mc_h, 2, function(x) calc_universality_score(x, return_pieces = TRUE))
consensus_sign <- apply(mc_h, 2, calc_consensus_sign)
all_scores <- rbind(all_scores,
                    data.frame(X1 = scores[1,],
                               X2 = scores[2,],
                               sign = consensus_sign,
                               dataset = "McMahon et al. (hypolimnion)",
                               n_subjects = length(unique(metadata$lake))))

# ------------------------------------------------------------------------------
#   David et al.
#
#   Note: This data set takes a long time to fit!
# ------------------------------------------------------------------------------

# Parse counts
counts <- read.delim(file.path("input", "david2014", "otu.table.ggref"),
                     header = TRUE,
                     sep = "\t")

metadata <- read.table(file.path("input", "david2014", "13059_2013_3286_MOESM18_ESM.csv"),
                       header = TRUE,
                       sep = ",")

sample_labels <- metadata$X
subj_labels <- metadata$AGE # subject A is 26, subject B is 36
day_labels <- metadata$COLLECTION_DAY
subj_labels[subj_labels == 26] <- "A"
subj_labels[subj_labels == 36] <- "B"

# Strip character from the count column names
sample_labels <- sapply(sample_labels, function(x) {
  str_replace(x, "\\.\\d+", "")
})
sample_labels <- unname(sample_labels)

# Remove saliva samples
keep_idx <- which(!sapply(sample_labels, function(x) {
  str_detect(x, "Saliva")
}))
subj_labels <- subj_labels[keep_idx]
sample_labels <- sample_labels[keep_idx]
day_labels <- day_labels[keep_idx]

# Unclude columns in counts that are in sample_labels
counts <- counts[,colnames(counts) %in% sample_labels]
keep_idx <- sample_labels %in% colnames(counts)
sample_labels <- sample_labels[keep_idx]
subj_labels <- subj_labels[keep_idx]
day_labels <- day_labels[keep_idx]

subj_A_idx <- c()
subj_A_days <- c()
subj_B_idx <- c()
subj_B_days <- c()
for(i in 1:length(colnames(counts))) {
  idx <- which(sample_labels == colnames(counts)[i])
  if(subj_labels[idx] == "A") {
    subj_A_idx <- c(subj_A_idx, i)
    subj_A_days <- c(subj_A_days, day_labels[idx])
  } else {
    subj_B_idx <- c(subj_B_idx, i)
    subj_B_days <- c(subj_B_days, day_labels[idx])
  }
}

subj_A_counts <- counts[,subj_A_idx]
subj_B_counts <- counts[,subj_B_idx]

reorder <- order(subj_A_days)
subj_A_days <- subj_A_days[reorder]
subj_A_counts <- subj_A_counts[,reorder]

reorder <- order(subj_B_days)
subj_B_days <- subj_B_days[reorder]
subj_B_counts <- subj_B_counts[,reorder]

counts <- cbind(subj_A_counts, subj_B_counts)

# Filter data
# counts <- counts[filter_taxa(counts),]

host_columns <- list(1:length(subj_A_days), (length(subj_A_days) + 1):ncol(counts)) # individuals 1, 2
host_dates <- list(subj_A_days, subj_B_days)

# Filter counts
counts <- filter_taxa2(host_columns, counts)$counts

# Stats
study_stats <- rbind(study_stats,
                     data.frame(name = "David et al.",
                                n_subjects = 2,
                                n_samples = ncol(counts),
                                n_taxa = nrow(counts),
                                tax_level = "ASV/OTU"))

saved_fn <- file.path("input", "david2014", "fitted_results.rds")
if(file.exists(saved_fn)) {
  david <- readRDS(saved_fn)
} else {
  david <- fit_model(counts, host_columns, host_dates, "David et al.", depth = 1)
  saveRDS(david, saved_fn)
}
david <- david$Sigmas[,,1]
scores <- apply(david, 2, function(x) calc_universality_score(x, return_pieces = TRUE))
consensus_sign <- apply(david, 2, calc_consensus_sign)
all_scores <- rbind(all_scores,
                    data.frame(X1 = scores[1,],
                               X2 = scores[2,],
                               sign = consensus_sign,
                               dataset = "David et al.",
                               n_subjects = 2))

# ------------------------------------------------------------------------------
#   Caporaso et al.
# ------------------------------------------------------------------------------

# Read count tables; OTUs all agree -- already checked
counts_1 <- read.table(file.path("input", "caporaso2011", "F4_feces_L6.txt"),
                       header = FALSE,
                       sep = "\t")
counts_1 <- counts_1[,2:ncol(counts_1)]

counts_2 <- read.table(file.path("input", "caporaso2011", "M3_feces_L6.txt"),
                       header = FALSE,
                       sep = "\t")
counts_2 <- counts_2[,2:ncol(counts_2)]

counts <- cbind(counts_1, counts_2)

# Filter low abundance taxa
# The first row consists of sample day identifiers
sample_days <- counts[1,]
# counts <- counts[c(TRUE, filter_taxa(counts[2:nrow(counts),])),]
counts <- counts[2:nrow(counts),]

host_columns <- list(1:ncol(counts_1), (ncol(counts_1)+1):ncol(counts)) # individuals 1, 2
host_dates <- list()
for(i in 1:2) {
  # host_dates[[i]] <- unname(unlist(counts[1,host_columns[[i]]]))
  host_dates[[i]] <- unname(unlist(sample_days[host_columns[[i]]]))
}

# Filter counts
counts <- filter_taxa2(host_columns, counts)$counts

# Stats
study_stats <- rbind(study_stats,
                     data.frame(name = "Caporaso et al.",
                                n_subjects = 2,
                                n_samples = ncol(counts),
                                n_taxa = nrow(counts),
                                tax_level = "unknown"))

saved_fn <- file.path("input", "caporaso2011", "fitted_results.rds")
if(file.exists(saved_fn)) {
  caporaso <- readRDS(saved_fn)
} else {
  caporaso <- fit_model(counts[2:nrow(counts),], host_columns, host_dates, "Caporaso et al.", depth = 1)
  saveRDS(caporaso, saved_fn)
}
caporaso <- caporaso$Sigmas[,,1]
scores <- apply(caporaso, 2, function(x) calc_universality_score(x, return_pieces = TRUE))
consensus_sign <- apply(caporaso, 2, calc_consensus_sign)
all_scores <- rbind(all_scores,
                    data.frame(X1 = scores[1,],
                               X2 = scores[2,],
                               sign = consensus_sign,
                               dataset = "Caporaso et al.",
                               n_subjects = 2))

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
#   QUESTION 1: Are strong relationships (large median unsigned association)
#               the most universal?
# ------------------------------------------------------------------------------

if(FALSE) {
  for(this_dataset in c("Johnson et al.", "DIABIMMUNE", "Amboseli")) {
    x <- all_scores[all_scores$dataset == this_dataset,]$X2
    y <- all_scores[all_scores$dataset == this_dataset,]$X1

    cat(paste0(toupper(this_dataset),
               " Spearman cor. betw. % agreement and median assoc.: ",
               round(cor(x, y, method = "spearman"), 3),
               "\n"))
  }
}

# ------------------------------------------------------------------------------
#   QUESTION 2: Are the top 1% to 2.5% most universal pairs also skewed towards
#               positive associations?
# ------------------------------------------------------------------------------

if(FALSE) {
  for(this_dataset in c("Johnson et al.", "DIABIMMUNE", "Amboseli")) {
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
}

# ------------------------------------------------------------------------------
#   QUESTION 3: Are the most universal pairs from the same families?
# ------------------------------------------------------------------------------

if(FALSE) {
  tax_map <- list("Johnson et al." = tax_johnson,
                  "DIABIMMUNE" = tax_diabimmune)
  for(this_dataset in c("Johnson et al.", "DIABIMMUNE")) {
    subset_scores <- all_scores %>%
      filter(dataset == this_dataset) %>%
      mutate(score = X1*X2) %>%
      arrange(desc(score))

    n_pairs <- nrow(subset_scores)
    top_2p5 <- round(n_pairs * 0.025)
    subset_scores$top <- c(rep(TRUE, top_2p5),
                           rep(FALSE, n_pairs - top_2p5))
    # Append families
    if(this_dataset == "Johnson et al.") {
      families <- unname(sapply(tax_johnson, function(x) {
        str_replace(str_split(x, pattern = ";")[[1]][5], "f__", "")
      }))
    } else {
      families <- unname(sapply(tax_diabimmune, function(x) {
        pieces <- str_split(x, pattern = "\\|")[[1]]
        if(length(pieces) < 5) {
          ""
        } else {
          str_replace(pieces[5], "f__", "")
        }
      }))
    }
    families[families %in% c("", "NA")] <- ""
    subset_scores$fam1 <- families[subset_scores$idx1]
    subset_scores$fam2 <- families[subset_scores$idx2]

    # Remove any with missing families
    subset_scores <- subset_scores %>%
      filter(fam1 != "" & fam2 != "")

    subset_scores$taxpair <- paste0(subset_scores$fam1, " - ", subset_scores$fam2)

    # ----------------------------------------------------------------------------
    #   Enrichment of family pairs
    # ----------------------------------------------------------------------------

    frequencies <- table(subset_scores$taxpair)
    frequencies_subset <- table(subset_scores %>% filter(top == TRUE) %>% pull(taxpair))

    signif <- c()
    for(fam in names(frequencies_subset)) {
      fam_in_sample <- unname(unlist(frequencies_subset[fam]))
      sample_size <- unname(unlist(sum(frequencies_subset)))
      fam_in_bg <- unname(unlist(frequencies[fam]))
      bg_size <- unname(unlist(sum(frequencies)))
      ctab <- matrix(c(fam_in_sample,
                       sample_size - fam_in_sample,
                       fam_in_bg,
                       bg_size - fam_in_bg),
                     2, 2, byrow = TRUE)
      prob <- fisher.test(ctab, alternative = "greater")$p.value
      if(prob < 0.05) {
        signif <- c(signif, fam)
        cat(paste0("ASV family: ", fam, ", p-value: ", round(prob, 3), "\n"))
      }
    }

    temp_palette <- generate_highcontrast_palette(length(signif))
    names(temp_palette) <- signif

    plot_enrichment(frequencies_subset1 = frequencies_subset,
                    frequencies = frequencies,
                    significant_families1 = signif,
                    plot_height = 6,
                    plot_width = 6,
                    legend_topmargin = 100,
                    use_pairs = TRUE,
                    rel_widths = c(1, 0.35, 1, 0.25, 3),
                    labels = c("overall", "top 2.5% pairs"),
                    palette = temp_palette,
                    save_name = paste0(this_dataset, "-enrichment-pair.png"))

    # ----------------------------------------------------------------------------
    #   Enrichment of families
    # ----------------------------------------------------------------------------

    frequencies <- table(c(subset_scores$fam1, subset_scores$fam2))
    frequencies_subset <- table(c(subset_scores %>% filter(top == TRUE) %>% pull(fam1),
                                  subset_scores %>% filter(top == TRUE) %>% pull(fam2)))

    signif <- c()
    for(fam in names(frequencies_subset)) {
      fam_in_sample <- unname(unlist(frequencies_subset[fam]))
      sample_size <- unname(unlist(sum(frequencies_subset)))
      fam_in_bg <- unname(unlist(frequencies[fam]))
      bg_size <- unname(unlist(sum(frequencies)))
      ctab <- matrix(c(fam_in_sample,
                       sample_size - fam_in_sample,
                       fam_in_bg,
                       bg_size - fam_in_bg),
                     2, 2, byrow = TRUE)
      prob <- fisher.test(ctab, alternative = "greater")$p.value
      if(prob < 0.05) {
        signif <- c(signif, fam)
        cat(paste0("ASV family: ", fam, ", p-value: ", round(prob, 3), "\n"))
      }
    }

    temp_palette <- generate_highcontrast_palette(length(signif))
    names(temp_palette) <- signif

    plot_enrichment(frequencies_subset1 = frequencies_subset,
                    frequencies = frequencies,
                    significant_families1 = signif,
                    plot_height = 6,
                    plot_width = 4.5,
                    legend_topmargin = 100,
                    use_pairs = FALSE,
                    rel_widths = c(1, 0.35, 1, 0.1, 1.75),
                    labels = c("overall", "top 2.5% pairs"),
                    palette = temp_palette,
                    save_name = paste0(this_dataset, "-enrichment.png"))
  }
}

# ------------------------------------------------------------------------------
#   VISUALIZATION VERSION 1: Plot densities together
# ------------------------------------------------------------------------------

# Define a color palette
dataset_palette <- c(brewer.pal(9, "Spectral")[c(1:4,6:9)], "#aaaaaa")
names(dataset_palette) <- sort(unique(all_scores$dataset))[c(2,1,3,4,6,8,5,7,9)]

# Density plot from:
# https://stackoverflow.com/questions/23437000/how-to-plot-a-contour-line-showing-where-95-of-values-fall-within-in-r-and-in

# Specify desired contour levels
prob <- c(0.95, 0.5)

# I'm just using the studies with 5+ subjects
dc_combined <- rbind(cbind(get_density_obj(all_scores %>% filter(dataset == "DIABIMMUNE") %>% dplyr::select(X1, X2)), type = "DIABIMMUNE"),
                     cbind(get_density_obj(all_scores %>% filter(dataset == "Amboseli") %>% dplyr::select(X1, X2)), type = "Amboseli"),
                     cbind(get_density_obj(all_scores %>% filter(dataset == "Johnson et al.") %>% dplyr::select(X1, X2)), type = "Johnson et al."))

dc_combined$type <- factor(dc_combined$type, levels = c("DIABIMMUNE", "Johnson et al.", "Amboseli"))

p <- ggplot(data = dc_combined,
            aes(x = Var1, y = Var2, z = prob, fill = type, alpha = ..level..)) +
  geom_contour_filled(breaks = c(prob, 0)) +
  scale_alpha_discrete(range = c(0.4, 0.7)) +
  scale_fill_manual(values = dataset_palette[c(2,7,8)]) +
  theme_bw() +
  xlim(c(0.5, 1)) +
  ylim(c(0, 1)) +
  labs(x = "proportion shared sign",
       y = "median association strength",
       fill = "Data set") +
  guides(alpha = "none",
         color = "none")

ggsave(file.path("output", "figures", "6-alt.png"),
       p,
       dpi = 100,
       units = "in",
       height = 5,
       width = 6.5)

# ------------------------------------------------------------------------------
#   VISUALIZATION VERSION 2: Stacked histograms
# ------------------------------------------------------------------------------

plot_scores <- all_scores %>%
  mutate(prop_subjects = n_subjects / max(n_subjects))

p1 <- ggplot(plot_scores %>% filter(!(dataset %in% c("McMahon et al. (epilimnion)",
                                                    "McMahon et al. (hypolimnion)",
                                                    "Grossart et al."))),
             aes(x = X2, y = reorder(dataset, n_subjects), alpha = prop_subjects, fill = dataset)) +
  geom_density_ridges(stat = "binline", bins = 20, scale = 0.95, draw_baseline = TRUE) +
  theme_bw() +
  # scale_alpha_continuous(range = c(0.1, 1.0)) +
  scale_fill_manual(values = dataset_palette) +
  scale_y_discrete(expand = expansion(add = c(0.35, 0.75))) +
  coord_cartesian(clip = "off") +
  guides(fill = "none",
         alpha = "none") +
  labs(x = "universality scores",
       y = "")

p2 <- ggplot(plot_scores %>% filter(dataset %in% c("McMahon et al. (epilimnion)",
                                                   "McMahon et al. (hypolimnion)",
                                                   "Grossart et al.")),
             aes(x = X2, y = reorder(dataset, n_subjects), alpha = prop_subjects, fill = dataset)) +
  geom_density_ridges(stat = "binline", bins = 20, scale = 0.5, draw_baseline = TRUE) +
  theme_bw() +
  scale_alpha_continuous(range = c(0.089, 0.161)) +
  scale_fill_manual(values = dataset_palette) +
  scale_y_discrete(expand = expansion(add = c(0.35, 0.75))) +
  coord_cartesian(clip = "off") +
  guides(fill = "none",
         alpha = "none") +
  labs(x = "universality scores",
       y = "")

p3 <- plot_grid(NULL, p1,
                ncol = 2,
                rel_widths = c(0.07, 1))

p <- plot_grid(p3, NULL, p2,
               ncol = 1,
               rel_heights = c(1, 0.1, 0.55))

ggsave(file.path("output", "figures", "6.png"),
       p,
       dpi = 100,
       units = "in",
       height = 7,
       width = 6)

