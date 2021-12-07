source("path_fix.R")

library(rulesoflife)
library(fido)
library(phyloseq)
library(tidyverse) # stringr, dplyr, ggplot2
library(MASS)
library(reshape2)
library(lubridate)
library(ggridges)

# ------------------------------------------------------------------------------
#   Functions
# ------------------------------------------------------------------------------

# This assumes a matrix of taxa (rows) x samples (columns)
# Identify taxa with > 1/10% average relative abundance
# TO DO: Remove rare in any individual; similar filtering to ABRP
filter_taxa <- function(counts) {
  proportions <- as.matrix(counts)
  proportions <- apply(proportions, 2, function(x) x / sum(x))
  mean_rel_abundance <- rowMeans(proportions)
  mean_rel_abundance >= 0.001
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

      rho <- calc_se_decay(min_correlation = 0.1, days_to_min_autocorrelation = 7)
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
        temp <- c(temp[upper.tri(temp, diag = F)])
        vectorized_Sigmas[subject,,k] <- temp
      }
    }
  }
  vectorized_Sigmas
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
counts <- counts[filter_taxa(counts[,2:ncol(counts)]),]

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
  johnson <- fit_model(counts, host_columns, host_dates, dataset_name = "Johnson et al.", depth = 1)[,,1]
  saveRDS(johnson, saved_fn)
}
scores <- apply(johnson, 2, function(x) calc_universality_score(x, return_pieces = TRUE))

all_scores <- rbind(all_scores,
                    data.frame(X1 = scores[1,],
                               X2 = scores[2,],
                               dataset = "Johnson et al.",
                               n_subjects = length(subjects)))

# ------------------------------------------------------------------------------
#   Dethlefsen & Relman
# ------------------------------------------------------------------------------

# Read count table
counts <- read.delim(file.path("input", "dethlefsen_relman2011", "sd01.txt"),
                     sep = "\t")

# Basic filtering
counts <- counts[filter_taxa(counts[,2:163]),]

# Read the sampling schedule metadata
metadata <- read.delim(file.path("input", "dethlefsen_relman2011", "sampling_schedule.txt"),
                       header = FALSE,
                       sep = "\t")
colnames(metadata) <- c("sample", "date")
metadata$date <- as.Date(mdy(metadata$date))

host_samples <- list(unname(which(sapply(metadata$sample, function(x) str_detect(x, "^D") ))),
                     unname(which(sapply(metadata$sample, function(x) str_detect(x, "^E") ))),
                     unname(which(sapply(metadata$sample, function(x) str_detect(x, "^F") ))))

host_dates <- list()
for(i in 1:3) {
  subject_dates <- metadata[host_samples[[i]],]$date
  host_dates[[i]] <- sapply(subject_dates, function(x) {
    difftime(x, subject_dates[1], units = "days")
  })
}

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
  dethlefsen_relman <- fit_model(counts, host_samples, host_dates, "Dethlefsen & Relman", depth = 1)[,,1]
  saveRDS(dethlefsen_relman, saved_fn)
}
scores <- apply(dethlefsen_relman, 2, function(x) calc_universality_score(x, return_pieces = TRUE))

all_scores <- rbind(all_scores,
                    data.frame(X1 = scores[1,],
                               X2 = scores[2,],
                               dataset = "Dethlefsen & Relman",
                               n_subjects = 3))

# ------------------------------------------------------------------------------
#   DIABIMMUNE
# ------------------------------------------------------------------------------

load(file.path("input", "DIABIMMUNE", "diabimmune_karelia_16s_data.rdata"))
load(file.path("input", "DIABIMMUNE", "DIABIMMUNE_Karelia_metadata.RData"))

counts <- t(data_16s)
rm(data_16s)
rownames(counts) <- NULL # omit taxonomy

# These aren't counts but relative abundances (summing to 5).
# Let's scale them into CPM.
counts <- counts*1e06/5
counts <- counts[filter_taxa(counts),]

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
  diabimmune <- fit_model(counts, host_columns, host_dates, "DIABIMMUNE", depth = 1)[,,1]
  saveRDS(diabimmune, saved_fn)
}
scores <- apply(diabimmune, 2, function(x) calc_universality_score(x, return_pieces = TRUE))

all_scores <- rbind(all_scores,
                    data.frame(X1 = scores[1,],
                               X2 = scores[2,],
                               dataset = "DIABIMMUNE",
                               n_subjects = length(use_subjects)))

# ------------------------------------------------------------------------------
#   Grossart et al.
# ------------------------------------------------------------------------------

# Read count table
counts <- readRDS(file.path("input", "grossart_lakes", "data.rds"))
counts <- otu_table(counts)@.Data

# Filter data
counts <- counts[filter_taxa(counts),]

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
  grossart <- fit_model(counts, host_columns, host_dates, "Grossart", depth = 1)[,,1]
  saveRDS(grossart, saved_fn)
}
scores <- apply(grossart, 2, function(x) calc_universality_score(x, return_pieces = TRUE))

all_scores <- rbind(all_scores,
                    data.frame(X1 = scores[1,],
                               X2 = scores[2,],
                               dataset = "Grossart et al.",
                               n_subjects = length(conditions)))

# ------------------------------------------------------------------------------
#   McMahon et al.
# ------------------------------------------------------------------------------

# Read count table
counts <- readRDS(file.path("input", "mcmahon_lakes", "data.rds"))
counts <- otu_table(counts)@.Data

# Filter data
counts <- counts[filter_taxa(counts),]

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

# Stats
study_stats <- rbind(study_stats,
                     data.frame(name = "McMahon et al.",
                                n_subjects = length(unique(metadata$lake)),
                                n_samples = ncol(counts),
                                n_taxa = nrow(counts),
                                tax_level = "unknown"))

saved_fn <- file.path("input", "mcmahon_lakes", "fitted_results.rds")
if(file.exists(saved_fn)) {
  mcmahon <- readRDS(saved_fn)
  mc_e <- mcmahon[[1]]
  mc_h <- mcmahon[[2]]
  rm(mcmahon)
} else {
  mc_e <- fit_model(counts, host_columns_E, host_dates_E, "McMahon (epilimnion; shallow water)", depth = 1)[,,1]
  mc_h <- fit_model(counts, host_columns_H, host_dates_H, "McMahon (hypolimnion; deep water)", depth = 1)[,,1]
  saveRDS(list(mc_e, mc_h), saved_fn)
}
scores <- apply(mc_e, 2, function(x) calc_universality_score(x, return_pieces = TRUE))
all_scores <- rbind(all_scores,
                    data.frame(X1 = scores[1,],
                               X2 = scores[2,],
                               dataset = "McMahon et al. (epilimnion)",
                               n_subjects = length(unique(metadata$lake))))
scores <- apply(mc_h, 2, function(x) calc_universality_score(x, return_pieces = TRUE))
all_scores <- rbind(all_scores,
                    data.frame(X1 = scores[1,],
                               X2 = scores[2,],
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
counts <- counts[filter_taxa(counts),]

host_columns <- list(1:length(subj_A_days), (length(subj_A_days) + 1):ncol(counts)) # individuals 1, 2
host_dates <- list(subj_A_days, subj_B_days)

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
  david <- fit_model(counts, host_columns, host_dates, "David et al.", depth = 1)[,,1]
  saveRDS(david, saved_fn)
}
scores <- apply(david, 2, function(x) calc_universality_score(x, return_pieces = TRUE))
all_scores <- rbind(all_scores,
                    data.frame(X1 = scores[1,],
                               X2 = scores[2,],
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
counts <- counts[c(TRUE, filter_taxa(counts[2:nrow(counts),])),]

host_columns <- list(1:ncol(counts_1), (ncol(counts_1)+1):ncol(counts)) # individuals 1, 2
host_dates <- list()
for(i in 1:2) {
  host_dates[[i]] <- unname(unlist(counts[1,host_columns[[i]]]))
}

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
  caporaso <- fit_model(counts[2:nrow(counts),], host_columns, host_dates, "Caporaso et al.", depth = 1)[,,1]
  saveRDS(caporaso, saved_fn)
}
scores <- apply(caporaso, 2, function(x) calc_universality_score(x, return_pieces = TRUE))
all_scores <- rbind(all_scores,
                    data.frame(X1 = scores[1,],
                               X2 = scores[2,],
                               dataset = "Caporaso et al.",
                               n_subjects = 2))

# ------------------------------------------------------------------------------
#   Amboseli baboons
# ------------------------------------------------------------------------------

amboseli <- summarize_Sigmas("asv_days90_diet25_scale1")
scores <- apply(amboseli$rug, 2, function(x) calc_universality_score(x, return_pieces = TRUE))
all_scores <- rbind(all_scores,
                    data.frame(X1 = scores[1,],
                               X2 = scores[2,],
                               dataset = "Amboseli",
                               n_subjects = 56))

# ------------------------------------------------------------------------------
#   VISUALIZATION VERSION 1: Plot densities together
# ------------------------------------------------------------------------------

# Density plot from:
# https://stackoverflow.com/questions/23437000/how-to-plot-a-contour-line-showing-where-95-of-values-fall-within-in-r-and-in

# Specify desired contour levels
prob <- c(0.95, 0.5)

# I'm just using the studies with 5+ subjects
dc_combined <- rbind(cbind(get_density_obj(all_scores %>% filter(dataset == "DIABIMMUNE") %>% select(X1, X2)), type = "A"),
                     cbind(get_density_obj(all_scores %>% filter(dataset == "McMahon et al. (hypolimnion)") %>% select(X1, X2)), type = "B"),
                     cbind(get_density_obj(all_scores %>% filter(dataset == "Amboseli") %>% select(X1, X2)), type = "C"),
                     cbind(get_density_obj(all_scores %>% filter(dataset == "Johnson et al.") %>% select(X1, X2)), type = "D"))

dc_combined$type <- factor(dc_combined$type)
levels(dc_combined$type) <- c("DIABIMMUNE human infants",
                              "McMahon et al. lakes",
                              "Amboseli baboons",
                              "Johnson et al. human adults")

p <- ggplot(data = dc_combined,
            aes(x = Var1, y = Var2, z = prob, color = type, fill = type, alpha = ..level..)) +
  geom_contour_filled(breaks = c(prob, 0)) +
  scale_alpha_discrete(range = c(0.4, 0.7)) +
  theme_bw() +
  xlim(c(0.45, 1.05)) +
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
       height = 6,
       width = 8)

# ------------------------------------------------------------------------------
#   VISUALIZATION VERSION 2: Stacked histograms
# ------------------------------------------------------------------------------

p <- ggplot(all_scores %>% filter(dataset != "McMahon et al. (epilimnion)"), aes(x = X2, y = reorder(dataset, n_subjects), alpha = n_subjects, fill = dataset)) +
  geom_density_ridges(stat = "binline", bins = 20, scale = 0.95, draw_baseline = TRUE) +
  # geom_density_ridges(scale = 1.5, draw_baseline = TRUE) +
  theme_bw() +
  # scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(expand = expansion(add = c(0.35, 0.75))) +
  coord_cartesian(clip = "off") +
  guides(fill = "none") +
  labs(x = "universality scores",
       y = "",
       alpha = "No. subjects")

ggsave(file.path("output", "figures", "6.png"),
       p,
       dpi = 100,
       units = "in",
       height = 6,
       width = 7)

