source("path_fix.R")

library(rulesoflife)
library(fido)
library(stringr)
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

      # fit <- fido::basset(Y, X, taxa_covariance$upsilon, Theta, Gamma, taxa_covariance$Xi,
      #                     n_samples = n_samples, ret_mean = ret_mean,
      #                     b2 = 0.98, step_size = 0.004, eps_f = 1e-11, eps_g = 1e-05,
      #                     max_iter = 10000L, optim_method = "adam")
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
#   Johnson et al.
# ------------------------------------------------------------------------------

saved_fn <- file.path("input", "johnson2019", "fitted_results.rds")
if(file.exists(saved_fn)) {
  johnson <- readRDS(saved_fn)
} else {
  # Read count table
  # Note that first column is taxonomy in the form k__XXX;p__YYY;c__ZZZ;...
  counts <- read.delim(file.path("input", "johnson2019", "taxonomy_counts_s.txt"),
                       sep = "\t" )

  # Read sample ID to subject ID mapping
  mapping <- read.delim(file.path("input", "johnson2019", "SampleID_map.txt"),
                        sep = "\t")

  # Basic filtering
  counts <- counts[filter_taxa(counts[,2:ncol(counts)]),]

  # cat(paste0("Johnson et al.\n\t", nrow(counts), " species x ", ncol(counts), " samples\n"))

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
  johnson <- fit_model(counts, host_columns, host_dates, dataset_name = "Johnson et al.", depth = 1)[,,1]
  saveRDS(johnson, saved_fn)
}
# plot_heatmap2(johnson, "species-species pairs", "Johnson et al.")
johnson_scores <- apply(johnson, 2, function(x) calc_universality_score(x, return_pieces = TRUE))
johnson_scores <- data.frame(X1 = johnson_scores[1,], X2 = johnson_scores[2,])

# ------------------------------------------------------------------------------
#   Dethlefsen & Relman
# ------------------------------------------------------------------------------

saved_fn <- file.path("input", "dethlefsen_relman2011", "fitted_results.rds")
if(file.exists(saved_fn)) {
  dethlefsen_relman <- readRDS(saved_fn)
} else {
  # Read count table
  counts <- read.delim(file.path("input", "dethlefsen_relman2011", "sd01.txt"),
                       sep = "\t")

  # Basic filtering
  counts <- counts[filter_taxa(counts[,2:163]),]

  # cat(paste0("Dethlefsen & Relman\n\t", nrow(counts), " ASVs x ", ncol(counts), " samples\n"))

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
  dethlefsen_relman <- fit_model(counts, host_samples, host_dates, "Dethlefsen & Relman", depth = 1)[,,1]
  saveRDS(dethlefsen_relman, saved_fn)
}
# plot_heatmap2(dethlefsen_relman, "species-species pairs", "Dethlefsen & Relman")
dr_scores <- apply(dethlefsen_relman, 2, function(x) calc_universality_score(x, return_pieces = TRUE))
dr_scores <- data.frame(X1 = dr_scores[1,], X2 = dr_scores[2,])

# ------------------------------------------------------------------------------
#   DIABIMMUNE
# ------------------------------------------------------------------------------

saved_fn <- file.path("input", "DIABIMMUNE", "fitted_results.rds")
if(file.exists(saved_fn)) {
  diabimmune <- readRDS(saved_fn)
} else {
  load(file.path("input", "DIABIMMUNE", "diabimmune_karelia_16s_data.rdata"))
  load(file.path("input", "DIABIMMUNE", "DIABIMMUNE_Karelia_metadata.RData"))

  counts <- t(data_16s)
  rm(data_16s)
  rownames(counts) <- NULL # omit taxonomy

  # These aren't counts but relative abundances (summing to 5).
  # Let's scale them into CPM.
  counts <- counts*1e06/5
  counts <- counts[filter_taxa(counts),]

  # cat(paste0("DIABIMMUNE\n\t", nrow(counts), " genera x ", ncol(counts), " samples\n"))

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
  diabimmune <- fit_model(counts, host_columns, host_dates, "DIABIMMUNE", depth = 1)[,,1]
  saveRDS(diabimmune, saved_fn)
}
# plot_heatmap2(diabimmune, "genus-genus pairs", "DIABIMMUNE")
diab_scores <- apply(diabimmune, 2, function(x) calc_universality_score(x, return_pieces = TRUE))
diab_scores <- data.frame(X1 = diab_scores[1,], X2 = diab_scores[2,])

# ------------------------------------------------------------------------------
#   Grossart et al.
# ------------------------------------------------------------------------------

saved_fn <- file.path("input", "grossart_lakes", "fitted_results.rds")
if(file.exists(saved_fn)) {
  grossart <- readRDS(saved_fn)
} else {
  # Read count table
  counts <- readRDS(file.path("input", "grossart_lakes", "data.rds"))
  counts <- otu_table(counts)@.Data

  # Filter
  counts <- counts[filter_taxa(counts),]

  cat(paste0("Grossart et al.\n\t", nrow(counts), " ASVs x ", ncol(counts), " samples\n"))

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
  grossart <- fit_model(counts, host_columns, host_dates, "Grossart", depth = 1)[,,1]
  saveRDS(grossart, saved_fn)
}
# plot_heatmap2(grossart, "ASV-ASV pairs", "Grossart")
gr_scores <- apply(grossart, 2, function(x) calc_universality_score(x, return_pieces = TRUE))
gr_scores <- data.frame(X1 = gr_scores[1,], X2 = gr_scores[2,])

# ------------------------------------------------------------------------------
#   McMahon et al.
# ------------------------------------------------------------------------------

saved_fn <- file.path("input", "mcmahon_lakes", "fitted_results.rds")
if(file.exists(saved_fn)) {
  mcmahon <- readRDS(saved_fn)
  mc_e <- mcmahon[[1]]
  mc_h <- mcmahon[[2]]
  rm(mcmahon)
} else {
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

  mc_e <- fit_model(counts, host_columns_E, host_dates_E, "McMahon (epilimnion; shallow water)", depth = 1)[,,1]
  mc_h <- fit_model(counts, host_columns_H, host_dates_H, "McMahon (hypolimnion; deep water)", depth = 1)[,,1]

  # plot_heatmap2(mc_e, "ASV-ASV pairs", "McMahon (epilimnion; shallow water)")
  # plot_heatmap2(mc_h, "ASV-ASV pairs", "McMahon (hypolimnion; deep water)")

  saveRDS(list(mc_e, mc_h), saved_fn)
}
mce_scores <- apply(mc_e, 2, function(x) calc_universality_score(x, return_pieces = TRUE))
mce_scores <- data.frame(X1 = mce_scores[1,], X2 = mce_scores[2,])
mch_scores <- apply(mc_h, 2, function(x) calc_universality_score(x, return_pieces = TRUE))
mch_scores <- data.frame(X1 = mch_scores[1,], X2 = mch_scores[2,])

# ------------------------------------------------------------------------------
#   David et al.
#
#   Note: This data set takes a long time to fit!
# ------------------------------------------------------------------------------

saved_fn <- file.path("input", "david2014", "fitted_results.rds")
if(file.exists(saved_fn)) {
  david <- readRDS(saved_fn)
} else {
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

  david <- fit_model(counts, host_columns, host_dates, "David et al.", depth = 1)[,,1]
  saveRDS(david, saved_fn)
}
# plot_heatmap2(david, "ASV-ASV pairs", "David et al.")
david_scores <- apply(david, 2, function(x) calc_universality_score(x, return_pieces = TRUE))
david_scores <- data.frame(X1 = david_scores[1,], X2 = david_scores[2,])

# ------------------------------------------------------------------------------
#   Caporaso et al.
# ------------------------------------------------------------------------------

saved_fn <- file.path("input", "caporaso2011", "fitted_results.rds")
if(file.exists(saved_fn)) {
  caporaso <- readRDS(saved_fn)
} else {
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
  caporaso <- fit_model(counts[2:nrow(counts),], host_columns, host_dates, "Caporaso et al.", depth = 1)[,,1]
  saveRDS(caporaso, saved_fn)
}
# plot_heatmap2(caporaso, "ASV-ASV pairs", "Caporaso et al.")
cap_scores <- apply(caporaso, 2, function(x) calc_universality_score(x, return_pieces = TRUE))
cap_scores <- data.frame(X1 = cap_scores[1,], X2 = cap_scores[2,])





quit()



# ------------------------------------------------------------------------------
#   Amboseli baboons
# ------------------------------------------------------------------------------

amboseli <- summarize_Sigmas("asv_days90_diet25_scale1")
abrp_scores <- apply(amboseli$rug, 2, function(x) calc_universality_score(x, return_pieces = TRUE))
abrp_scores <- data.frame(X1 = abrp_scores[1,], X2 = abrp_scores[2,])



# TBD: Arrange data into a save-able object for quick parsing later




# ------------------------------------------------------------------------------
#   Plot densities together
# ------------------------------------------------------------------------------

# Density plot from:
# https://stackoverflow.com/questions/23437000/how-to-plot-a-contour-line-showing-where-95-of-values-fall-within-in-r-and-in

# specify desired contour levels
prob <- c(0.9, 0.5)
prob <- c(0.9)

dc1 <- get_density_obj(johnson_scores)
dc2 <- get_density_obj(dr_scores)
dc3 <- get_density_obj(abrp_scores)
dc4 <- get_density_obj(diab_scores)
dc5 <- get_density_obj(gr_scores)
dc6 <- get_density_obj(mce_scores)
dc7 <- get_density_obj(mch_scores)

dc_combined <- rbind(cbind(dc1, type = "J"),
                     cbind(dc2, type = "DR"),
                     cbind(dc3, type = "ABRP"),
                     cbind(dc4, type = "DIAB"),
                     cbind(dc5, type = "Grossart"),
                     cbind(dc6, type = "McMahon-E"),
                     cbind(dc7, type = "McMahon-H"))

dc_combined$type <- factor(dc_combined$type, levels = c("DR", "DIAB", "J", "ABRP", "Grossart", "McMahon-E", "McMachon-H"))

p <- ggplot(data = dc_combined,
            aes(x = Var1, y = Var2, z = prob, fill = type, alpha = ..level..)) +
  geom_contour_filled(breaks = c(prob, 0)) +
  scale_alpha_discrete(range = c(0.5, 0.6)) +
  theme_bw() +
  xlim(c(0.45, 1.05)) +
  ylim(c(0, 1)) +
  labs(x = "proportion shared sign",
       y = "median association strength",
       fill = "Data set") +
  guides(alpha = "none")
print(p)

# Stacked densities
scores_combined <- rbind(data.frame(x = johnson_scores$X2, type = "J"),
                         data.frame(x = dr_scores$X2, type = "DR"),
                         data.frame(x = abrp_scores$X2, type = "ABRP"),
                         data.frame(x = diab_scores$X2, type = "DIAB"),
                         data.frame(x = gr_scores$X2, type = "Grossart"),
                         data.frame(x = mce_scores$X2, type = "McMahon (shallow)"),
                         data.frame(x = mch_scores$X2, type = "McMahon (deep)"))

ggplot(scores_combined, aes(x = x, y = type, fill = type)) +
  geom_density_ridges(stat = "binline", bins = 20, scale = 0.9, draw_baseline = TRUE) +
  theme_bw()












