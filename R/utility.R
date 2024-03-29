#' Create output directory/ies that may not already exist
#'
#' @param path_to_dir a character vector of subdirectories to create (if not
#' already in existence)
#' @return the path as a joined string
#' @export
check_dir <- function(path_to_dir) {
  for(i in 1:length(path_to_dir)) {
    dir.create(do.call(file.path, as.list(path_to_dir[1:i])),
               showWarnings = FALSE)
  }
  return(do.call(file.path, as.list(path_to_dir)))
}

#' Return a human-readable taxon label
#'
#' @param tax taxonomy data.frame
#' @param coord index of logratio to label
#' @param coord_system_label one of c("CLR", "ALR")
#' @return label
#' @export
get_tax_label <- function(tax, coord, coord_system_label) {
  coord_system_label <- toupper(coord_system_label)
  if(coord_system_label %in% c("CLR", "ALR")) {
    if(coord > nrow(tax)) {
      stop("Invalid taxonomy index")
    }
    tax_pieces <- tax[coord,]
    max_level <- max(which(!is.na(tax_pieces)))
    if(coord_system_label == "CLR") {
      label <- paste0(coord_system_label, "(ASVs in ", names(tax)[max_level], " ", tax_pieces[max_level], ")")
    } else {
      max_level_ref <- max(which(!is.na(tax[nrow(tax),])))
      # ALR
      numerator <- paste0("ASVs in ", names(tax)[max_level], " ", tax_pieces[max_level])
      denominator <- paste0("ASVs in ", names(tax)[max_level_ref], " ", tax[nrow(tax),max_level])
      label <- paste0("ALR(", numerator, "/", denominator, ")")
    }
    return(label)
  } else {
    stop("Unrecognized coordinate system!")
  }
}

#' Calculate a score (0.0 to 1.0) for the "universality" of given association
#' pair
#'
#' @param x vector of correlations between a pair of CLR microbes across hosts
#' @param return_pieces if TRUE, returns the two components of the universality
#' score: the proportion agreement in direction and the median magnitude of
#' correlation in that set
#' @return numeric "score"
#' @export
calc_universality_score <- function(x, return_pieces = FALSE) {
  x.sign <- sapply(x, sign)
  neg.idx <- which(x.sign < 0)
  pos.idx <- which(x.sign > 0)
  neg.abs.median <- NA
  pos.abs.median <- NA
  if(length(neg.idx) > 0) {
    neg.abs.median <- abs(median(x[neg.idx]))
  }
  if(length(pos.idx) > 0) {
    pos.abs.median <- abs(median(x[pos.idx]))
  }
  score <- 0
  # R doesn't do short circuit evaluation? yes, it does: update with &&
  if(!is.na(neg.abs.median) & !is.na(pos.abs.median)) {
    if(length(pos.idx) > length(neg.idx)) {
      # use this as majority direction
      p1 <- length(pos.idx)/length(x)
      p2 <- pos.abs.median
    } else {
      p1 <- length(neg.idx)/length(x)
      p2 <- neg.abs.median
    }
  } else {
    if(is.na(neg.abs.median)) {
      p1 <- length(pos.idx)/length(x)
      p2 <- pos.abs.median
    } else {
      p1 <- length(neg.idx)/length(x)
      p2 <- neg.abs.median
    }
  }
  if(return_pieces) {
    return(c(p1, p2))
  } else {
    return (p1*p2)
  }
}

#' Calculate a score (0.0 to 1.0) for the "universality" of association pairs
#' for each taxon
#'
#' @param rug_obj output of the summarize_Sigmas function
#' @param collapse if value is "signed", computes an average
#' universality score for a taxon split up for consensus positive and negative
#' associations; if "unsigned" (default), computes na overall average universality score
#' for a taxon
#' @return list of numeric "scores"
#' @export
calc_universality_score_taxon <- function(rug_obj, collapse = "unsigned") {
  scores <- apply(rug_obj$rug, 2, calc_universality_score)
  signs <- apply(rug_obj$rug, 2, calc_consensus_sign)
  df <- data.frame(score = scores,
                   sign = signs,
                   tax1 = rug_obj$tax_idx1,
                   tax2 = rug_obj$tax_idx2)
  df <- df %>%
    filter(sign != 0)
  if(collapse == "signed") {
    df %>%
      group_by(tax1, sign) %>%
      mutate(mean_score = mean(score)) %>%
      select(tax1, sign, mean_score) %>%
      distinct()
  } else {
    df %>%
      group_by(tax1) %>%
      mutate(mean_score = mean(score)) %>%
      select(tax1, mean_score) %>%
      distinct()
  }
}

#' Calculate a consensus CLR correlation sign for a given association pair
#'
#' @param x vector of correlations between a pair of CLR microbes across hosts
#' @return positive or negative consensus sign
#' @export
calc_consensus_sign <- function(x) {
  sign(sum(sign(x)))
}

#' Reorder columns of the "rug" from taxa with the least to most difference in
#' mean CLR abundance
#'
#' @param rug_obj output of the summarize_Sigmas function
#' @param counts count table at the same taxonomic level as results in rug_obj
#' @return named list of column ordering and absolute CLR differences
#' @import driver
#' @export
order_rug_cols_mean_abundance <- function(rug_obj, counts) {
  clr.means <- rowMeans(clr_array(counts + 0.5, parts = 1))
  clr.diff <- sapply(1:length(rug_obj$tax_idx1), function(i) {
    abs(clr.means[rug_obj$tax_idx1[i]] - clr.means[rug_obj$tax_idx2[i]])
  })
  new_order <- order(clr.diff)
  return(list(order = new_order, difference = clr.diff))
}

#' Reorder rows of the "rug" based on Aitchison distance between host average
#' compositions
#'
#' @param rug_obj output of the summarize_Sigmas function
#' @param counts count table at the same taxonomic level as results in rug_obj
#' @param metadata metadata at the same taxonomic level as results in rug_obj
#' @return row ordering
#' @import driver
#' @export
order_rug_row_baseline <- function(rug_obj, counts, metadata) {
  clr.counts <- clr_array(counts + 0.5, parts = 1)
  # Get each host's "baseline" composition (avg. CLR composition)
  hosts <- rug_obj$hosts
  host_baselines <- matrix(NA, length(hosts), nrow(clr.counts))
  for(i in 1:length(hosts)) {
    host <- hosts[i]
    host_data <- clr.counts[, metadata$sname == host]
    host_baselines[i,] <- rowMeans(host_data)
  }
  host_dist <- dist(host_baselines)
  hc <- hclust(host_dist)
  return(list(order = hc$order,
              hc = hc))
}

#' Reorder rows of the "rug" based on average taxonomic average host alpha
#' diversity (low to high Shannon index)
#'
#' @param rug_obj output of the summarize_Sigmas function
#' @return row ordering
#' @import phyloseq
#' @import vegan
#' @export
order_rug_row_diversity <- function(rug_obj) {
  # Load raw data; we'll calculate alpha-diversity from this rather than the
  # processed data, which consolidates to taxa-in-common across hosts
  raw_data <- readRDS(file.path("input", "ps0.rds"))
  # Subset to hosts in question and calculate average Shannon index
  host_diversity <- c()
  hosts <- rug_obj$hosts
  for(i in 1:length(hosts)) {
    host <<- hosts[i]
    host_samples <- subset_samples(raw_data, sname == host)
    host_counts <- otu_table(host_samples)@.Data
    host_diversity[i] <- mean(unname(apply(host_counts, 1, diversity)))
  }
  host_dist <- dist(host_diversity)
  hc <- hclust(host_dist)
  return(list(order = hc$order,
              hc = hc))
}

#' Reorder rows of the "rug" based on primary social group
#'
#' @param rug_obj output of the summarize_Sigmas function
#' @return row ordering
#' @import dplyr
#' @export
order_rug_row_group <- function(rug_obj) {
  grp_assignments <- get_host_social_groups(rug_obj$hosts)
  grp_assignments <- grp_assignments %>%
    arrange(grp)
  temp <- data.frame(sname = rug_obj$hosts)
  temp <- temp %>%
    left_join(grp_assignments, temp, by = "sname") %>%
    mutate(new_label = paste0(sname, " (", grp, ")"))
  labels <- temp %>%
    pull(new_label)
  temp <- temp %>%
    arrange(grp)
  map <- data.frame(sname = rug_obj$hosts,
                    index = 1:length(rug_obj$hosts))
  index <- temp %>%
    right_join(map, by = "sname") %>%
    pull(index)
  return(list(order = index,
              label = labels))
}

#' Reorder rows of the "rug" based on known pedigree
#'
#' @param rug_obj output of the summarize_Sigmas function
#' @return row ordering
#' @export
order_rug_row_pedigree <- function(rug_obj) {
  pedigree <- readRDS(file.path("input", "pedigree_56hosts.RDS"))

  # The pedigree has an arbitrary(?) ordering of hosts. We need to make sure this
  # order matches our current indexing in the rug (alphabetical).
  host_order_rug <- data.frame(host = rug_obj$hosts,
                               index_rug = 1:length(rug_obj$hosts))
  host_order_ped <- data.frame(host = rownames(pedigree),
                               index_ped = 1:nrow(pedigree))
  host_reordering <- left_join(host_order_rug, host_order_ped, by = "host")$index_ped

  # `host_reodering` is the mapping of indices from the pedigree back to
  # alphabetical. Perform this ordering of the pedigree matrix.
  pedigree2 <- pedigree[host_reordering,host_reordering]

  # Use 1 - % genes shared as the distance to cluster rows on.
  hc <- hclust(as.dist(1-pedigree2))
  return(list(order = hc$order,
              hc = hc))
}

#' Calculate phylogenetic distances on pairs of taxa
#'
#' @param rug_obj output of the summarize_Sigmas function
#' @param taxonomy taxonomy appropriate for data at the same level as the rug
#' object
#' @param as_matrix if TRUE, returns the distance matrix, not the vectorized
#' result
#' @return row ordering
#' @import Biostrings
#' @import phangorn
#' @import DECIPHER
#' @export
rug_phylogenetic_distances <- function(rug_obj, taxonomy, as_matrix = FALSE) {
  seqs <- taxonomy$OTU
  seqs <- seqs[1:(length(seqs)-1)]
  alignment <- AlignSeqs(DNAStringSet(seqs), anchor = NA)
  phang.align <- phyDat(as(alignment, "matrix"), type = "DNA")
  dm <- dist.ml(phang.align)
  dm <- as.matrix(dm)
  if(as_matrix) {
    return(dm)
  }
  tax_distances <- sapply(1:length(rug_obj$tax_idx1), function(i) {
    dm[rug_obj$tax_idx1[i], rug_obj$tax_idx2[i]]
  })
  return(tax_distances)
}

#' Reorder columns of the "rug" based on phylogenetic distance
#'
#' @param rug_obj output of the summarize_Sigmas function
#' @param taxonomy taxonomy appropriate for data at the same level as the rug
#' object
#' @return row ordering
#' @import Biostrings
#' @import phangorn
#' @import DECIPHER
#' @export
order_rug_col_phylogenetic <- function(rug_obj, taxonomy) {
  tax_distances <- rug_phylogenetic_distances(rug_obj, taxonomy)
  return(list(order = order(tax_distances)))
}

#' Get the mean CLR abundance for all taxa
#'
#' @param tax_level taxonomic level at which to calculate average CLR abundance
#' @return vector of mean CLR abundances
#' @import driver
#' @export
get_mean_clr_abundance <- function(tax_level = "ASV") {
  data <- load_data(tax_level = tax_level)
  clr.counts <- clr_array(data$counts + 0.5, parts = 1)
  rowMeans(clr.counts)
}

#' Get the mean CLR abundance for all taxa
#'
#' @param tax_level taxonomic level at which to calculate average CLR abundance
#' @return vector of mean CLR abundances
#' @import driver
#' @export
get_mean_clr_abundance <- function(tax_level = "ASV") {
  data <- load_data(tax_level = tax_level)
  clr.counts <- clr_array(data$counts + 0.5, parts = 1)
  rowMeans(clr.counts)
}

#' Converts a percent to a relative count
#'
#' @param percent percent
#' @param n_features total number of features
#' @return relative count
#' @export
percent_to_k <- function(percent, n_features) {
  round(n_features*(percent/100))
}

#' Predict mean trajectories (Lambda) for all CLR taxa in a given host
#'
#' @param fit an optional fitted model object
#' @param host an optional host short name
#' @param output_dir an optional model fit directory
#' @param interpolation allowed values are "mean", "linear", or "none"
#' @return named list with inferred trajectories and time span
#' @details Either `fit` or `host` and `output_dir` must be provided. If all are
#' provided, `fit` overrides all. Note also: the runtime on this function is
#' around 10s per host at the ASV-level.
#' @export
predict_GP_mean <- function(fit = NULL, host = NULL, output_dir = NULL,
                            interpolation = "linear") {
  if(!(interpolation %in% c("mean", "linear", "none"))) {
    stop(paste0("Disallowed interpolation method: ", interpolation, "\n"))
  }
  if(is.null(fit)) {
    if(is.null(host)) {
      stop("Missing host short name!")
    }
    if(is.null(output_dir)) {
      stop("Missing model output directory!")
    }
    fit_filename <- file.path("output", "model_fits", output_dir,
                              "full_posterior", paste0(host, ".rds"))
    if(!file.exists(fit_filename)) {
      stop(paste0("No such file exists: ", fit_filename, "\n"))
    }
    fit <- readRDS(fit_filename)
  }

  # Observed data
  X_o <- fit$X

  # Build unobserved data set
  first_day <- min(X_o[1,])
  last_day <- max(X_o[1,])
  n_days <- last_day
  span <- seq(from = first_day, to = last_day)
  X_u <- matrix(NA, nrow(X_o), length(span))
  X_u[1,] <- span
  # Linearly interpolate the diet PCs for lack of a better option
  x <- fit$X[1,]
  for(cov_idx in 2:4) {
    y <- fit$X[cov_idx,]
    if(interpolation == "linear") {
      X_u[cov_idx,] <- approx(x = x, y = y, xout = span, ties = "ordered")$y
    } else if(interpolation == "mean") {
      X_u[cov_idx,] <- mean(y)
      X_u[cov_idx,x] <- y
    } else {
      X_u[cov_idx,] <- 0
    }
  }

  Gamma <- fit$Gamma(cbind(X_o, X_u))
  obs <- c(rep(TRUE, ncol(fit$X)), rep(FALSE, ncol(X_u)))

  # Predict Lambda
  Gamma_oo <- Gamma[obs, obs, drop=F]
  Gamma_ou <- Gamma[obs, !obs, drop=F]

  Theta_o <- fit$Theta(X_o)
  Theta_u <- fit$Theta(X_u)
  Lambda_o <- apply(fit$Lambda, c(1,2), mean)
  mean_Lambda_pred <- Theta_u + (Lambda_o-Theta_o)%*%solve(Gamma_oo, Gamma_ou)

  # convert to CLR
  mean_Lambda.clr <- clr_array(alrInv_array(mean_Lambda_pred, fit$alr_base, 1), 1)

  return(list(predictions = mean_Lambda.clr, span = span))
}

#' Pull the predicted host x taxon series from a data.frame of predictions
#'
#' @param prediction_obj a long tibble containing aligned centered trajectories
#' @param sname host short name
#' @param idx taxon index
#' @return filtered tibble
#' @import dplyr
#' @export
pull_series <- function(prediction_obj, sname, idx) {
  prediction_obj %>%
    filter(host == sname) %>%
    arrange(day) %>%
    filter(coord == idx) %>%
    pull(centered_clr)
}

#' Sample to construct between- and within-host distributions of correlation
#' for a given pair of taxa. The "within" distribution characterizes correlation
#' of a pair of taxa within selected hosts; the "between" distribution
#' characterizes the correlation of a given taxon with its own series across
#' a samples of hosts.
#'
#' @param pair index of taxon pair
#' @param coord1 index of taxon 1
#' @param coord2 index of taxon 2
#' @param predictions predictions data.frame subsetted to the taxa of interest
#' @param selected_hosts subset of full host short name list to use
#' @param verbose if TRUE, prints status comments
#' @return data.frame of between- and within- correlation distribution samples
#' @export
build_between_within_distributions <- function(pair, coord1, coord2,
                                               predictions, selected_hosts,
                                               verbose = TRUE) {
  if(verbose) {
    cat(paste0("Evaluating pair #", pair, "...\n"))
    cat("Subsetting predictions to coordinates of interest...\n")
    cat("Building within-host distribution...\n")
  }
  within_distro <- c()
  for(host in selected_hosts) {
    series1 <- pull_series(predictions, host, coord1)
    series2 <- pull_series(predictions, host, coord2)
    within_distro <- c(within_distro, cor(series1, series2))
  }

  if(verbose) {
    cat("Building between-host distribution...\n")
  }
  host_combos <- combn(selected_hosts, m = 2)
  if(subset_hosts) {
    host_combos <- host_combos[,sample(1:ncol(host_combos), size = 100)]
  }
  between_distros <- sapply(1:ncol(host_combos), function(i) {
    c(cor(pull_series(predictions, host_combos[1,i], coord1),
          pull_series(predictions, host_combos[2,i], coord1)),
      cor(pull_series(predictions, host_combos[1,i], coord2),
          pull_series(predictions, host_combos[2,i], coord2)))
  })

  temp_df <- data.frame(correlation = within_distro,
                        type = "within hosts")
  temp_df <- bind_rows(temp_df,
                       data.frame(correlation = between_distros[1,],
                                  type = "between hosts (1)"))
  temp_df <- bind_rows(temp_df,
                       data.frame(correlation = between_distros[2,],
                                  type = "between hosts (2)"))
  temp_df$pair <- pair
  temp_df$tax1 <- coord1
  temp_df$tax2 <- coord2
  temp_df
}

#' Convert an ALR covariance matrix to a CLR correlation matrix
#'
#' @param Sigma covariance matrix
#' @return named list of CLR covariance and CLR correlation matrices
#' @export
convert_alr_Sigma_clr <- function(Sigma) {
  # Convert to CLR
  D <- nrow(Sigma) + 1
  Sigma.clr <- alrvar2clrvar(Sigma, d1 = D)
  Sigma.clr.cor <- cov2cor(Sigma.clr)
  return(list(Sigma.clr = Sigma.clr, Sigma.clr.cor = Sigma.clr.cor))
}

#' Compute sequence distances across ASVs
#'
#' @param distance_type optional; default is "N" (number of mismatches)
#' @return matrix of ASV-ASV distances in terms of their 16S sequences
#' @import stringr
#' @import ape
#' @export
sequence_distance <- function(distance_type = "N") {
  data <- load_data(tax_level = "ASV")
  OTUs <- data$taxonomy$OTU
  OTUs <- OTUs[1:(length(OTUs)-1)]
  OTUs <- unname(sapply(OTUs, tolower))

  # Split these up into a list of character vectors; that's apparently the input
  # format dist.dna() wants
  OTU_list <- list()
  for(i in 1:length(OTUs)) {
    OTU_list[[i]] <- str_split(OTUs[i], "")[[1]][1:252]
  }
  OTU_list <- as.DNAbin(OTU_list)

  if(distance_type == "raw") {
    d <- dist.dna(OTU_list, model = "raw")
    d <- 1 - as.matrix(d)
  } else {
    d <- as.matrix(dist.dna(OTU_list, model = distance_type))
  }
  d
}

#' Gather all host ASV-ASV (MAP) correlation matrices into a 3D array
#'
#' @param output_dir specifies the subdirectory of the model fits directory in
#' which to look for fitted model output
#' @return D x D x H array of host ASV-ASV correlation matrices
#' @import fido
#' @export
pull_Sigmas <- function(output_dir) {
  # Get all fitted model objects
  output_dir_full <- check_dir(c("output", "model_fits", output_dir, "MAP"))
  file_list <- list.files(path = output_dir_full, pattern = "*.rds")
  if(length(file_list) == 0) {
    stop("No model output found!")
  }
  # Get taxa number and posterior sample number
  fit <- readRDS(file.path(output_dir_full, file_list[1]))
  D <- fit$D
  Sigmas <- array(NA, dim = c(D-1, D-1, length(file_list)))
  for(f in 1:length(file_list)) {
    file <- file_list[f]
    fit <- readRDS(file.path(output_dir_full, file))
    # Convert to CLR
    if(fit$coord_system != "clr") {
      fit <- to_clr(fit)
    }
    Sigma_correlation <- cov2cor(fit$Sigma[,,1])
    Sigma_correlation <- Sigma_correlation[1:(D-1),1:(D-1)]
    Sigmas[,,f] <- Sigma_correlation
  }
  return(Sigmas)
}

# Internal functions

# Return a D-length logical vector indicating whether a given taxon participates
# in any non-filtered pairs. We'll later use this to excluded (in labeling) any
# taxa whose associations with other taxa are fully filtered out.
#' @import dplyr
#' @import tidyr
#' @export
represented_taxa <- function(filtered_pairs) {
  temp <- filtered_pairs %>%
    select(idx1, idx2, threshold) %>%
    pivot_longer(cols = c(idx1, idx2), names_to = "source") %>%
    group_by(value) %>%
    mutate(passes_threshold = any(threshold)) %>%
    ungroup() %>%
    select(index = value, passes_threshold) %>%
    distinct()

  D <- nrow(temp)
  keep <- logical(D)
  for(i in 1:D) {
    if(temp$passes_threshold[i]) {
      keep[i] <- TRUE
    }
  }
  keep
}

# Reindex remaining taxa after excluding those who do not participate in any
# unfiltered pairs
#' @import dplyr
#' @import tidyr
#' @export
renumber_taxon <- function(representation, taxon_idx) {
  if(representation[taxon_idx]) {
    data.frame(orig_index = 1:length(representation),
               represented = representation) %>%
      group_by(represented) %>%
      mutate(new_index = case_when(
        represented ~ row_number(),
        T ~ NA_integer_
      )) %>%
      ungroup() %>%
      as.data.frame() %>%
      filter(row_number() == taxon_idx) %>%
      pull(new_index)
  } else {
    NA
  }
}

