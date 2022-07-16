#' Description TBD
#'
#' @param seq_similarity_threshold percent upper limit on sequence similarity;
#' for ASV pairs strictly more similar than this, the less abundant partner in
#' the pair will be removed from the data as it's likely these partners ID the
#' same taxon
#' @param host_sample_min minimum sample number for host inclusion in the
#' filtered data set
#' @param count_threshold minimum count for taxon inclusion in the filtered data
#' set
#' @param sample_threshold minimum proportion of samples within each host at
#' which a taxon must be observed at or above count_threshold
#' @details Typically, ASV pairs with >99% similarity differ by 1-2 bases only
#' @return phyloseq object
#' @import driver
#' @import phyloseq
#' @importFrom phyloseq sample_data
#' @import dplyr
#' @import ape
#' @import magrittr
#' @importFrom stringr str_split
#' @export
prune_data <- function(seq_similarity_threshold = 0.99) {
  filename <- file.path("input", "ps0.rds")
  if(file.exists(filename)) {
    data <- readRDS(filename)
  } else {
    stop(paste0("Input data file ", filename, " does not exist!\n"))
  }

  tax <- tax_table(data)

  # 5) Check for extremely closely related sequences (>99% sequence identity)
  # between pairs of taxa. We'll remove the less abundant partner in these pairs
  # as it's possible they represent hits on the same ASVs.

  OTUs <- rownames(tax)
  OTUs <- OTUs[1:length(OTUs)]
  OTUs <- unname(sapply(OTUs, tolower))

  # Split these up into a list of character vectors; that's apparently the input
  # format dist.dna() wants
  OTU_list <- list()
  for(i in 1:length(OTUs)) {
    OTU_list[[i]] <- str_split(OTUs[i], "")[[1]][1:252]
  }
  OTU_list <- as.DNAbin(OTU_list)
  d <- dist.dna(OTU_list, model = "raw") # takes ~1 min.
  d <- 1 - as.matrix(d) # similarity

  closely_related <- as.data.frame(which(d > seq_similarity_threshold,
                                         arr.ind = T)) %>%
    filter(row != col) %>%
    arrange(row, col)
  if(nrow(closely_related) > 0) {
    rownames(closely_related) <- NULL
    # Get unique pairs
    for(i in 1:nrow(closely_related)) {
      if(closely_related[i,1] > closely_related[i,2]) {
        x <- closely_related[i,1]
        closely_related[i,1] <- closely_related[i,2]
        closely_related[i,2] <- x
      }
    }
    closely_related %<>%
      distinct()

    # Get average abundance of each partner in these pairs; we'll axe the less
    # abundant partner
    counts <- otu_table(data)
    mean_abundance <- apply(counts, 2, mean)
    closely_related %<>%
      mutate(remove = case_when(
        mean_abundance[row] < mean_abundance[col] ~ row,
        TRUE ~ col
      ))

    # Prune taxa
    kept_taxa <- rep(TRUE, nrow(tax))
    kept_taxa[closely_related$remove] <- FALSE
    pruned_data <- prune_taxa(kept_taxa, data)
  } else {
    pruned_data <- data
  }

  filename <- file.path("input", paste0("pruned_",
                                        round(seq_similarity_threshold*100),
                                        ".rds"))

  saveRDS(pruned_data, file = filename)
  return(pruned_data)
}

#' This function loads the raw ABRP data (as a phyloseq object) and 1) filters
#' to a set of best-sampled hosts, 2) agglomerates taxa to the level specified,
#' and 3) merges taxa below a specified minimum abundance across hosts
#'
#' @param tax_level taxonomic level at which to agglomerate data
#' @param host_sample_min minimum sample number for host inclusion in the
#' filtered data set
#' @param count_threshold minimum count for taxon inclusion in the filtered data
#' set
#' @param sample_threshold minimum proportion of samples within each host at
#' which a taxon must be observed at or above count_threshold
#' @details Together count_threshold and sample_threshold specify a minimum
#' representation for a taxon. Taxa below this threshold will be grouped
#' together into an <NA> category.
#' @return phyloseq object
#' @import driver
#' @import phyloseq
#' @importFrom phyloseq sample_data
#' @import dplyr
#' @import ape
#' @import magrittr
#' @importFrom stringr str_split
#' @export
filter_data <- function(tax_level = "ASV", host_sample_min = 75,
                      count_threshold = 1, sample_threshold = 0.2,
                      seq_similarity_threshold = 0.99) {
  if(is.null(tax_level)) {
    tax_level <- "ASV"
  }
  if(!(tax_level %in% c("domain", "phylum", "class", "order", "family", "genus",
                        "ASV"))) {
    stop("Unrecognized tax_level specified!")
  }

  # Load pruned data
  filename <- file.path("input", paste0("pruned_",
                                        round(seq_similarity_threshold*100),
                                        ".rds"))
  if(file.exists(filename)) {
    data <- readRDS(filename)
  } else {
    data <- prune_data(seq_similarity_threshold)
  }

  # 1) Filter to hosts above a threshold
  metadata <- sample_data(data)
  metadata <- suppressWarnings(metadata %>%
    group_by(sname) %>%
    filter(n() >= host_sample_min))
  filter_sample_ids <<- metadata$sample_id
  data <- subset_samples(data, sample_id %in% filter_sample_ids)

  cat("Subsetting to", nsamples(data), "samples from",
      length(unique(metadata$sname)), "hosts...\n")

  # 2) Agglomerate taxa; from ASV to family this takes ~5 min.
  if(tax_level == "ASV") {
    agglomerated_data <- data
  } else {
    agglomerated_data <- tax_glom(data, taxrank = tax_level, NArm = FALSE)
  }

  cat("Agglomerated from", ntaxa(data), "to", ntaxa(agglomerated_data), "...\n")

  # 3) Filter to taxa above a minimum prevalence; note: this takes 5+ min. at
  # the ASV level!
  counts <- otu_table(agglomerated_data)@.Data
  binary_counts <- apply(counts, c(1,2), function(x) {
    ifelse(x >= count_threshold, 1, 0)
  })
  # colnames(binary_counts) <- NULL

  sname_map <- metadata[metadata$sample_id %in% rownames(counts),]$sname
  # These are snames (baboon host short names) in order of sample IDs (rows) in
  # the count table

  taxa_passed <- rep(TRUE, ntaxa(agglomerated_data))
  for(sname in unique(sname_map)) {
    host_bin_counts <- binary_counts[sname_map == sname,,drop=F]
    host_sample_sums <- colSums(host_bin_counts) / nrow(host_bin_counts)
    taxa_passed <- taxa_passed & host_sample_sums >= sample_threshold
  }

  # 4) It doesn't look like there are instances of family Mitochondria or order
  # Chloroplast in this data set. Remove archaea though.
  tax <- tax_table(agglomerated_data)@.Data
  archaea_rowidx <- which(tax[taxa_passed,] == "Archaea", arr.ind = TRUE)
  if(length(archaea_rowidx) > 0) {
    archaea_rowidx <- unname(archaea_rowidx[1,1])
    taxa_passed[archaea_rowidx] <- FALSE
  }

  # The `merge_taxa` function in phyloseq will dump all the "other" taxa into
  # a new taxon in a new final index.

  merge_idx <- unname(which(!taxa_passed)[1])
  merged_data <- merge_taxa(agglomerated_data, which(!taxa_passed), archetype = 1)

  # Tag the merged taxons as "other" so we can identify them later
  tax <- tax_table(merged_data)@.Data
  merge_OTU <- rownames(tax)[merge_idx]

  dim(otu_table(merged_data))

  save_obj <- list(merged_data = pruned_data, merge_OTU = merge_OTU)

  filename <- file.path("input", paste0("filtered_",
                                        tax_level,
                                        "_",
                                        count_threshold,
                                        "_",
                                        round(sample_threshold*100),
                                        ".rds"))

  saveRDS(save_obj, file = filename)
  return(save_obj)
}

#' This function wrangles the filtered ABRP data and metadata
#'
#' @param tax_level taxonomic level at which to agglomerate data
#' @param host_sample_min minimum sample number for host inclusion in the
#' filtered data set
#' @param count_threshold minimum count for taxon inclusion in the filtered data
#' set
#' @param sample_threshold minimum proportion of samples within each host at
#' which a taxon must be observed at or above count_threshold
#' @param alr_median if TRUE, uses the taxon with the median coefficient of
#' variation in relative abundance as the ALR reference taxon, otherwise uses
#' the taxon in the last index in the count table
#' @details Together count_threshold and sample_threshold specify a minimum
#' representation for a taxon. Taxa below this threshold will be grouped
#' together into an <NA> category.
#' @return a named list of count table, taxonomy, and metadata components
#' @import phyloseq
#' @importFrom phyloseq psmelt
#' @import dplyr
#' @importFrom tidyr pivot_wider pivot_longer
#' @export
load_data <- function(tax_level = "ASV", host_sample_min = 75,
                      count_threshold = 1, sample_threshold = 0.2,
                      alr_median = FALSE) {
  processed_filename <- file.path("input", paste0("processed_",
                                                  tax_level,
                                                  "_",
                                                  count_threshold,
                                                  "_",
                                                  round(sample_threshold*100),
                                                  ".rds"))
  if(file.exists(processed_filename)) {
    return(readRDS(processed_filename))
  }
  # File not found; process data afresh

  filtered_filename <- file.path("input", paste0("filtered_",
                                        tax_level,
                                        "_",
                                        count_threshold,
                                        "_",
                                        round(sample_threshold*100),
                                        ".rds"))
  if(file.exists(filtered_filename)) {
    data <- readRDS(filtered_filename)
  } else {
    cat("Filtered data file not found. Filtering data now...\n")
    data <- filter_data(tax_level = tax_level,
                        host_sample_min = host_sample_min,
                        count_threshold = count_threshold,
                        sample_threshold = sample_threshold)
  }

  cat("Wrangling data and metadata...\n")
  # Pull taxonomy -> data.frame
  long_data <- psmelt(data$merged_data)

  ordered_data <- long_data %>%
    arrange(sname, collection_date) %>%
    select(c("OTU",
           "Sample",
           "Abundance",
           "plate",
           "extract_dna_conc_ng",
           "sample_status",
           "sname",
           "matgrp",
           "grp",
           "sex",
           "age",
           "collection_date",
           "season",
           "sample_id"))

  ordered_data <- pivot_wider(ordered_data,
                              names_from = "OTU",
                              values_from = "Abundance")

  # Pull off metadata
  metadata <- as.data.frame(ordered_data[,1:12])

  # Pull separate count table (taxa x samples)
  counts <- as.data.frame(ordered_data[,13:ncol(ordered_data)])
  tax_sequences <- colnames(counts)
  colnames(counts) <- NULL
  rownames(counts) <- NULL
  counts <- t(counts)

  # Pull matching taxonomy
  all_tax <- c("domain",
               "phylum",
               "class",
               "order",
               "family",
               "genus")
  if(tax_level == "ASV") {
    use_tax <- all_tax[1:length(all_tax)]
  } else {
    use_tax <- all_tax[1:which(all_tax == tax_level)]
  }

  tax <- as.data.frame(long_data %>%
                         select(c("OTU", all_of(use_tax))) %>%
                         group_by(OTU) %>%
                         distinct())
  tax_sequences <- data.frame(OTU = tax_sequences)
  tax <- left_join(tax_sequences, tax, by = "OTU")

  # Update "other" taxon as other
  merge_idx <- which(tax$OTU == data$merge_OTU)
  tax$OTU[merge_idx] <- "other"

  if(alr_median) {
    # Reshuffle stablest relative abundance (in terms of median coefficient of
    # variation of relative abundances) to the end
    rel_ab <- apply(counts, 2, function(x) x/sum(x))
    c_of_v <- apply(rel_ab, 1, function(x) sd(x)/mean(x))
    median_taxon <- order(c_of_v)[round(length(c_of_v)/2)]
    new_order <- c(setdiff(1:nrow(counts), median_taxon), median_taxon)
  } else {
    # Use the "other" bucket as the reference (last index in the count table)
    new_order <- c(setdiff(1:nrow(counts), merge_idx), merge_idx)
  }

  counts <- counts[new_order,]
  tax <- tax[new_order,]

  # Finally, pull in additional metadata if available
  covariate_filename <- file.path("input", "ps_w_covs.RDS")
  if(file.exists(covariate_filename)) {
    covs <- readRDS(covariate_filename)
    cov_md <- sample_data(covs) %>%
      data.frame() %>%
      select(sample_id, diet_PC1, diet_PC2, diet_PC3, diet_PC4)
    metadata <- left_join(metadata, cov_md, by = "sample_id")
  } else {
    cat("Additional covariates not found...\n")
  }

  processed_data <- list(counts = counts, taxonomy = tax, metadata = metadata)
  saveRDS(processed_data, file = processed_filename)
  return(processed_data)
}

#' Load scrambled data
#'
#' @param tax_level taxonomic level at which to agglomerate data
#' @param host_sample_min minimum sample number for host inclusion in the
#' filtered data set
#' @param count_threshold minimum count for taxon inclusion in the filtered data
#' set
#' @param sample_threshold minimum proportion of samples within each host at
#' which a taxon must be observed at or above count_threshold
#' @return a named list of count table, taxonomy, and metadata components
#' @export
load_scrambled_data <- function(tax_level = "ASV", host_sample_min = 75,
                                    count_threshold = 1, sample_threshold = 0.2) {
  scrambled_filename <- file.path("input", paste0("scrambled_",
                                                  tax_level,
                                                  "_",
                                                  count_threshold,
                                                  "_",
                                                  round(sample_threshold*100),
                                                  ".rds"))
  if(file.exists(scrambled_filename)) {
    data <- readRDS(scrambled_filename)
  } else {
    data <- generate_scrambled_data(tax_level = tax_level,
                                    host_sample_min = host_sample_min,
                                    count_threshold = count_threshold,
                                    sample_threshold = sample_threshold)
  }
  return(data)
}

#' Permute the data within samples, maintaining relative abundances but
#' scrambling patterns of variation within taxa
#'
#' @param tax_level taxonomic level at which to agglomerate data
#' @param host_sample_min minimum sample number for host inclusion in the
#' filtered data set
#' @param count_threshold minimum count for taxon inclusion in the filtered data
#' set
#' @param sample_threshold minimum proportion of samples within each host at
#' which a taxon must be observed at or above count_threshold
#' @return a named list of count table, taxonomy, and metadata components
#' @export
generate_scrambled_data <- function(tax_level = "ASV", host_sample_min = 75,
                                    count_threshold = 1, sample_threshold = 0.2) {
  processed_filename <- file.path("input", paste0("processed_",
                                                  tax_level,
                                                  "_",
                                                  count_threshold,
                                                  "_",
                                                  round(sample_threshold*100),
                                                  ".rds"))
  if(!file.exists(processed_filename)) {
    stop("Processed data file does not exist!")
  }
  data <- readRDS(processed_filename)
  new_counts <- data$counts
  for(j in 1:ncol(new_counts)) {
    new_order <- sample(1:nrow(new_counts))
    new_counts[,j] <- new_counts[new_order,j]
  }
  data$counts <- new_counts
  data$taxonomy <- NULL
  scrambled_filename <- file.path("input", paste0("scrambled_",
                                                  tax_level,
                                                  "_",
                                                  count_threshold,
                                                  "_",
                                                  round(sample_threshold*100),
                                                  ".rds"))
  saveRDS(data, file = scrambled_filename)
  return(data)
}

#' Pull a fitted model object from the specified directory
#'
#' @param sname host short name indicating which baboon's series to fit
#' @param MAP flag indicating whether or not to return (single) MAP estimate
#' @param output_dir specifies the subdirectory of the model fits directory in
#' which to save the predictions
#' @return a bassetfit object or NULL if no such object exists
#' @export
load_fit <- function(sname, MAP, output_dir) {
  if(MAP) {
    fit_dir <- "MAP"
  } else {
    fit_dir <- "full_posterior"
  }
  filename <- file.path("output", "model_fits", output_dir, fit_dir, paste0(sname, ".rds"))
  if(file.exists(filename)) {
    return(readRDS(filename))
  } else {
    return(NULL)
  }
}

#' Pull a prediction object from the specified directory
#'
#' @param sname host short name indicating which baboon's series to fit
#' @param output_dir specifies the subdirectory of the model fits directory in
#' which to save the predictions
#' @param generate flag indicating whether or not to generate prediction data
#' if not found
#' @return a named list of predicted values or NULL if no such object exists
#' @export
load_predictions <- function(sname, output_dir, generate = FALSE) {
  filename <- file.path("output", "model_fits", output_dir, "predictions", paste0(sname, ".rds"))
  if(file.exists(filename)) {
    return(readRDS(filename))
  } else {
    if(generate) {
      fit_pred_obj <- predict_trajectory(sname, output_dir)
      return(fit_pred_obj$predictions)
    } else {
      return(NULL)
    }
  }
}

#' Assign host social group as social group in which a given host has the
#' majority of its samples
#'
#' @param host_list optional list of host short names
#' @return list of host short name to social group mappings
#' @import dplyr
#' @export
get_host_social_groups <- function(host_list = NULL) {
  metadata <- load_data()$metadata
  results <- as.data.frame(metadata %>%
                             select(sname, grp) %>%
                             count(sname, grp) %>%
                             group_by(sname) %>%
                             slice_max(order_by = n, n = 1) %>%
                             arrange(sname) %>%
                             select(sname, grp))
  if(!is.null(host_list)) {
    results <- results %>%
      filter(sname %in% host_list)
  }
  # grp_assignments <- results$grp
  # names(grp_assignments) <- results$sname
  return(results)
}

#' Find host with largest numbers of samples in each social group
#'
#' @param host_sample_min minimum sample number for host inclusion in the
#' filtered data set
#' @param n_per_group number of hosts to return per group
#' @return vector of host short names
#' @import dplyr
#' @export
get_reference_hosts <- function(host_sample_min = 75, n_per_group = 1) {
  metadata <- load_data(host_sample_min = host_sample_min)$metadata
  # For each group, select the host with the largest number of  samples in that
  # group ()
  grp_assignments <- as.data.frame(metadata %>%
                                     select(grp, sname) %>%
                                     group_by(grp) %>%
                                     count(sname) %>%
                                     slice_max(order_by = n, n = n_per_group) %>%
                                     filter(n > host_sample_min) %>%
                                     arrange(grp) %>%
                                     select(grp, sname))
  return(grp_assignments)
}
