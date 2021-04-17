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
#' @import phyloseq
#' @importFrom phyloseq sample_data
#' @import dplyr
#' @export
filter_data <- function(tax_level = "ASV", host_sample_min = 75,
                      count_threshold = 5, sample_threshold = 0.2) {
  if(is.null(tax_level)) {
    tax_level <- "ASV"
  }
  if(!(tax_level %in% c("domain", "phylum", "class", "order", "family", "genus",
                        "species", "ASV"))) {
    stop("Unrecognized tax_level specified!")
  }
  filename <- file.path("input", "ps0.rds")
  if(file.exists(filename)) {
    data <- readRDS(filename)
  } else {
    stop(paste0("Input data file ", filename, " does not exist!\n"))
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

  # 2) Agglomerate taxa; from ASV level to family level this takes ~5 min (FYI)
  if(tax_level == "ASV") {
    agglomerated_data <- data
  } else {
    agglomerated_data <- tax_glom(data, taxrank = tax_level, NArm = FALSE)
  }

  cat("Agglomerated from", ntaxa(data), "to", ntaxa(agglomerated_data), "...\n")

  # 3) Filter to taxa above a minimum prevalence
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
    host_bin_counts <- binary_counts[sname_map == sname,]
    host_sample_sums <- colSums(host_bin_counts) / nrow(host_bin_counts)
    taxa_passed <- taxa_passed & host_sample_sums >= sample_threshold
  }

  # It doesn't look like there are instances of family Mitochondria or order
  # Chloroplast in this data set. Remove archaea though.

  # The `merge_taxa` function in phyloseq will dump all the "other" taxa into
  # a new taxon in a new final index.

  merge_idx <- which(!taxa_passed)[1]
  merged_data <- merge_taxa(agglomerated_data, which(!taxa_passed))

  filename <- file.path("input", paste0("filtered_",
                                        tax_level,
                                        "_",
                                        count_threshold,
                                        "_",
                                        round(sample_threshold*100),
                                        ".rds"))
  saveRDS(merged_data, file = filename)
  return(merged_data)
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
#' @details Together count_threshold and sample_threshold specify a minimum
#' representation for a taxon. Taxa below this threshold will be grouped
#' together into an <NA> category.
#' @return A named list of count table, taxonomy, and metadata components
#' @import phyloseq
#' @importFrom phyloseq psmelt
#' @import dplyr
#' @importFrom tidyr pivot_wider pivot_longer
#' @export
load_data <- function(tax_level = "ASV", host_sample_min = 75,
                        count_threshold = 5, sample_threshold = 0.2) {
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
  long_data <- psmelt(data)

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
               "genus",
               "species")
  use_tax <- all_tax[1:which(all_tax == tax_level)]

  tax <- as.data.frame(long_data %>%
                         select(c("OTU", all_of(use_tax))) %>%
                         group_by(OTU) %>%
                         distinct())
  tax_sequences <- data.frame(OTU = tax_sequences)
  tax <- left_join(tax_sequences, tax, by = "OTU")

  # Reshuffle stablest relative abundance (in terms of median coefficient of
  # variation of relative abundances) to the end

  rel_ab <- apply(counts, 2, function(x) x/sum(x))
  c_of_v <- apply(rel_ab, 1, function(x) sd(x)/mean(x))
  median_taxon <- order(c_of_v)[round(length(c_of_v)/2)]

  new_order <- c(setdiff(1:nrow(counts), median_taxon), median_taxon)
  counts <- counts[new_order,]
  tax <- tax[new_order,]

  processed_data <- list(counts = counts, taxonomy = tax, metadata = metadata)
  saveRDS(processed_data, file = processed_filename)
  return(processed_data)
}
