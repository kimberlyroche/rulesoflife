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
#' @return numeric "score"
#' @export
calc_universality_score <- function(x) {
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
      score <- (length(pos.idx)/length(x))*(pos.abs.median)
    } else {
      score <- (length(neg.idx)/length(x))*(neg.abs.median)
    }
  } else {
    if(is.na(neg.abs.median)) {
      score <- (length(pos.idx)/length(x))*(pos.abs.median)
    } else {
      score <- (length(neg.idx)/length(x))*(neg.abs.median)
    }
  }
  score
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
  return(hclust(host_dist)$order)
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
  return(order(host_diversity))
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
  temp <- data.frame(sname = rug_obj$hosts, index = 1:length(rug_obj$hosts))
  temp <- left_join(grp_assignments, temp, by = "sname")
  return(temp$index)
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
  row_order <- hclust(as.dist(1-pedigree2))$order
  return(row_order)
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
order_rug_row_pedigree <- function(rug_obj, taxonomy) {
  # Calculate distances across sequences associated with the ASVs we've retained,
  # post-filtering.
  seqs <- taxonomy$OTU
  seqs <- seqs[1:(length(seqs)-1)]
  alignment <- AlignSeqs(DNAStringSet(seqs), anchor = NA)
  phang.align <- phyDat(as(alignment, "matrix"), type = "DNA")
  dm <- dist.ml(phang.align)
  dm <- as.matrix(dm)
  # Reorder rug based on these distances
  tax_distances <- sapply(1:length(rug_obj$tax_idx1), function(i) {
    dm[rug_obj$tax_idx1[i], rug_obj$tax_idx2[i]]
  })
  return(order(tax_distances))
}
