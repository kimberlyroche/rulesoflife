source("path_fix.R")

library(tidyverse)
library(rulesoflife)

# Here first - collapse these to family level so stringency of filtering at
# least sort of makes sense
# Should these be collapsed to comparable taxonomic levels -- family?

output_dir <- "fam_days90_diet25_scale1"

if(file.exists(file.path("output", "rug_fam.rds"))) {
  rug_obj <- readRDS(file.path("output", "rug_fam.rds"))
} else {
  rug_obj <- summarize_Sigmas(output_dir)
}

scores_ABRP <- apply(rug_obj$rug, 2, calc_universality_score)

scores_df <- data.frame(scores = scores_ABRP,
                        dataset = "ABRP",
                        tax_level = "ASV")

ggplot(scores_df, aes(x = dataset, y = scores, fill = tax_level)) +
  geom_boxplot()

# ------------------------------------------------------------------------------
#   Johnson et al. 2019
# ------------------------------------------------------------------------------

# Parse read count table
counts <- read.table(file.path("input", "johnson2019", "taxonomy_counts_s.txt"),
                     header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# Parse read sample ID to subject ID mapping
mapping <- read.table(file.path("input", "johnson2019", "SampleID_map.txt"),
                      header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Chop out taxonomy (first column)
counts <- counts[,2:ncol(counts)]


# Left off here
# Need to evaluate similarity of filtering between ABRP and Johnson et al. data
# set


# Filter to some minimum relative average abundance
retain_idx <- filter_taxa_2D_summary(counts)
counts <- counts[retain_idx,]

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

# # Testing
# for(i in 1:length(host_columns)) {
#   limit <- min(5, length(host_columns[[i]]))
#   host_columns[[i]] <- host_columns[[i]][1:limit]
#   host_dates[[i]] <- host_dates[[i]][1:limit]
# }

vectorized_Sigmas <- fit_or_visualize_sub(counts, host_columns, host_dates, visualize, data_file, GP)
